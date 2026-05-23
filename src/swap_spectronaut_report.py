"""Gold-standard precursor-intensity swap for Spectronaut reports.

Pairs proteins from the high- and low-abundance tails with matching precursor
counts, then in a randomly chosen "G2" half of each dilution's runs swaps the
intensities (fragment / precursor / peptide / protein levels) rank-for-rank
between paired precursors. Output is schema-identical to the input Spectronaut
TSV minus dropped rows for surplus precursors / peptides / fragments.

See /Users/wolski/.claude/plans/just-to-make-sure-cheeky-sprout.md for the
algorithmic plan this implements.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import polars as pl

# ---------------------------------------------------------------------------
# Column conventions in the Spectronaut report
# ---------------------------------------------------------------------------

COL_RUN = "R.FileName"
COL_COND = "R.Condition"
COL_PG = "PG.ProteinGroups"
COL_PEP = "PEP.StrippedSequence"
COL_PREC = "EG.PrecursorId"

# Intensity columns to swap, grouped by level.
F_COLS = ["F.PeakArea", "F.PeakHeight", "F.NormalizedPeakArea", "F.NormalizedPeakHeight"]
FG_COLS = ["FG.Quantity", "FG.MS2RawQuantity"]
PEP_COLS = ["PEP.Quantity"]
PG_COLS = ["PG.Quantity"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--report", required=True, type=Path)
    p.add_argument("--annotation", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--swap-fraction", type=float, default=0.05,
                   help="Fraction of proteins to make Positive pair members. "
                        "target_n_pairs = round(swap_fraction * n_proteins).")
    p.add_argument("--min-precursors", type=int, default=2,
                   help="Drop proteins with fewer precursors than this from "
                        "the entire eligible universe (Positives and Negatives).")
    p.add_argument("--min-log2fc", type=float, default=0.5,
                   help="Minimum |log2(prot_mean_hi/prot_mean_lo)| a pair must "
                        "exceed; no upper bound. The pair-builder prefers pairs "
                        "closest to this minimum and expands outward as needed.")
    p.add_argument("--good-rule", default=None,
                   choices=("label_good", "neat_only"),
                   help="Restrict the reference run set used by "
                        "compute_protein_stats. label_good: keep annotation "
                        "rows with Label == 'Good'. neat_only: keep R.FileName "
                        "matching 'NeatCSF'. Default (None): use all non-blank "
                        "runs. Mirrors `src/build_subsets.py --good-rule` so "
                        "the pair-level abundance gap is measured on the same "
                        "subset that build_subsets/good_data uses downstream.")
    p.add_argument("--blank-condition", default="blank",
                   help="R.Condition value identifying blank runs (excluded from stats and output groups).")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def load_report(path: Path) -> pl.DataFrame:
    print(f"[load] reading {path}", file=sys.stderr)
    df = pl.read_csv(
        path,
        separator="\t",
        null_values=["NaN", ""],
        infer_schema_length=20000,
        low_memory=False,
    )
    print(f"[load] {df.height:,} rows × {df.width} cols", file=sys.stderr)
    return df


# ---------------------------------------------------------------------------
# Step 1 — per-protein statistics on the reference set
# ---------------------------------------------------------------------------

def filter_runs_by_good_rule(annotation: pl.DataFrame, good_rule: str) -> pl.DataFrame:
    """Filter an annotation to the 'good' runs. Mirrors
    `src/build_subsets.py:filter_good` so the pair-level abundance gap can
    be measured on the same subset the downstream good_data benchmark uses.
    """
    if good_rule == "label_good":
        if "Label" not in annotation.columns:
            raise ValueError(
                "--good-rule=label_good requires a 'Label' column in the annotation."
            )
        return annotation.filter(pl.col("Label") == "Good")
    if good_rule == "neat_only":
        return annotation.filter(pl.col(COL_RUN).str.contains("NeatCSF"))
    raise ValueError(f"Unknown good_rule: {good_rule}")


def compute_protein_stats(df: pl.DataFrame, blank_cond: str,
                            min_precursors: int = 2,
                            reference_runs: list[str] | None = None) -> pl.DataFrame:
    """One row per protein with prot_mean_log2, n_precursors, n_peptides.

    Proteins with `n_precursors < min_precursors` are dropped from the
    eligible universe entirely (so they are neither Positives nor
    Negatives in the downstream swap_list).

    If `reference_runs` is given (a list of R.FileName values), the
    abundance statistics are computed only on those runs. Otherwise all
    non-blank runs are used. Use this to make the pair-level log2fc
    match the realized effect size on a specific benchmark subset.
    """
    df_ref = df.filter(pl.col(COL_COND) != blank_cond)
    if reference_runs is not None:
        df_ref = df_ref.filter(pl.col(COL_RUN).is_in(list(reference_runs)))

    # Mean precursor intensity (over runs) per (protein, precursor)
    prec_mean = (
        df_ref.group_by(COL_PG, COL_PREC)
              .agg(prec_mean=pl.col("FG.Quantity").mean())
              .drop_nulls("prec_mean")
              .filter(pl.col("prec_mean") > 0)
    )

    # Peptide -> peptide mean (proxy by FG.Quantity averaged within peptide)
    pep_mean = (
        df_ref.group_by(COL_PG, COL_PEP)
              .agg(pep_mean=pl.col("FG.Quantity").mean())
              .drop_nulls("pep_mean")
    )

    n_peptides = pep_mean.group_by(COL_PG).agg(n_peptides=pl.col(COL_PEP).n_unique())

    prot_stats = (
        prec_mean.group_by(COL_PG)
                 .agg(
                     prot_mean_log2=pl.col("prec_mean").log(2).mean(),
                     n_precursors=pl.col(COL_PREC).n_unique(),
                 )
                 .join(n_peptides, on=COL_PG, how="inner")
                 .filter(pl.col("n_precursors") >= min_precursors)
    )
    return prot_stats, prec_mean, pep_mean


# ---------------------------------------------------------------------------
# Step 2 + 3 + 4 — pools, candidate pairs, greedy matching
# ---------------------------------------------------------------------------

def build_pairs(
    prot_stats: pl.DataFrame,
    swap_fraction: float,
    min_log2fc: float,
    seed: int,
) -> pl.DataFrame:
    """Stratified-by-`n_precursors` greedy pair selection.

    For each distinct `n_precursors` value (processed high → low), this
    routine targets `round(swap_fraction * bin_pop)` pairs, where
    `bin_pop` is the number of proteins with that `n_precursors` in
    `prot_stats`. Within each bin, candidates are sorted by:

        1. (log2fc - min_log2fc) ascending — pairs closest to the
           minimum-effect threshold are picked first, then the matcher
           expands outward through larger log2fc as needed.
        2. abs(n_peptides_hi - n_peptides_lo) ascending — prefer pairs
           with similar peptide counts as a tiebreaker.

    Any quota a bin cannot fill (no FC-eligible candidates whose two
    proteins are still unused) **spills down** to the next bin
    immediately below. This keeps the Positive `n_precursors`
    distribution proportional to the protein-universe distribution in
    expectation, while keeping the total pair count close to
    `round(swap_fraction * n_proteins)`.

    Bin candidates are pre-shuffled with `seed` so that equal-quality
    rows are picked uniformly within ties.
    """
    n_proteins = prot_stats.height
    target_n_pairs = max(1, round(swap_fraction * n_proteins))
    print(f"[pairs] target_n_pairs = {target_n_pairs} (= {swap_fraction:.0%} of "
          f"{n_proteins} eligible proteins)", file=sys.stderr)

    # Self-join on n_precursors, keep ordered pairs where left is higher-abundance.
    left = prot_stats.rename({c: f"{c}_hi" for c in prot_stats.columns if c != "n_precursors"})
    right = prot_stats.rename({c: f"{c}_lo" for c in prot_stats.columns if c != "n_precursors"})
    candidates = (
        left.join(right, on="n_precursors", how="inner")
            .filter(pl.col(f"{COL_PG}_hi") != pl.col(f"{COL_PG}_lo"))
            .with_columns(log2fc=pl.col("prot_mean_log2_hi") - pl.col("prot_mean_log2_lo"))
            .filter(pl.col("log2fc") >= min_log2fc)
            .with_columns(
                pep_count_diff=(pl.col("n_peptides_hi") - pl.col("n_peptides_lo")).abs(),
                fc_above_min=pl.col("log2fc") - min_log2fc,
            )
    )
    print(f"[pairs] {candidates.height} candidate ordered pairs with log2fc >= "
          f"{min_log2fc}", file=sys.stderr)

    # Per-n_precursors bin populations (count proteins, NOT candidate rows).
    bin_pop = (prot_stats.group_by("n_precursors")
                          .agg(pl.len().alias("bin_pop"))
                          .sort("n_precursors", descending=True))

    used: set[str] = set()
    rows: list[dict] = []
    spill = 0
    for bin_row in bin_pop.iter_rows(named=True):
        n = bin_row["n_precursors"]
        pop = bin_row["bin_pop"]
        bin_quota_self = int(round(swap_fraction * pop))
        bin_quota = bin_quota_self + spill

        # Candidates restricted to this n_precursors bin; shuffle for tie-break,
        # then sort so the smallest-fc-above-min comes first.
        bin_cands = (
            candidates.filter(pl.col("n_precursors") == n)
                       .sample(fraction=1.0, shuffle=True, seed=seed)
                       .sort(["fc_above_min", "pep_count_diff"])
        )

        filled = 0
        for row in bin_cands.iter_rows(named=True):
            a, b = row[f"{COL_PG}_hi"], row[f"{COL_PG}_lo"]
            if a in used or b in used:
                continue
            used.add(a); used.add(b)
            rows.append(row)
            filled += 1
            if filled >= bin_quota:
                break

        unmet = max(0, bin_quota - filled)
        print(f"[pairs]   n_precursors={n:>3}  bin_pop={pop:>5}  "
              f"quota={bin_quota_self}+spill={spill}={bin_quota}  "
              f"filled={filled}  spill_out={unmet}", file=sys.stderr)
        spill = unmet

    if not rows:
        return pl.DataFrame()
    pairs = pl.DataFrame(rows).with_row_index(name="pair_id")
    print(f"[pairs] {pairs.height} matched pairs after stratified selection "
          f"(trailing unmet spill: {spill})", file=sys.stderr)
    return pairs


# ---------------------------------------------------------------------------
# Step 5 — rank tables (precursor and peptide) per pair
# ---------------------------------------------------------------------------

def build_rank_pairs(
    pairs: pl.DataFrame,
    item_mean: pl.DataFrame,
    item_col: str,
    mean_col: str,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Rank-align items within paired proteins.

    Returns (matched, dropped). `matched` has columns
    `pair_id, rank, item_hi, item_lo`. `dropped` lists items dropped because
    their pair partner has fewer items at that rank.
    """
    ranked = item_mean.with_columns(
        rank=pl.col(mean_col).rank(method="ordinal", descending=True)
                              .over(COL_PG).cast(pl.Int32),
    )

    hi = (pairs.select(["pair_id", f"{COL_PG}_hi"])
                .join(ranked.rename({COL_PG: f"{COL_PG}_hi", item_col: f"{item_col}_hi"}),
                      on=f"{COL_PG}_hi", how="inner")
                .select("pair_id", "rank", f"{item_col}_hi"))
    lo = (pairs.select(["pair_id", f"{COL_PG}_lo"])
                .join(ranked.rename({COL_PG: f"{COL_PG}_lo", item_col: f"{item_col}_lo"}),
                      on=f"{COL_PG}_lo", how="inner")
                .select("pair_id", "rank", f"{item_col}_lo"))

    joined = hi.join(lo, on=["pair_id", "rank"], how="full", coalesce=True)
    matched = joined.filter(
        pl.col(f"{item_col}_hi").is_not_null() & pl.col(f"{item_col}_lo").is_not_null()
    )
    dropped_hi = (joined.filter(pl.col(f"{item_col}_lo").is_null())
                        .select("pair_id", pl.col(f"{item_col}_hi").alias("item"),
                                pl.lit("hi").alias("side")))
    dropped_lo = (joined.filter(pl.col(f"{item_col}_hi").is_null())
                        .select("pair_id", pl.col(f"{item_col}_lo").alias("item"),
                                pl.lit("lo").alias("side")))
    dropped = pl.concat([dropped_hi, dropped_lo])
    return matched, dropped


# ---------------------------------------------------------------------------
# Step 6 — random G1/G2 split per dilution
# ---------------------------------------------------------------------------

def assign_groups(df: pl.DataFrame, blank_cond: str, seed: int) -> pl.DataFrame:
    runs = (
        df.filter(pl.col(COL_COND) != blank_cond)
          .select(COL_RUN, COL_COND)
          .unique()
          .sort(COL_RUN)
    )
    # Reproducible random group assignment within each condition.
    runs = runs.with_columns(
        rnd=pl.lit(0).cast(pl.Int32),
    )
    # Use polars' shuffle with seed per condition by sampling row indices.
    out = []
    for cond in runs[COL_COND].unique().sort().to_list():
        sub = runs.filter(pl.col(COL_COND) == cond).sample(fraction=1.0, shuffle=True, seed=seed)
        sub = sub.with_columns(
            group=pl.Series(["G1" if i % 2 == 0 else "G2" for i in range(sub.height)])
        )
        out.append(sub)
    df_groups = pl.concat(out).drop("rnd")
    # Blanks stay ungrouped (untouched, kept in output as-is).
    return df_groups


# ---------------------------------------------------------------------------
# Step 7 — apply the swap
# ---------------------------------------------------------------------------

def apply_swap(
    df: pl.DataFrame,
    groups: pl.DataFrame,
    pairs: pl.DataFrame,
    prec_match: pl.DataFrame,
    pep_match: pl.DataFrame,
    dropped_precursors: pl.DataFrame,
    dropped_peptides: pl.DataFrame,
) -> pl.DataFrame:
    # Tag every row with its group ("G1" / "G2" / null for blanks).
    df = df.join(groups.select(COL_RUN, "group"), on=COL_RUN, how="left")

    # Drop surplus precursors / peptides entirely (across all runs, all dilutions).
    # TODO: scope this drop to the paired proteins only. Currently the filter
    # removes the precursor/peptide IDs everywhere — if any dropped ID is
    # shared with a non-pair (Negative) protein, that protein loses rows
    # unfairly. Tryptic peptides are usually protein-specific so the
    # real-world impact is small, but the filter should restrict to rows
    # where PG.ProteinGroups is one of the pair members involved.
    drop_prec_ids = set(dropped_precursors["item"].to_list())
    drop_pep_ids = set(dropped_peptides["item"].to_list())
    if drop_prec_ids:
        df = df.filter(~pl.col(COL_PREC).is_in(list(drop_prec_ids)))
    if drop_pep_ids:
        df = df.filter(~pl.col(COL_PEP).is_in(list(drop_pep_ids)))
    print(f"[swap] dropped {len(drop_prec_ids)} surplus precursors, "
          f"{len(drop_pep_ids)} surplus peptides", file=sys.stderr)

    # --- Build directional swap tables (each maps self -> partner) ---
    # Precursor (and FG-level) swap: rank-matched precursor -> partner precursor.
    prec_swap = pl.concat([
        prec_match.select(
            pl.col(f"{COL_PREC}_hi").alias("self_prec"),
            pl.col(f"{COL_PREC}_lo").alias("partner_prec"),
        ),
        prec_match.select(
            pl.col(f"{COL_PREC}_lo").alias("self_prec"),
            pl.col(f"{COL_PREC}_hi").alias("partner_prec"),
        ),
    ]).unique()

    pep_swap = pl.concat([
        pep_match.select(
            pl.col(f"{COL_PEP}_hi").alias("self_pep"),
            pl.col(f"{COL_PEP}_lo").alias("partner_pep"),
        ),
        pep_match.select(
            pl.col(f"{COL_PEP}_lo").alias("self_pep"),
            pl.col(f"{COL_PEP}_hi").alias("partner_pep"),
        ),
    ]).unique()

    pg_swap = pl.concat([
        pairs.select(
            pl.col(f"{COL_PG}_hi").alias("self_pg"),
            pl.col(f"{COL_PG}_lo").alias("partner_pg"),
        ),
        pairs.select(
            pl.col(f"{COL_PG}_lo").alias("self_pg"),
            pl.col(f"{COL_PG}_hi").alias("partner_pg"),
        ),
    ]).unique()

    # --- Fragment rank within (run, precursor) ---
    df = df.with_columns(
        frag_rank=pl.col("F.PeakArea")
                    .rank(method="ordinal", descending=True)
                    .over([COL_RUN, COL_PREC])
                    .cast(pl.Int32),
    )

    # --- Build per-run partner intensity lookups ---
    # Fragment- and FG-level partner table: keyed on (run, partner_prec, frag_rank)
    df_partner_frag = df.select(
        pl.col(COL_RUN),
        pl.col(COL_PREC).alias("partner_prec"),
        pl.col("frag_rank"),
        *[pl.col(c).alias(c + "_NEW") for c in F_COLS + FG_COLS],
    )

    # Peptide-level partner: keyed on (run, partner_pep). PEP.Quantity is constant
    # within (run, pep), so take first per group.
    df_partner_pep = (
        df.group_by([COL_RUN, COL_PEP]).agg(pl.col("PEP.Quantity").first())
          .rename({COL_PEP: "partner_pep", "PEP.Quantity": "PEP.Quantity_NEW"})
    )

    # Protein-level partner: keyed on (run, partner_pg).
    df_partner_pg = (
        df.group_by([COL_RUN, COL_PG]).agg(pl.col("PG.Quantity").first())
          .rename({COL_PG: "partner_pg", "PG.Quantity": "PG.Quantity_NEW"})
    )

    # --- Apply the swap: only G2 rows whose precursor is in prec_swap ---
    df = (
        df.join(prec_swap, left_on=COL_PREC, right_on="self_prec", how="left")
          .join(pep_swap, left_on=COL_PEP, right_on="self_pep", how="left")
          .join(pg_swap, left_on=COL_PG, right_on="self_pg", how="left")
    )

    is_g2 = pl.col("group") == "G2"

    # Fragment + FG: join partner data at (run, partner_prec, frag_rank).
    df = df.join(df_partner_frag, on=[COL_RUN, "partner_prec", "frag_rank"], how="left")

    # Drop fragment rows that are G2-swapped but have no partner at their rank
    # (their partner precursor has fewer fragments in that run).
    mask_drop_frag = is_g2 & pl.col("partner_prec").is_not_null() & pl.col("F.PeakArea_NEW").is_null()
    df = df.filter(~mask_drop_frag)

    # Peptide partner
    df = df.join(df_partner_pep, on=[COL_RUN, "partner_pep"], how="left")
    # Protein partner
    df = df.join(df_partner_pg, on=[COL_RUN, "partner_pg"], how="left")

    # Apply: for G2 swapped rows, replace original with _NEW.
    for c in F_COLS + FG_COLS + PEP_COLS + PG_COLS:
        new_c = c + "_NEW"
        df = df.with_columns(
            pl.when(is_g2 & pl.col(new_c).is_not_null())
              .then(pl.col(new_c))
              .otherwise(pl.col(c))
              .alias(c)
        )

    # Clean up helper columns
    helper_cols = ["group", "partner_prec", "partner_pep", "partner_pg", "frag_rank",
                   *[c + "_NEW" for c in F_COLS + FG_COLS + PEP_COLS + PG_COLS]]
    df = df.drop([c for c in helper_cols if c in df.columns])
    return df


# ---------------------------------------------------------------------------
# Step 8 — ground truth & true-positive list
# ---------------------------------------------------------------------------

def build_ground_truth(
    pairs: pl.DataFrame,
    prec_match: pl.DataFrame,
    pep_match: pl.DataFrame,
    dropped_precursors: pl.DataFrame,
    dropped_peptides: pl.DataFrame,
    groups: pl.DataFrame,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    n_prec_used = prec_match.group_by("pair_id").agg(n_precursors_used=pl.len())
    n_pep_used = pep_match.group_by("pair_id").agg(n_peptides_used=pl.len())
    dropped_prec_per_pair = (
        dropped_precursors.group_by("pair_id")
                          .agg(dropped_precursors=pl.col("item").unique().sort())
    )
    dropped_pep_per_pair = (
        dropped_peptides.group_by("pair_id")
                        .agg(dropped_peptides=pl.col("item").unique().sort())
    )

    g_runs = (
        groups.group_by(COL_COND, "group")
              .agg(runs=pl.col(COL_RUN).sort())
              .pivot(on="group", values="runs", index=COL_COND)
    )

    gt = (
        pairs.join(n_prec_used, on="pair_id", how="left")
             .join(n_pep_used, on="pair_id", how="left")
             .join(dropped_prec_per_pair, on="pair_id", how="left")
             .join(dropped_pep_per_pair, on="pair_id", how="left")
             .with_columns(
                 expected_log2fc=pl.col("log2fc"),
                 # The "hi" protein in the pair becomes DOWN-regulated in G2;
                 # the "lo" protein becomes UP-regulated.
             )
             .select(
                 "pair_id",
                 pl.col(f"{COL_PG}_hi").alias("PG_high"),
                 pl.col(f"{COL_PG}_lo").alias("PG_low"),
                 "n_precursors",
                 "n_peptides_hi",
                 "n_peptides_lo",
                 "n_precursors_used",
                 "n_peptides_used",
                 "expected_log2fc",
                 "dropped_precursors",
                 "dropped_peptides",
             )
    )

    tp_up = gt.select(
        pl.col("PG_low").alias(COL_PG),
        pl.lit("up").alias("role"),
        "pair_id",
        pl.col("PG_high").alias("partner_PG"),
        pl.col("expected_log2fc").alias("expected_log2fc"),
        "n_precursors_used",
        "n_peptides_used",
    )
    tp_down = gt.select(
        pl.col("PG_high").alias(COL_PG),
        pl.lit("down").alias("role"),
        "pair_id",
        pl.col("PG_low").alias("partner_PG"),
        (-pl.col("expected_log2fc")).alias("expected_log2fc"),
        "n_precursors_used",
        "n_peptides_used",
    )
    tp = pl.concat([tp_up, tp_down]).sort([COL_PG])
    return gt, tp, g_runs


# ---------------------------------------------------------------------------
# End-of-run summary log
# ---------------------------------------------------------------------------

def summary_log(df_swapped: pl.DataFrame,
                  groups: pl.DataFrame,
                  swap_list: pl.DataFrame,
                  pairs: pl.DataFrame,
                  blank_cond: str) -> None:
    """Print a one-line-per-metric summary to stderr.

    Reports n_TP, n_TN, median |log2FC| over selected pairs, and the
    median within-G1 / within-G2 SD of log2 protein mean intensity
    for TP and TN proteins. SDs are computed on the **whole** swapped
    `Report.tsv` universe (no good_data filtering), per the user's
    spec in `TODO/TODO_updated_swap.md`.
    """
    n_tp = int(swap_list.filter(pl.col("Label") == "Positive").height)
    n_tn = int(swap_list.filter(pl.col("Label") == "Negative").height)
    median_log2fc = float(pairs["log2fc"].median())

    # Aggregate fragment-level FG.Quantity to one (Run, Protein) log2 value.
    pld = (
        df_swapped.filter(pl.col(COL_COND) != blank_cond)
                   .filter(pl.col("FG.Quantity").is_not_null())
                   .filter(pl.col("FG.Quantity") > 0)
                   .group_by([COL_RUN, COL_PG, COL_PREC])
                   .agg(pl.col("FG.Quantity").first())  # precursor-level, constant
                   .group_by([COL_RUN, COL_PG])
                   .agg(log2I=pl.col("FG.Quantity").log(2).mean())
    )
    pld_g = pld.join(groups.select(COL_RUN, "group"), on=COL_RUN, how="inner")
    pld_lab = pld_g.join(swap_list.rename({"Protein": COL_PG}), on=COL_PG, how="inner")

    sd_per_prot = (
        pld_lab.group_by([COL_PG, "Label", "group"])
                .agg(sd_log=pl.col("log2I").std(), n=pl.len())
                .filter(pl.col("n") >= 2)
    )

    def med(label: str, g: str) -> float:
        s = sd_per_prot.filter(
            (pl.col("Label") == label) & (pl.col("group") == g)
        )["sd_log"]
        if s.len() == 0:
            return float("nan")
        return float(s.median())

    sd_tp_g1 = med("Positive", "G1")
    sd_tp_g2 = med("Positive", "G2")
    sd_tn_g1 = med("Negative", "G1")
    sd_tn_g2 = med("Negative", "G2")

    print("", file=sys.stderr)
    print(f"[summary] n_TP                       = {n_tp}", file=sys.stderr)
    print(f"[summary] n_TN                       = {n_tn}", file=sys.stderr)
    print(f"[summary] target |log2FC|, median    = {median_log2fc:.3f}",
          file=sys.stderr)
    print(f"[summary] within-G1 SD of TP, median = {sd_tp_g1:.3f}", file=sys.stderr)
    print(f"[summary] within-G2 SD of TP, median = {sd_tp_g2:.3f}", file=sys.stderr)
    print(f"[summary] within-G1 SD of TN, median = {sd_tn_g1:.3f}", file=sys.stderr)
    print(f"[summary] within-G2 SD of TN, median = {sd_tn_g2:.3f}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    report_stem = args.report.stem
    out_report = args.out_dir / args.report.name
    out_gt = args.out_dir / f"{report_stem}_swap_ground_truth.tsv"
    out_tp = args.out_dir / f"{report_stem}_swap_true_positives.tsv"
    out_grp = args.out_dir / f"{report_stem}_swap_group_annotation.csv"

    df = load_report(args.report)
    original_cols = df.columns
    original_rows = df.height

    # Optionally restrict the abundance-reference run set so the pair-level
    # log2fc is measured on the same subset that build_subsets/good_data
    # picks downstream. Without this, the script measures abundance across
    # the full dilution series, which can make pair-level log2fc and the
    # realized per-protein Cond1-Cond2 effect diverge.
    reference_runs: list[str] | None = None
    if args.good_rule is not None:
        ann_for_filter = pl.read_csv(args.annotation)
        good_runs = filter_runs_by_good_rule(ann_for_filter, args.good_rule)
        reference_runs = good_runs[COL_RUN].to_list()
        print(f"[stats] --good-rule={args.good_rule}: restricting reference "
              f"runs to {len(reference_runs)} of {ann_for_filter.height}",
              file=sys.stderr)

    prot_stats, prec_mean, pep_mean = compute_protein_stats(
        df, args.blank_condition, min_precursors=args.min_precursors,
        reference_runs=reference_runs,
    )
    print(f"[stats] {prot_stats.height} proteins after n_precursors >= "
          f"{args.min_precursors} filter", file=sys.stderr)

    pairs = build_pairs(prot_stats, args.swap_fraction,
                        args.min_log2fc, args.seed)
    if pairs.height == 0:
        print("[error] no pairs survived selection — adjust --min-log2fc or "
              "--min-precursors", file=sys.stderr)
        return 1

    prec_match, dropped_precursors = build_rank_pairs(pairs, prec_mean, COL_PREC, "prec_mean")
    pep_match, dropped_peptides = build_rank_pairs(pairs, pep_mean, COL_PEP, "pep_mean")
    print(f"[ranks] matched precursor pairs: {prec_match.height} "
          f"(dropped {dropped_precursors.height})", file=sys.stderr)
    print(f"[ranks] matched peptide pairs:   {pep_match.height} "
          f"(dropped {dropped_peptides.height})", file=sys.stderr)

    groups = assign_groups(df, args.blank_condition, args.seed)

    df_swapped = apply_swap(df, groups, pairs, prec_match, pep_match,
                            dropped_precursors, dropped_peptides)

    # Restore column order to match the input exactly.
    df_swapped = df_swapped.select(original_cols)
    print(f"[out] {df_swapped.height:,} rows ({original_rows - df_swapped.height:,} dropped)",
          file=sys.stderr)

    gt, tp, g_runs = build_ground_truth(pairs, prec_match, pep_match,
                                         dropped_precursors, dropped_peptides, groups)

    print(f"[write] {out_report}", file=sys.stderr)
    df_swapped.write_csv(out_report, separator="\t", null_value="NaN",
                          include_header=True, quote_style="never",
                          line_terminator="\r\n")

    print(f"[write] {out_gt}", file=sys.stderr)
    def _list_to_json(x: object) -> str:
        # Polars list-column cells come through as pl.Series (not list) under
        # polars >= 1.x; coerce to a Python list before json.dumps. `x or []`
        # used to work here but raises on a Series ("truth value is ambiguous").
        if x is None:
            return "[]"
        return json.dumps(list(x))
    gt.with_columns(
        pl.col("dropped_precursors").map_elements(_list_to_json, return_dtype=pl.String),
        pl.col("dropped_peptides").map_elements(_list_to_json, return_dtype=pl.String),
    ).write_csv(out_gt, separator="\t")

    print(f"[write] {out_tp}", file=sys.stderr)
    tp.write_csv(out_tp, separator="\t")

    print(f"[write] {out_grp}", file=sys.stderr)
    groups.write_csv(out_grp)

    # Canonical companions for the Makefile pipeline.
    # IMPORTANT: protein-swap's design contrast is G1 vs G2 (the random
    # split within each dilution), NOT the original dilution-based Condition.
    # So we rewrite the annotation's Condition column from the swap group
    # assignment: G1 -> Condition1, G2 -> Condition2. Other columns
    # (BioReplicate, Order, Label) come from the input annotation unchanged.
    out_annotation = args.out_dir / "annotation.csv"
    ann_orig = pl.read_csv(args.annotation)
    ann_new = (
        groups.select(["R.FileName", "group"])
              .with_columns(
                  pl.when(pl.col("group") == "G1").then(pl.lit("Condition1"))
                    .when(pl.col("group") == "G2").then(pl.lit("Condition2"))
                    .otherwise(pl.col("group"))
                    .alias("Condition")
              )
              .drop("group")
              .join(ann_orig.drop("Condition"), on="R.FileName", how="left")
    )
    print(f"[write] {out_annotation}", file=sys.stderr)
    ann_new.write_csv(out_annotation)

    # CSF_protein_swap_list.csv: every protein in prot_stats (universe used
    # by pairing) gets a Positive/Negative label. Positives are pair members.
    positives = set(pairs[f"{COL_PG}_hi"].to_list()) | set(pairs[f"{COL_PG}_lo"].to_list())
    swap_list = (
        prot_stats.select(pl.col(COL_PG).alias("Protein"))
        .unique()
        .with_columns(
            pl.when(pl.col("Protein").is_in(list(positives)))
              .then(pl.lit("Positive"))
              .otherwise(pl.lit("Negative"))
              .alias("Label")
        )
    )
    out_swap_list = args.out_dir / "CSF_protein_swap_list.csv"
    print(f"[write] {out_swap_list}", file=sys.stderr)
    swap_list.write_csv(out_swap_list)

    # End-of-run summary: TP/TN counts, median |log2FC|, median
    # within-group SDs for TP and TN. Computed on the swapped report
    # universe (whole Report.tsv, not the good_data subset).
    summary_log(df_swapped, groups, swap_list, pairs, args.blank_condition)

    return 0


if __name__ == "__main__":
    sys.exit(main())
