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

import sys
from pathlib import Path
from typing import Literal

import cyclopts
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


app = cyclopts.App(name="swap-spectronaut-report", help=__doc__,
                    help_on_error=True)


GoodRule = Literal["label_good", "neat_only"]


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

# Identity columns that define "which protein / peptide / precursor a row
# belongs to". The swap relabels these (in G2 only) to the rank-aligned
# partner's identity; everything else on the row — intensities, RT, q-values,
# fragment-ion annotations — rides along unchanged. Only the columns that
# actually exist in the report are relabeled.
IDENTITY_COLS = [
    COL_PG, "PG.ProteinAccessions",
    COL_PEP, "PEP.GroupingKey",
    COL_PREC, "EG.ModifiedSequence", "FG.Charge",
]


def build_prec_swap(prec_match: pl.DataFrame) -> pl.DataFrame:
    """Bidirectional precursor map: self_prec -> partner_prec."""
    return pl.concat([
        prec_match.select(
            pl.col(f"{COL_PREC}_hi").alias("self_prec"),
            pl.col(f"{COL_PREC}_lo").alias("partner_prec"),
        ),
        prec_match.select(
            pl.col(f"{COL_PREC}_lo").alias("self_prec"),
            pl.col(f"{COL_PREC}_hi").alias("partner_prec"),
        ),
    ]).unique()


def apply_swap(
    df: pl.DataFrame,
    groups: pl.DataFrame,
    prec_match: pl.DataFrame,
) -> pl.DataFrame:
    """Apply the swap by RELABELLING identity columns in G2.

    Instead of overwriting intensity columns in place (and dropping rows
    where the rank-aligned partner had fewer fragments), this moves whole
    rows between rank-aligned partner precursors by rewriting only their
    identity columns:

      * A G2 row that physically measured precursor Y (protein B) has its
        identity columns rewritten to Y's rank-aligned partner X
        (protein A). The row keeps Y's intensities, RT, q-values, fragment
        annotations — it is now attributed to X.
      * Symmetrically, X's G2 rows are relabelled to Y.

    Because the swap is a relabelling of existing rows, it is exact and
    symmetric:

      * Nothing is created or destroyed — total row count is preserved and
        every intensity value survives, just attributed to the partner.
      * NA holes fill naturally: if X is undetected in a G2 run but Y is
        detected, Y's rows get relabelled to X, so A gains data there and B
        loses it — the symmetric outcome. No "drop because the partner was
        NaN" asymmetry.
      * Fragment-count mismatches need no special handling: a precursor's
        whole set of fragment rows moves as a unit.

    G1 rows are never touched. `prec_match` must contain the rank-aligned
    `(EG.PrecursorId_hi, EG.PrecursorId_lo)` pairs.
    """
    id_cols = [c for c in IDENTITY_COLS if c in df.columns]
    other_id = [c for c in id_cols if c != COL_PREC]

    df = df.join(groups.select(COL_RUN, "group"), on=COL_RUN, how="left")

    prec_swap = build_prec_swap(prec_match)

    # One row per precursor with its (constant) identity-column values.
    prec_identity = df.group_by(COL_PREC).agg([pl.first(c).alias(c) for c in other_id])

    # For each self_prec, the partner precursor's identity, suffixed __new.
    partner_identity = (
        prec_swap
        .join(prec_identity, left_on="partner_prec", right_on=COL_PREC, how="inner")
        .rename({c: f"{c}__new" for c in other_id})
        .with_columns(pl.col("partner_prec").alias(f"{COL_PREC}__new"))
        .select(["self_prec", *[f"{c}__new" for c in id_cols]])
    )

    df = df.join(partner_identity, left_on=COL_PREC, right_on="self_prec", how="left")

    is_g2 = pl.col("group") == "G2"
    n_swapped_rows = df.filter(is_g2 & pl.col(f"{COL_PREC}__new").is_not_null()).height

    # Relabel identity columns on G2 swap rows. Order matters: read all __new
    # columns before overwriting (they were materialised by the join, so the
    # per-column with_columns below is safe).
    for c in id_cols:
        nc = f"{c}__new"
        df = df.with_columns(
            pl.when(is_g2 & pl.col(nc).is_not_null())
              .then(pl.col(nc))
              .otherwise(pl.col(c))
              .alias(c)
        )

    print(f"[swap] relabelled {n_swapped_rows:,} G2 rows to their partner "
          f"precursor identity (no rows dropped)", file=sys.stderr)

    helper_cols = ["group", *[f"{c}__new" for c in id_cols]]
    df = df.drop([c for c in helper_cols if c in df.columns])
    return df


# ---------------------------------------------------------------------------
# Step 8 — ground truth & true-positive list
# ---------------------------------------------------------------------------

def build_ground_truth(
    pairs: pl.DataFrame,
    prec_match: pl.DataFrame,
    groups: pl.DataFrame,
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """Per-pair ground truth + flat true-positive list + G1/G2 run map.

    The relabelling swap drops nothing, so there are no dropped-precursor
    or dropped-peptide lists to report — every paired precursor's rows are
    moved to its partner. `n_precursors_used` therefore equals the pair's
    `n_precursors`.
    """
    n_prec_used = prec_match.group_by("pair_id").agg(n_precursors_used=pl.len())

    g_runs = (
        groups.group_by(COL_COND, "group")
              .agg(runs=pl.col(COL_RUN).sort())
              .pivot(on="group", values="runs", index=COL_COND)
    )

    gt = (
        pairs.join(n_prec_used, on="pair_id", how="left")
             .with_columns(expected_log2fc=pl.col("log2fc"))
             .select(
                 "pair_id",
                 pl.col(f"{COL_PG}_hi").alias("PG_high"),
                 pl.col(f"{COL_PG}_lo").alias("PG_low"),
                 "n_precursors",
                 "n_peptides_hi",
                 "n_peptides_lo",
                 "n_precursors_used",
                 "expected_log2fc",
             )
    )

    # Direction (see apply_swap): in G2 the HIGH protein receives the LOW
    # partner's intensities, so its Cond1-Cond2 difference is +log2fc
    # ("down" in G2); the LOW protein receives the HIGH partner's
    # intensities → "up" in G2 with -log2fc as the signed Cond1-Cond2 value.
    tp_up = gt.select(
        pl.col("PG_low").alias(COL_PG),
        pl.lit("up").alias("role"),
        "pair_id",
        pl.col("PG_high").alias("partner_PG"),
        pl.col("expected_log2fc").alias("expected_log2fc"),
        "n_precursors_used",
    )
    tp_down = gt.select(
        pl.col("PG_high").alias(COL_PG),
        pl.lit("down").alias("role"),
        "pair_id",
        pl.col("PG_low").alias("partner_PG"),
        (-pl.col("expected_log2fc")).alias("expected_log2fc"),
        "n_precursors_used",
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

@app.default
def main(
    *,
    report: Path,
    annotation: Path,
    out_dir: Path,
    swap_fraction: float = 0.05,
    min_precursors: int = 2,
    min_log2fc: float = 0.5,
    good_rule: GoodRule | None = None,
    blank_condition: str = "blank",
    seed: int = 42,
) -> int:
    """Gold-standard precursor-intensity swap for a Spectronaut report.

    Parameters
    ----------
    report
        Spectronaut Report.tsv to read.
    annotation
        Annotation CSV with R.FileName + Condition (+ optional Label).
    out_dir
        Output directory; created if missing.
    swap_fraction
        Fraction of proteins to make Positive pair members. The matcher
        targets `round(swap_fraction * n_proteins)` pairs.
    min_precursors
        Drop proteins with fewer precursors than this from the eligible
        universe (Positives and Negatives).
    min_log2fc
        Minimum |log2(prot_mean_hi / prot_mean_lo)| a pair must exceed.
        No upper bound; the matcher prefers pairs closest to this floor
        and expands outward as needed.
    good_rule
        Restrict the abundance-reference run set used by
        `compute_protein_stats`. `label_good` keeps annotation rows with
        `Label == "Good"`; `neat_only` keeps R.FileName matching
        "NeatCSF". Default `None` uses all non-blank runs. Mirrors
        `src/build_subsets.py --good-rule` so the pair-level log2fc is
        measured on the same subset the good_data benchmark uses.
    blank_condition
        R.Condition value identifying blank runs (excluded from stats
        and output groups).
    seed
        RNG seed for the matcher's tie-break shuffle and the G1/G2 split.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    report_stem = report.stem
    out_report = out_dir / report.name
    out_gt = out_dir / f"{report_stem}_swap_ground_truth.tsv"
    out_tp = out_dir / f"{report_stem}_swap_true_positives.tsv"
    out_grp = out_dir / f"{report_stem}_swap_group_annotation.csv"

    df = load_report(report)
    original_cols = df.columns
    original_rows = df.height

    # Optionally restrict the abundance-reference run set so the pair-level
    # log2fc is measured on the same subset that build_subsets/good_data
    # picks downstream. Without this, the script measures abundance across
    # the full dilution series, which can make pair-level log2fc and the
    # realized per-protein Cond1-Cond2 effect diverge.
    reference_runs: list[str] | None = None
    if good_rule is not None:
        ann_for_filter = pl.read_csv(annotation)
        good_runs = filter_runs_by_good_rule(ann_for_filter, good_rule)
        reference_runs = good_runs[COL_RUN].to_list()
        print(f"[stats] --good-rule={good_rule}: restricting reference "
              f"runs to {len(reference_runs)} of {ann_for_filter.height}",
              file=sys.stderr)

    prot_stats, prec_mean, _pep_mean = compute_protein_stats(
        df, blank_condition, min_precursors=min_precursors,
        reference_runs=reference_runs,
    )
    print(f"[stats] {prot_stats.height} proteins after n_precursors >= "
          f"{min_precursors} filter", file=sys.stderr)

    pairs = build_pairs(prot_stats, swap_fraction, min_log2fc, seed)
    if pairs.height == 0:
        print("[error] no pairs survived selection — adjust --min-log2fc or "
              "--min-precursors", file=sys.stderr)
        return 1

    prec_match, _dropped_precursors = build_rank_pairs(pairs, prec_mean, COL_PREC, "prec_mean")
    print(f"[ranks] matched precursor pairs: {prec_match.height}", file=sys.stderr)

    groups = assign_groups(df, blank_condition, seed)

    df_swapped = apply_swap(df, groups, prec_match)

    # Restore column order to match the input exactly.
    df_swapped = df_swapped.select(original_cols)
    print(f"[out] {df_swapped.height:,} rows "
          f"({original_rows - df_swapped.height:,} delta vs input — expect 0)",
          file=sys.stderr)

    gt, tp, g_runs = build_ground_truth(pairs, prec_match, groups)

    print(f"[write] {out_report}", file=sys.stderr)
    df_swapped.write_csv(out_report, separator="\t", null_value="NaN",
                          include_header=True, quote_style="never",
                          line_terminator="\r\n")

    print(f"[write] {out_gt}", file=sys.stderr)
    gt.write_csv(out_gt, separator="\t")

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
    out_annotation = out_dir / "annotation.csv"
    ann_orig = pl.read_csv(annotation)
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
    out_swap_list = out_dir / "CSF_protein_swap_list.csv"
    print(f"[write] {out_swap_list}", file=sys.stderr)
    swap_list.write_csv(out_swap_list)

    # End-of-run summary: TP/TN counts, median |log2FC|, median
    # within-group SDs for TP and TN. Computed on the swapped report
    # universe (whole Report.tsv, not the good_data subset).
    summary_log(df_swapped, groups, swap_list, pairs, blank_condition)

    return 0


if __name__ == "__main__":
    app()
