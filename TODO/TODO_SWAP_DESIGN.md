# Gold-standard precursor-intensity swap — algorithmic plan (polars)

> **Status (2026-05-22).** This document is the historical algorithmic plan
> that `src/swap_spectronaut_report.py` was built from. It has been edited
> in place to **match the actual implementation**. Sections that describe
> an idea that did not survive contact with the data (the quantile-based
> abundance pools in former Step 2) are flagged "DROPPED" with a note on
> why. New behaviour that the script does but the plan did not anticipate
> (annotation Condition rewrite, `CSF_protein_swap_list.csv`) is added.
> The CLI block, file paths, and "Decisions captured" reflect the actual
> argparse and output layout.

## Context

Why: Existing CSF gold-standard datasets are built by physically spiking
proteins between samples. We want an **in silico** ground-truth dataset built
directly from a real Spectronaut report by swapping intensities between
high- and low-abundance proteins, so we get controlled, known fold changes on
top of a realistic measurement matrix. The output must be **schema-identical**
to the original Spectronaut TSV so every existing Spectronaut reader in this
repo (MSstats, MSstats+, DEqMS, limma, limpa, mapDIA, msqrob2) consumes it
unchanged.

Inputs:
- `RMSV000000701.3-rerun/quant/CSF_Spectronaut/20250130_163144_CSF dilutions Jan 2025 no normalization_Report.tsv` (8.84M rows × 61 cols, fragment-level).
- `RMSV000000701.3-rerun/quant/CSF_Spectronaut/CSF_annotation.csv` (R.FileName → Condition mapping).

Outputs:
- `RMSV000000701.3-rerun/quant/CSF_Spectronaut_swap/<report>.tsv` — modified report, same schema, minus dropped-precursor/peptide rows.
- `<report>_swap_ground_truth.tsv` — per-pair ground truth.
- `<report>_swap_group_annotation.csv` — G1/G2 assignment per run.

Same script later re-run on `K562_Spectronaut`.

## Language / library decision

Python + **polars** (lazy where possible). Rationale: 8.8M-row TSV streamed via
`pl.scan_csv`, joins/group-bys/window functions are the natural primitives for
this problem, and intermediate dataframes are easy to inspect.

## Key columns in the report

Identifier columns (never touched):
- `R.FileName`, `R.Condition` (dilution)
- `PG.ProteinGroups`, `PG.ProteinAccessions`
- `PEP.StrippedSequence`, `PEP.GroupingKey`
- `EG.PrecursorId`, `EG.ModifiedSequence`, `FG.Charge`
- `FG.Id`, `F.FrgIon`, `F.FrgType`, `F.FrgNum`, `F.Charge`, `F.FrgLossType`

Quantitative columns (subject to swap):
- Fragment: `F.PeakArea`, `F.PeakHeight`, `F.NormalizedPeakArea`, `F.NormalizedPeakHeight`
- Precursor: `FG.Quantity`, `FG.MS2RawQuantity`
- Peptide: `PEP.Quantity`
- Protein: `PG.Quantity`

All other columns (Q-values, RT, m/z, mass accuracy, library) pass through
unchanged.

## Algorithm (intermediate dataframes & joins)

The algorithm is built as a sequence of polars dataframes. Each `df_*` below is
a named, inspectable intermediate.

### Step 0 — Load and key the report

```
df_raw     = pl.scan_csv(report_tsv, sep="\t", infer_schema_length=10000)
df_ann     = pl.read_csv(annotation_csv)              # R.FileName, Condition, BioReplicate
df_keys    = df_raw.select(
                 R.FileName, R.Condition,
                 PG.ProteinGroups, PEP.StrippedSequence,
                 EG.PrecursorId, FG.Quantity
             ).unique()                                # ~precursor×run rows
```

Use `EG.PrecursorId` as the precursor key (peptide+charge+modifications encoded
in it). `PG.ProteinGroups` as the protein key. `PEP.StrippedSequence` as the
peptide key.

### Step 1 — Compute per-protein statistics on the reference set

Reference set = **all non-blank dilutions pooled** (`R.Condition != "blank"`).
Rationale (user choice): pooling gives the most data for stable ranking; we
exclude blanks because they have no real signal.

```
df_ref            = df_keys.filter(R.Condition != "blank")

df_prec_mean      = df_ref.group_by(PG, EG.PrecursorId).agg(
                        prec_mean    = pl.col(FG.Quantity).mean(),
                        prec_n_runs  = pl.col(FG.Quantity).len()
                    )                                      # one row per (protein, precursor)

df_prot_stats     = df_prec_mean.group_by(PG).agg(
                        prot_mean_log2  = pl.col(prec_mean).log(2).mean(),
                        n_precursors    = pl.col(EG.PrecursorId).n_unique(),
                        n_peptides      = ...                            # via peptide group_by, see below
                    )

df_pep_mean       = df_ref.group_by(PG, PEP.StrippedSequence).agg(
                        pep_mean = pl.col(FG.Quantity).mean()           # peptide intensity proxy
                    )
df_n_peptides     = df_pep_mean.group_by(PG).agg(
                        n_peptides = pl.col(PEP.StrippedSequence).n_unique()
                    )
df_prot_stats     = df_prot_stats.join(df_n_peptides, on=PG)
```

Filter: keep only proteins with `n_precursors >= 2` and `n_peptides >= 1`.

### Step 2 — Define abundance pools — **DROPPED**

> **Not implemented.** The original plan was to restrict pairing to proteins
> in the top and bottom 5 % of the abundance distribution. In practice this
> produced fold changes that were too large (and too rare a candidate pool
> to greedy-match cleanly). The implementation in
> `src/swap_spectronaut_report.py` pairs **across the entire abundance
> distribution** under the constraint that the two members are
> `log2fc_min … log2fc_max` apart (default `[1.4, 1.8]`) and share
> `n_precursors`. There are no `--abundance-quantile-high` /
> `--abundance-quantile-low` CLI flags.

### Step 3 — Build pairwise candidate table

This is the heart of the algorithm: a candidate table of every (left, right)
ordered protein pair drawn from the **full** abundance distribution that
shares precursor count and satisfies the FC window.

```
# Self-join the full prot_stats on n_precursors (exact match); keep ordered
# pairs where the left member is higher-abundance than the right.
left  = prot_stats.rename({c: f"{c}_hi" for c in prot_stats.columns if c != "n_precursors"})
right = prot_stats.rename({c: f"{c}_lo" for c in prot_stats.columns if c != "n_precursors"})
df_candidates = (
    left.join(right, on="n_precursors", how="inner")
        .filter(pl.col("PG.ProteinGroups_hi") != pl.col("PG.ProteinGroups_lo"))
        .with_columns(log2fc = pl.col("prot_mean_log2_hi") - pl.col("prot_mean_log2_lo"))
        .filter(pl.col("log2fc").is_between(log2fc_min, log2fc_max))
        .with_columns(
            pep_count_diff = (pl.col("n_peptides_hi") - pl.col("n_peptides_lo")).abs(),
            fc_distance    = (pl.col("log2fc") - (log2fc_min + log2fc_max) / 2).abs(),
        )
)
```

`df_candidates` columns: `PG.ProteinGroups_hi, PG.ProteinGroups_lo,
n_precursors, n_peptides_hi, n_peptides_lo, pep_count_diff,
prot_mean_log2_hi, prot_mean_log2_lo, log2fc, fc_distance`.

Note: the self-join makes the same unordered pair appear twice (once with A
on the left, once with B on the left). The greedy matcher below handles this
by tracking `used` proteins; we deliberately do **not** dedupe upfront, so
that the matcher can pick the orientation that survives.

### Step 4 — Greedy bipartite matching

Goal: pick a one-to-one matching (each high protein used at most once, each low
protein used at most once) of size approximately
`swap_fraction × n_proteins_total / 2`.

```
# Random shuffle first so equal-quality candidates are picked uniformly
# across the abundance distribution, then sort by quality.
df_candidates_sorted = (
    df_candidates
        .sample(fraction=1.0, shuffle=True, seed=seed)
        .sort(["pep_count_diff", "fc_distance"])
)

used: set[str] = set()       # single set; protein can be used at most once,
                              # regardless of which side
rows = []
for row in df_candidates_sorted.iter_rows(named=True):
    a, b = row["PG.ProteinGroups_hi"], row["PG.ProteinGroups_lo"]
    if a in used or b in used:
        continue
    used.add(a); used.add(b)
    rows.append(row)
    if len(rows) >= target_n_pairs:
        break

df_pairs = pl.DataFrame(rows).with_row_index(name="pair_id")
```

Target size: `target_n_pairs = max(1, round(swap_fraction * n_proteins))`,
not a hard cap on the candidate side as the plan originally framed it.

> **Fallback (±20 % precursor-count tolerance) — DROPPED.** The plan
> sketched an inequality-join fallback when the exact-`n_precursors` pool
> ran dry. The implementation does not include it; `n_precursors` is
> always an exact match. There is no `match_kind` column on `df_pairs`.

### Step 5 — Build precursor-rank and peptide-rank tables for each pair

For each pair (A, B), we need rank-aligned precursors and rank-aligned peptides
so we can swap "biggest with biggest", "second biggest with second biggest",
etc. Rank computed on the reference dilution.

```
df_prec_rank = df_prec_mean.with_columns(
    rank = pl.col(prec_mean).rank(descending=True).over(PG).cast(int)
)
# join twice — once for the high protein, once for the low
df_prec_pairs = df_pairs.select(pair_id, PG_high, PG_low).join(
                    df_prec_rank.rename({PG: PG_high, EG.PrecursorId: prec_high, prec_mean: mean_high, rank: rank}),
                    on=PG_high
                ).join(
                    df_prec_rank.rename({PG: PG_low, EG.PrecursorId: prec_low, prec_mean: mean_low, rank: rank}),
                    on=[PG_low, rank]
                )
```

Result `df_prec_pairs` has one row per matched precursor pair, columns:
`pair_id, PG_high, PG_low, rank, prec_high, prec_low, mean_high, mean_low`.

Precursor counts are equal by construction (exact-match join in Step 3),
so the precursor-rank table always matches one-to-one. **Peptide** counts
may still differ between the two pair members — `df_dropped_peptides`
captures the unmatched ranks on whichever side has the larger peptide
count. Same code path also produces a (currently empty in practice)
`df_dropped_precursors` in case of upstream nulls.

### Step 6 — Assign G1 / G2 per dilution

```
df_groups = df_ann.with_columns(
    group = pl.col(R.FileName).shuffle(seed=42).rank().over(R.Condition).mod(2).replace({0:"G1",1:"G2"})
)
```

Balanced split, seeded, per dilution. Persisted as
`<report>_swap_group_annotation.csv`.

### Step 7 — Apply the swap

Strategy: build a **swap-lookup dataframe** that, for each row in the original
report, tells us what intensity values to use. Then join it back.

#### 7a. Precursor level (FG.Quantity, FG.MS2RawQuantity)

For each G2 run, swap the two scalars between rank-matched precursors via a
self-join on `df_prec_pairs` and `R.FileName`. Straightforward.

#### 7b. Fragment level (F.PeakArea, F.PeakHeight, F.NormalizedPeakArea, F.NormalizedPeakHeight)

User directive: **ignore ion identity** (we already break the protein/sequence
↔ intensity coupling, so doing the same at the fragment level is consistent).

Rule: within each G2 run, for each precursor pair (A, B):
- Take A's fragment rows for this run, rank them by `F.PeakArea` descending.
- Take B's fragment rows for this run, rank them likewise.
- Swap the four `F.*` intensity columns rank-for-rank.
- If A has more fragment rows than B in this run, **drop the surplus rows
  from the output** (consistent with the surplus-precursor/peptide rule:
  drop rather than NA, no imputation contamination).
- Symmetric for B-surplus.

Concretely:

```
df_frag_rank = df_raw.filter(group == "G2").with_columns(
    frag_rank = pl.col("F.PeakArea").rank(descending=True)
                  .over(["R.FileName", "EG.PrecursorId"]).cast(int)
)

df_frag_swap = df_frag_rank.join(
    df_prec_pairs,
    left_on=["EG.PrecursorId"],
    right_on=["prec_side_A"]
).join(
    df_frag_rank.rename({...partner cols...}),
    left_on=["R.FileName", "prec_side_B", "frag_rank"],
    right_on=["R.FileName", "EG.PrecursorId", "frag_rank"]
).with_columns(
    F.PeakArea = F.PeakArea_partner,
    F.PeakHeight = F.PeakHeight_partner,
    F.NormalizedPeakArea = F.NormalizedPeakArea_partner,
    F.NormalizedPeakHeight = F.NormalizedPeakHeight_partner,
)
```

Inner join naturally drops surplus rows on either side (they have no partner
at that rank).

#### 7c. Peptide level

`PEP.Quantity` is constant within `(R.FileName, PG, PEP.StrippedSequence)`.
Join `df_pep_pairs` analogously and overwrite `PEP.Quantity` in G2 rows
belonging to swapped peptides.

#### 7d. Protein level

`PG.Quantity` is constant within `(R.FileName, PG)`. Build a one-row-per-
`(R.FileName, PG)` swap-target table from `df_pairs`, join, overwrite.

#### 7e. Reassemble

```
df_out = pl.concat([
    df_raw.filter(group == "G1"),                  # untouched
    df_raw.filter(group == "G2", PG not in pairs), # untouched
    df_swapped                                     # modified
]).filter(
    pl.col(EG.PrecursorId).is_in(dropped_precursors).not_()
    & pl.col(PEP.StrippedSequence).is_in(dropped_peptides).not_()
)
```

Verify: `df_out` row count == `df_raw` row count − (rows of dropped
precursors/peptides across all runs). Verify column order identical.

### Step 8 — Ground truth

`df_ground_truth` (one row per pair):
- `pair_id, match_kind (exact|fallback)`
- `PG_high, PG_low`
- `n_precursors_high, n_precursors_low, n_precursors_used`
- `n_peptides_high, n_peptides_low, n_peptides_used`
- `expected_log2fc` = mean over matched precursors of `log2(mean_low / mean_high)` (computed on G1 runs only — the unperturbed baseline)
- `realized_log2fc_per_dilution` (JSON, computed post-swap from G2/G1 on the output table — verification step)
- `dropped_precursors` (JSON list)
- `dropped_peptides` (JSON list)
- `g1_runs`, `g2_runs` per dilution (JSON)

### Step 9 — Write output

```
df_out.sink_csv(out_tsv, separator="\t", quote_style="never", null_value="")
df_ground_truth.write_csv(out_truth, separator="\t")
df_groups.write_csv(out_group_ann)
```

Verify empty-cell convention against the input file before writing
(`null_value` may need to be `"NaN"` or empty).

## Verification

1. **Schema parity:** `head -1` of input and output produce identical tab-split column lists.
2. **Row count:** `wc -l` of output = input − sum of dropped fragment rows. Logged.
3. **Spot check:** pick one pair, one run; verify G1 rows numerically identical to input; verify G2 rows for protein A carry protein B's original fragment intensities at matched ranks.
4. **Statistical sanity:** compute per-protein log2(mean G2 / mean G1) on `FG.Quantity`; swapped proteins should cluster near ±1.6, unswapped near 0.
5. **End-to-end:** run `CSF_Spectronaut_processing.R` against the new TSV — must complete without schema errors and produce a usable input for the downstream tool list.

## Files actually created (script output)

The script is at **[`src/swap_spectronaut_report.py`](../src/swap_spectronaut_report.py)**.
Invoked via the Makefile's `prep-protein-swap` target, it writes into the
`--out-dir` argument (typically `CSF_Spectronaut_protein_swap/`):

- `<report-basename>.tsv` — swapped Spectronaut report (schema-identical
  to the input, minus dropped rows). The file is **written under the
  same basename as the input**, so the symlinked `Report.tsv` in the
  swap folder points to it via the existing prep stamp logic.
- `<report-stem>_swap_ground_truth.tsv` — per-pair ground truth
  (one row per pair, with `dropped_precursors` / `dropped_peptides`
  serialised as JSON strings).
- `<report-stem>_swap_true_positives.tsv` — flat true-positive list,
  one row per touched protein, with role ∈ {`up`, `down`} and signed
  `expected_log2fc`. See "True-positive deliverable" below.
- `<report-stem>_swap_group_annotation.csv` — per-run G1/G2 assignment.
- **`annotation.csv`** — the canonical MSstats annotation used by all
  downstream model adapters. The script **rewrites the Condition column
  from G1/G2 → Condition1/Condition2** so that downstream code (which
  expects a binary Condition1 vs Condition2 contrast) works without
  modification. Other columns (`BioReplicate`, `Order`, `Label`) are
  joined in from the input annotation unchanged.
- **`CSF_protein_swap_list.csv`** — `Protein, Label` table where
  `Label = "Positive"` for any protein that is a member of a swap pair
  and `Label = "Negative"` otherwise. This is the gold-standard label
  vector that the diagnostics qmd reads as `truth_path` for
  `truth_kind=protein_swap`.

In the rerun layout, the same script is also invoked by the
sample-swap folder's prep via `src/swap_spectronaut_report_samples.py`
(a sibling script with the same shape).

CLI (actual argparse, mirrors `src/swap_spectronaut_report.py`):

```
python src/swap_spectronaut_report.py \
  --report          <input.tsv> \
  --annotation      <annotation.csv> \
  --out-dir         <dir> \
  --swap-fraction   0.05 \
  --target-log2fc-min 1.4 \
  --target-log2fc-max 1.8 \
  --blank-condition blank \
  --seed            42
```

> Removed vs the original plan: `--reference-dilution`,
> `--abundance-quantile-high`, `--abundance-quantile-low`,
> `--precursor-count-tol`. The implementation pools **all non-blank
> dilutions** as the reference set (no `--reference-dilution` knob),
> pairs across the **full abundance distribution** (no quantile knobs),
> and uses **exact `n_precursors` matching** with no inequality fallback
> (no tol knob).

## Decisions captured

- **Reference set:** all non-blank dilutions pooled. Configurable via
  `--blank-condition` only (default `blank`).
- **Swap fraction:** CLI parameter `--swap-fraction` (default 0.05).
  `target_n_pairs = max(1, round(swap_fraction * n_proteins))`. Pairs
  are drawn from the **full abundance distribution**, not from quantile
  tails.
- **Matching constraint:** equal `n_precursors` (exact); peptide counts
  may differ and surplus peptides are dropped from the output.
- **Fragment ion identity:** ignored. Fragment intensities swapped
  rank-for-rank within each (precursor, run), with surplus fragment rows
  dropped — consistent with the precursor/peptide drop rule.

## Known limitations (carry over to the next iteration)

- **`dropped_precursors` / `dropped_peptides` drop scope is global, not
  pair-scoped.** See the in-code TODO at
  [`src/swap_spectronaut_report.py:261-266`](../src/swap_spectronaut_report.py#L261).
  If a dropped precursor/peptide ID is shared with a non-pair (Negative)
  protein, that protein loses rows it should keep. Tryptic peptides are
  usually protein-specific so the real-world impact is small; the fix
  is to scope the filter to rows where `PG.ProteinGroups` is one of the
  pair members involved.
- **No `match_kind` column** on the ground truth (because the ±20 %
  fallback isn't implemented and every pair is exact-match).
- **`expected_log2fc`** on the ground-truth file is the *target* fold
  change (`prot_mean_log2_hi - prot_mean_log2_lo` measured pre-swap on
  the pooled reference set). The script does not currently re-compute
  the realised post-swap fold change as a verification field — that's
  the `realized_log2fc_per_dilution` mentioned in §"Step 8" of this
  plan and never landed in the code.

## True-positive deliverable

In addition to the per-pair ground-truth TSV, the script writes a flat
**true-positive list** — one row per touched protein — usable directly as
the gold-standard label vector by downstream evaluation scripts.

File: `<report>_swap_true_positives.tsv` with columns:
- `PG.ProteinGroups`
- `role` ∈ {`up`, `down`} — `up` if this protein received intensities from a
  higher-abundance partner (so it appears up-regulated in G2 vs G1),
  `down` otherwise.
- `pair_id` — links back to the per-pair ground-truth file.
- `partner_PG` — the other member of the pair.
- `expected_log2fc` — signed (positive for `up`, negative for `down`).
- `n_precursors_used`, `n_peptides_used`.

All proteins NOT in this list are true negatives (expected log2FC = 0).
