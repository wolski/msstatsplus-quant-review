# Plan: stratified pair-selection by `n_precursors`

> Status (2026-05-22): draft, **revised after user comments**. Not yet
> implemented. Companion diagnostic lives in
> [`vignettes/swap_pairs_before_after.qmd`](../vignettes/swap_pairs_before_after.qmd)
> under "Why is the Positive baseline already higher than Negative?".

## Resolved from your comments

- **`--min-precursors` CLI parameter** (default 2). Drops proteins below
  this threshold from the *entire eligible universe* (both Positives and
  Negatives), shrinks the search space, and is exposed so the threshold
  is not hidden in code.
- **Single `--min-log2fc` CLI parameter** (default 0.5) replaces the
  current `--target-log2fc-min` / `--target-log2fc-max` pair. **No upper
  bound.** This is the minimum effect size you want the benchmark to
  test.
- **Smallest-log2fc-first preference**: within each bin, the matcher
  sorts candidates by `(log2fc - min_log2fc)` ascending so pairs near
  the detection threshold are picked first; if a bin's quota cannot be
  satisfied near the minimum, the matcher keeps expanding outward
  through increasing `log2fc` until the quota is met or candidates are
  exhausted.
- **Bin spill-over** is enabled. Unfilled quota from a high-`n` bin
  rolls down to the next bin below, so the total pair count stays close
  to `round(swap_fraction × n_proteins)`. Per-bin sizes will still
  follow the universe distribution in expectation; only the bins that
  literally cannot fill their quota lose mass.
- **End-of-run logging**: after writing the swap, the script prints
  - number of TP (Positive proteins),
  - number of TN (Negative proteins),
  - median target |log2FC| over selected pairs,
  - median within-G2 SD of TP (log2 protein mean across G2 runs),
  - median within-G2 SD of TN (same metric),
  - median within-G1 SD of TP and TN (same metric).
- **Unit tests**: add a `tests/` folder under `quant/`, register pytest
  in `pyproject.toml`, and cover the core functions
  (`compute_protein_stats`, `build_pairs`, `build_rank_pairs`,
  `apply_swap`) with synthetic polars frames.

## Motivation (unchanged)

The current pair-builder in
[`src/swap_spectronaut_report.py`](../src/swap_spectronaut_report.py)
exhausts a global greedy queue and ends up enriched 3× for low-precursor
proteins:

| | Positives | Negatives |
|---|---|---|
| median `n_precursors`   | **3** | **9** |
| median log2 abundance   | 5.48  | 6.73  |

Because protein-level intensity per run is the mean of `log2(FG.Quantity)`
over precursors, a protein with ~3 precursors has a noisier per-run mean
than one with ~9. The Positive class therefore carries a higher
within-condition SD already in the **pre-swap** data on the same 15 runs.
That's a selection effect, not a swap effect. The fix is to make the
Positive `n_precursors` distribution match the protein-universe
distribution (which is also the Negative distribution, since Negatives =
"everyone not selected").

## What the current `build_pairs` does

```python
target_n_pairs = round(swap_fraction * n_proteins)
candidates = self_join_on_n_precursors(prot_stats)
candidates = candidates.filter(log2fc in [1.4, 1.8])
candidates = candidates.sort_by(pep_count_diff, fc_distance)  # random tie-break
# Single greedy match across ALL bins:
for row in candidates:
    if either protein already used: skip
    add row
    if len(pairs) >= target_n_pairs: break
```

Because the most-populous bins are the low-`n_precursors` bins, and the
matcher exhausts the global target before reaching the high-`n` end,
Positives end up enriched for `n_precursors = 2-3`.

## Proposed `build_pairs` (revised)

```python
# CLI parameters (defaults shown):
#   --min-precursors 2
#   --min-log2fc     0.5
#   --swap-fraction  0.05

prot_stats = compute_protein_stats(df, blank_cond, min_precursors=min_precursors)
n_proteins = prot_stats.height
target_n_pairs = round(swap_fraction * n_proteins)

candidates = self_join_on_n_precursors(prot_stats)
candidates = candidates.filter(log2fc >= min_log2fc)           # no upper bound
candidates = candidates.with_columns(
    fc_above_min = log2fc - min_log2fc,                        # 0 at minimum
)

pairs = []
spill = 0
for n in sorted(unique_n_values, descending=True):             # n >= min_precursors
    bin_pop   = number of proteins in prot_stats with n_precursors == n
    bin_quota = round(swap_fraction * bin_pop) + spill         # spill from previous bin
    bin_cands = (candidates[n_precursors == n]
                   .sample(fraction=1.0, shuffle=True, seed=seed)  # tie-break
                   .sort_by(["fc_above_min", "pep_count_diff"]))   # smallest FC first
    filled = 0
    for row in bin_cands:
        if either protein already used: continue
        pairs.append(row); used.update(...)
        filled += 1
        if filled >= bin_quota: break
    spill = max(0, bin_quota - filled)                         # roll unmet quota down
# Trailing spill (very low-n bins couldn't fill) is reported but not retried.
```

Effect:

- Positive `n_precursors` distribution matches the universe
  distribution at each `n` value, **in expectation** — exact match only
  when every bin's quota can be filled.
- Within each bin, pairs nearest the **minimum-effect threshold** are
  picked first. The benchmark therefore stresses methods at the
  detection boundary, not at an artificially-set 1.6-log2 midpoint.
- Unfilled quota rolls down to the next bin, so the total pair count
  stays at `≈ round(swap_fraction × n_proteins)` unless the search
  space genuinely runs out.
- Bins are processed high → low, so high-`n` bins (which are rare and
  more constrained by the FC requirement) get their quota first.

## Things this DOES change

- New CLI parameter `--min-precursors` (default 2). Affects both
  Positive and Negative universes via the `compute_protein_stats`
  filter.
- `--target-log2fc-min` / `--target-log2fc-max` retired in favour of
  `--min-log2fc` (default 0.5). The candidate filter becomes
  `log2fc >= min_log2fc` (no upper limit).
- `build_pairs` algorithm switched to stratified-by-`n_precursors` with
  spill-over.
- Selection preference inside a bin switched from "closest to the FC
  window midpoint" to "smallest log2FC above the minimum threshold".

## Things this does NOT change

- The downstream apply-swap, ground-truth, annotation rewrite — same.
- `--swap-fraction` semantics: `target_n_pairs = round(swap_fraction ×
  n_proteins)` ⇒ Positives ≈ `2 × swap_fraction × n_proteins`. Same as
  today.
- `apply_swap`, `build_rank_pairs`, `assign_groups`, `build_ground_truth`
  are untouched.

## Side effects you should know about

- **Invalidates the current protein-swap `Report.tsv` and everything
  downstream.** Cell stamps and bundle stamps that depend on
  `CSF_Spectronaut_protein_swap/.prep.stamp` would need to be regenerated
  (`make clean-prep` of the ps folder + re-run prep + cells). The
  currently in-flight build is unaffected as long as `prep-protein-swap`
  is not re-run.
- A few very high-`n` bins may still not fill their quota even with
  spill-over enabled — they're rare proteins that simply have no
  partner satisfying `log2fc >= min_log2fc`. The slack rolls down to
  lower-`n` bins until consumed or exhausted.
- Total pair count should be **≈ unchanged** (with spillover) compared
  to today.

## End-of-run logging (new)

After the script writes `Report.tsv`, `_swap_ground_truth.tsv`,
`_swap_true_positives.tsv`, `_swap_group_annotation.csv`,
`annotation.csv`, `CSF_protein_swap_list.csv`, it prints a summary block
to stderr:

```
[summary] n_TP                       = 224
[summary] n_TN                       = 2020
[summary] target |log2FC|, median    = 0.61
[summary] within-G1 SD of TP, median = 0.34
[summary] within-G2 SD of TP, median = 0.35
[summary] within-G1 SD of TN, median = 0.31
[summary] within-G2 SD of TN, median = 0.31
```

These let the user verify, without re-running the vignettes, that
**after** the stratified change Positive and Negative SD medians are
close, and that the median FC is near `min_log2fc` (because the matcher
now prefers small FCs).

## Files to touch

1. `src/swap_spectronaut_report.py`
   - Rewrite `build_pairs` (per-bin + spillover + min-FC sort).
   - Add `--min-precursors` and `--min-log2fc` CLI args; remove the two
     log2fc bound args.
   - Thread `min_precursors` through `compute_protein_stats`.
   - Add end-of-run summary log block.
2. `TODO/TODO_SWAP_DESIGN.md`
   - Update Step 3 (candidate table) and Step 4 (greedy matching) to
     describe the stratified, spill-over, smallest-FC-first design.
   - Replace the CLI block with the new flags.
   - Add a "Decisions captured" entry: Positive `n_precursors`
     distribution matches the universe distribution.
3. `vignettes/swap_pairs_before_after.qmd`
   - Update the n_precursors-table commentary so a future reader knows
     to check that Positive and Negative medians are close after
     re-running prep.
4. `tests/` (new directory under `quant/`)
   - `tests/__init__.py`
   - `tests/conftest.py` (fixtures: synthetic Spectronaut frames)
   - `tests/test_compute_protein_stats.py`
   - `tests/test_build_pairs.py`
   - `tests/test_apply_swap.py`
5. `pyproject.toml`
   - Add `pytest` as a dev dependency.
   - Register `tests` as the test root (or rely on default discovery).

## Unit-test outline

Test fixtures: build small synthetic polars frames of shape similar to
Spectronaut report (5-10 proteins, 6-10 precursors total, 4-6 runs).
Mark some intensities as null to stress the NaN-partner path.

- **`compute_protein_stats`**
  - returns columns `PG.ProteinGroups, prot_mean_log2, n_precursors, n_peptides`
  - drops blank-condition rows
  - respects `min_precursors` (no protein with `n_precursors < min`
    survives)
- **`build_pairs`** (the core change)
  - all returned pairs have `n_precursors >= min_precursors`
  - all returned pairs have `log2fc >= min_log2fc`
  - no protein appears in two pairs
  - per-bin counts are `≤ round(swap_fraction × bin_pop) + carry`
  - smallest-`fc_above_min` rows in a bin are picked before larger ones
  - spillover: simulate one bin with no eligible candidates and verify
    the unmet quota rolls into the next bin
- **`build_rank_pairs`**
  - rank-aligned matched table has equal HI/LO rows
  - dropped table contains the asymmetric ranks
- **`apply_swap`** (sanity, not exhaustive)
  - G1 rows are byte-equal to the input
  - G2 rows for a Positive precursor have FG.Quantity replaced with the
    partner's pre-swap FG.Quantity (matches the `self.after ==
    partner.before` check in the vignette)
  - row count = input count minus expected dropped-fragment rows

## Order of operations once approved

1. Edit `src/swap_spectronaut_report.py`.
2. Add `tests/` + extend `pyproject.toml`.
3. Run tests, make them pass.
4. Update `TODO/TODO_SWAP_DESIGN.md`.
5. Update vignette commentary.
6. Commit (does **not** run prep — live build is untouched).
7. When ready to rebuild from the new pair set:
   - `make clean-prep` (or only the protein-swap subset of it)
   - `make prep-protein-swap`
   - re-run cells (the swap-folder bundle + msstats stamps).
8. After the re-run, `swap_pairs_before_after.qmd` panels (especially
   the G2-only SD density and the n_precursors table) will validate the
   change: Positive and Negative SD densities should overlap, and the
   median `n_precursors` should match.

## Remaining open questions (after your comments)

1. **Zero-quota bins.** `round(swap_fraction × bin_pop) = 0` for very
   small bins (e.g. `n_precursors = 50+` where only one protein lives
   in the bin). Should the spill-over logic still consider that bin
   when receiving carry from above, or skip it entirely? Default
   proposal: include it (single-protein bins can't form a pair anyway,
   so the carry just continues falling through).
COMMENT: include 
2. **End-of-run SD computation scope.** Should the within-condition SD
   medians be computed on the **whole** input universe (using whatever
   protein/condition pairs exist) or only on the **good_data subset**
   used by the vignettes? Default proposal: whole `Report.tsv` universe
   for the script log; the vignette will still produce its own
   good_data-subset version for visualisation.
COMMENT: Whole.

3. **Test data fixtures.** I'll synthesise small polars frames inline.
   If you'd prefer a tiny on-disk fixture (e.g. a 100-row pinned
   Spectronaut TSV) checked into `tests/fixtures/`, say so and I'll do
   that instead.
COMMENT: polars frames inline is good.

## Notes from the user
<!-- Add further comments below this line. -->
