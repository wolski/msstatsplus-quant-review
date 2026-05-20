# §4 CSF_Spectronaut_protein_swap restructure & sweep

Plan: TODO §4 in [../TODO/TODO_step_by_step.md](../TODO/TODO_step_by_step.md).

## What was done

1. **Renamed folder** on disk: `CSF_Spectronaut_swap/` → `CSF_Spectronaut_protein_swap/` via `git mv`. Internal references updated in `README.md` and the two fish drivers; legacy script filenames left as-is (those scripts are superseded by `R/`, not retained as the active path).
2. **`good_data` definition corrected**: `subset_good_data` in [R/dataset_subsets.R](../R/dataset_subsets.R) now filters by `annotation$Label == "Good"` (annotation flag = `neat + 1/2` runs only), not by `Condition`. The previous read of `Condition` was wrong — `Condition` in this annotation is the `Condition1`/`Condition2` design factor, not the raw dilution.
3. **MSstats wrapper fixes** in [R/models_msstats.R](../R/models_msstats.R):
   - Aliased `R.FileName → Run` on the annotation before `SpectronauttoMSstatsFormat` (the function requires `Run`).
   - Pass `runOrder = unique(annotation[, .(Run, Order)])` in the MSstats+ branch (anomaly-score temporal features need it).
4. **Path adjustment** in [R/run_cell.R](../R/run_cell.R): protein_swap report path → the actual on-disk filename (`20250130_..._Report.tsv`).
5. **Legacy folder README** at [CSF_Spectronaut_protein_swap/README.md](../CSF_Spectronaut_protein_swap/README.md) now points users at `R/run_cell.R`; legacy scripts marked as superseded.
6. **`noswap` branch** — N/A in protein_swap. The legacy scripts under this folder only ever ran on the already-swapped TSV; there was never a `noswap` branch here. The "disable noswap" item in TODO §4.2 was a no-op for this folder and belongs to §3 (CSF_Spectronaut) where both swapped and pre-swapped runners coexist.

## Sweep results — 84 cells (4 datasets × 3 normalizations × 7 packages)

| Status | Count | Detail |
|---|---|---|
| ok | **79** | Wrote `<pkg>_model.csv` to `CSF_Spectronaut_protein_swap/<dataset>/<normalization>/swap/<pkg>/` |
| known-bad | 4 | MSstats+ × `quantile` × all four datasets — internal MSstats `dataProcess` crash after summarization (`merge.data.table ... missing from y: [FEATURE, RUN, cen]`). Reproduces in legacy: `all_dilutions/v3_quantile/MSstats+/` contains only the input CSV, no model — same bug pre-rename. Bug is in `MSstats::dataProcess(normalization = "quantile")`, not in our wrapper. |
| known-bad | 1 | DEqMS × `small_good_data` × `log2` — `spectraCounteBayes` calls `loess` which fails on 4-vs-4 data (`NA/NaN/Inf in foreign function call`). Real DEqMS edge case: the variance-shrinkage step needs more replicates than 4 per group. Other small_good_data normalizations work; only `log2` triggers it. |

The 5 known-bad cells produce no model CSV; downstream Table 1 builders treat them as missing rows.

## Numerical regression check vs legacy `V1_log2/MSstats/`

`Rscript` diff against `CSF_Spectronaut_protein_swap/all_dilutions/V1_log2/MSstats/MSstats_model.csv`:

```
rows new=2673 old=2673 common=2673
max |logFC.new - log2FC.old|         = 0.000e+00
max |pvalue.new - pvalue.old|        = 0.000e+00
max |adj.pvalue.new - adj.pvalue.old| = 0.000e+00
```

Bit-for-bit identical across all 2673 proteins. The rewrite preserves the legacy semantics exactly.

## Sweep timings

Initial sweep wall-clock: 10:55 → 14:37 ≈ 3h 41m for 84 cells (single-threaded; MSstats `n_cores = 1L`).

Approximate per-package medians from the log:
- MSstats/MSstats+ on `all_data`: 8–22 min/cell (the bottleneck).
- limma/DEqMS/prolfqua: 5–60 sec/cell.
- msqrob2: 1–4 min/cell.
- limpa: 1–2 min/cell.

If we ever need to compress the wall-clock, the route is `n_cores > 1` for MSstats; everything else is already single-threaded by design.

## Files written

- 79 × `CSF_Spectronaut_protein_swap/<dataset>/<normalization>/swap/<package>/<pkg>_model.csv` (+ timing CSV).
- Per-cell stdout/stderr in [results/logs/cell_*.log](logs/) for every attempt.
- Sweep summary logs: [results/logs/protein_swap_sweep.log](logs/protein_swap_sweep.log), [results/logs/protein_swap_msstats_retry.log](logs/protein_swap_msstats_retry.log), [results/logs/protein_swap_msstatsplus_retry.log](logs/protein_swap_msstatsplus_retry.log).

## Known open items

- **MSstats+ × quantile failure** — defer to a future MSstats upgrade or report upstream. For the manuscript Table 1 we still have MSstats+ × `log2` and `median` working; the `quantile` cells will show "missing" in the table.
- **DEqMS × small_good_data × log2** — flag this in the §7 review as a real method limitation on tiny datasets, not a bug.
- **§0.2 `good_data` rule revisit (deferred)** — the user flagged that the `<good_data>` definition for protein_swap (neat-only vs neat+1/2 via Label) is still under discussion. Current implementation uses `Label == "Good"` (neat+1/2). If the rule changes to neat-only, only `dataset_subsets.R::subset_good_data` needs to change; the sweep can be re-run for `good_data` and `small_good_data` cells only.

## Next item

§6 — bring `Mix_of_Proteome` under the shared pipeline (one dataset × 3 normalizations × 7 packages = 21 cells).
