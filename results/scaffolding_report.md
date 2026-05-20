# §1 Shared scaffolding — what was built

Plan: [/Users/wolski/.claude/plans/jolly-petting-giraffe.md](/Users/wolski/.claude/plans/jolly-petting-giraffe.md).

## Files created

```
quant/
├── R/
│   ├── README_normalization.md        (§0.5 audit, pre-existing)
│   ├── build_manifest.R
│   ├── comparison_table.R
│   ├── dataset_subsets.R
│   ├── figures.R
│   ├── ground_truth.R
│   ├── models_deqms.R
│   ├── models_limpa.R
│   ├── models_maxlfq_limma.R
│   ├── models_msqrob2.R
│   ├── models_msstats.R       (handles both MSstats and MSstats+)
│   ├── models_prolfqua.R
│   ├── normalize.R
│   ├── paths.R
│   ├── preprocess.R
│   ├── run_cell.R
│   └── timing.R
├── vignettes/
│   ├── diagnostics.qmd
│   ├── review.qmd             (skeleton, populated in §7)
│   └── swap_visualization.qmd
├── results/
│   ├── normalization_report.md (§0.5 short)
│   └── scaffolding_report.md   (this file)
├── manifest.csv               (210 rows)
└── run_all.fish               (driver, executable)
```

## What was lifted vs. new

| Module | Lifted from | New? |
|---|---|---|
| `paths.R` | — | new (§0.1 path scheme codified) |
| `normalize.R` | [run_step_common.R:46–83](../CSF_Spectronaut_swap/run_step_common.R#L46-L83) | helpers lifted, `apply_normalization` / `msstats_normalization_arg` dispatchers new |
| `preprocess.R` | [benchmark_experiments_functions.R:113–403](../benchmark_experiments_functions.R#L113-L403) | helpers lifted, in-script swap **removed** (now upstream of the pipeline) |
| `dataset_subsets.R` | — | new (§0.2 table codified) |
| `ground_truth.R` | `label_proteins` from [run_step_common.R:114–118](../CSF_Spectronaut_swap/run_step_common.R#L114-L118) | three `truth_*` constructors new |
| `models_msstats.R` | [CSF_Spectronaut_swap/run_msstats.R](../CSF_Spectronaut_swap/run_msstats.R) | wrapped; canonical 6-column output |
| `models_maxlfq_limma.R` | [run_nonmsstats.R:138–190](../CSF_Spectronaut_swap/run_nonmsstats.R#L138-L190) | wrapped |
| `models_deqms.R` | [run_nonmsstats.R:259–327](../CSF_Spectronaut_swap/run_nonmsstats.R#L259-L327) | wrapped, writes extra `deqms_pep_count.csv` |
| `models_msqrob2.R` | [run_nonmsstats.R:28–136](../CSF_Spectronaut_swap/run_nonmsstats.R#L28-L136) | wrapped; peptide-level `center.median` retained |
| `models_prolfqua.R` | [run_prolfqua_step.R](../CSF_Spectronaut_swap/run_prolfqua_step.R) | wrapped; vsn branch removed |
| `models_limpa.R` | [run_nonmsstats.R:192–257](../CSF_Spectronaut_swap/run_nonmsstats.R#L192-L257) | wrapped; `tryCatch` around `dpc()` kept |
| `comparison_table.R` | [CSF_Spectronaut_swap_comparison_table.R](../CSF_Spectronaut_swap/CSF_Spectronaut_swap_comparison_table.R) | `tpr_ppv` lifted, table builder generalised |
| `figures.R` | density/heatmap snippets from [Mix_of_Proteome/Mix_of_Proteome_visualization.qmd](../Mix_of_Proteome/Mix_of_Proteome_visualization.qmd) | wrappers new |
| `timing.R` | `tic`/`toc`/`write_timing` from [run_step_common.R:125–140](../CSF_Spectronaut_swap/run_step_common.R#L125-L140) | lifted |
| `vignettes/*.qmd` | [V1_log_diagnostics.qmd](../CSF_Spectronaut/all_dilutions/V1_log_diagnostics.qmd), [CSF_swap_visualization.qmd](../CSF_Spectronaut/CSF_swap_visualization.qmd) | parameterised |
| `run_all.fish` + `run_cell.R` + `build_manifest.R` | — | new driver |

## Uniform contract for every `run_<pkg>()`

```r
run_<pkg>(merged_input, annotation, normalization, out_path)   # non-MSstats
run_msstats(raw_input, annotation, normalization, out_path, plus = FALSE, n_cores = 1L)
```

Returns `list(model, preprocess_seconds, model_seconds)` and writes:

- `<pkg>_model.csv` with the canonical 6-column schema `Protein, logFC, SE, DF, pvalue, adj.pvalue`. (msqrob2 also carries a `source` column = "intensity" or "count" for the hurdle blend; DEqMS additionally writes a `deqms_pep_count.csv` companion.)
- `<pkg>_timing.csv` with `method, preprocess_seconds, model_seconds`.

Output goes to `paths.R::out_dir(csffolder, dataset, normalization, package)`, which always produces `<csffolder>/<dataset>/<normalization>/swap/<package>/`.

## Verification done

- All 16 R files parse cleanly via `Rscript -e "parse('<file>')"`.
- `run_all.fish` parses cleanly via `fish -n`.
- `Rscript R/build_manifest.R` produces `manifest.csv` with **210 rows** = (4 folders × 10 datasets total) × 3 normalizations × 7 packages.
- Dry run shape check on the manifest: every row's `out_dir(...)` resolves to the §0.1 path.

## Verification not done in §1 (deferred to §3–§6 execution)

- End-to-end run of any cell — none of the modelling code was executed against real data in §1. The first time the pipeline actually runs is §4 (protein_swap restructure).
- Numerical agreement against legacy `CSF_Spectronaut_swap/all_dilutions/V1_log2/MSstats/MSstats_model.csv` — to be done when the first §4 run completes.
- Rendering `vignettes/diagnostics.qmd` — needs at least one cell with `<pkg>_model.csv` files present.

## Known limitations / open items

1. **`prepare_data_for_*` no longer swap.** The in-script swap that the original `benchmark_experiments_functions.R` performed inside each `prepare_*` is removed here. The shared pipeline assumes the TSV is **already** swapped (Python script for protein-swap, §2 for sample-swap). For the `CSF_Spectronaut` folder, the author-faithful scripts under `CSF_Spectronaut/` continue to source `../benchmark_experiments_functions.R` and use the original swap-in-script helpers — they are untouched.
2. **MSstats numericCores.** Set to `1` by default in `run_msstats.R` for reproducibility; bump via the `n_cores` argument when running large cells.
3. **No vsn branch.** §0.5 settled on three normalization labels (`log2`, `median`, `quantile`). The vsn path from the legacy scripts is not exposed in the shared pipeline; it can be added later as a fourth label without breaking anything.
4. **Manifest entry for `CSF_Spectronaut`.** §1 emits manifest rows for `CSF_Spectronaut × {all_data, good_data}`, but `CSF_Spectronaut` is supposed to run author-faithful code (§3). The driver will still hit `run_cell.R` for those rows. §3 will replace this with a thin wrapper that calls the author scripts and only redirects their output paths. For now the rows are present so the manifest reflects the full grid.
5. **Per-folder input paths in `run_cell.R`.** Hard-coded TSV/annotation names per folder. These are correct for the current folder contents but should be confirmed for the new `CSF_Spectronaut_sample_swap/` once §2 emits its TSV.

## Next item

§2 — write `quant/src/swap_spectronaut_report_samples.py` (TSV-level sample-swap script) and relocate the existing `swap_spectronaut_report.py` into `quant/src/`.
