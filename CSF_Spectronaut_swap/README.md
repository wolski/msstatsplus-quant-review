# CSF Spectronaut precursor-swap benchmark

In-silico gold-standard built from the CSF dilution Spectronaut report by
swapping precursor intensities rank-for-rank between matched protein pairs
inside a randomly chosen "G2" half of each dilution's runs. Methods are
then benchmarked on the G1 vs G2 contrast, with the un-swapped TSV serving
as a null reference.

## Running the analysis

All commands are run from `RMSV000000701.3-rerun/quant/`.

### 1. Generate the swapped report

```bash
.venv/bin/python swap_spectronaut_report.py \
  --report     "CSF_Spectronaut/20250130_163144_CSF dilutions Jan 2025 no normalization_Report.tsv" \
  --annotation "CSF_Spectronaut/CSF_annotation.csv" \
  --out-dir    "CSF_Spectronaut_swap" \
  --swap-fraction 0.05 \
  --target-log2fc-min 1.4 \
  --target-log2fc-max 1.8 \
  --seed 42
```

Emits into `CSF_Spectronaut_swap/`:

- modified Spectronaut TSV (schema-identical to the input),
- `..._swap_ground_truth.tsv` — per-pair record,
- `..._swap_true_positives.tsv` — flat TP list used downstream,
- `..._swap_group_annotation.csv` — G1/G2 assignment per run.

The annotation file `CSF_annotation.csv` and `CSF_protein_swap_list.csv` are
derived from these and are committed in `CSF_Spectronaut_swap/`.

### 2. Fit all methods × all normalizations × {post-swap, pre-swap}

```bash
cd CSF_Spectronaut_swap

# all 39 non-blank runs, write outputs into all_dilutions/
fish run_split_benchmark.fish all_dilutions ""

# drop the highest dilutions, write outputs into no_high_dilutions/
fish run_split_benchmark.fish no_high_dilutions "1to32,1to64"
```

The runner orchestrates six cells (V1_log2 / v2_vsn / v3_quantile × post / pre).
Within each cell, the MSstats and non-MSstats halves run **in parallel**.
Outputs land under `<OUT_TAG>/<variant>/<method>{,_preswap}/`. A per-method
`<method>_model.csv` and `<method>_timing.csv` are written; the latter
contains `preprocess_seconds` and `model_seconds`.

### 3. Build the comparison tables

Pure summary step — reads the model + timing CSVs, no re-fitting:

```bash
fish CSF_Spectronaut_swap_comparison_table.fish
# or for a single run:
env OUT_TAG=all_dilutions Rscript CSF_Spectronaut_swap_comparison_table.R
```

Writes `<OUT_TAG>/comparison_table.{csv,txt}`. Columns:
`Method, Variant, SwapState, TPR, PPV, preprocess_seconds, model_seconds`.

## Normalization variants

| Variant       | MSstats methods         | Non-MSstats methods |
|---------------|-------------------------|---------------------|
| `V1_log2`     | `normalization=FALSE`   | log2 only |
| `v2_vsn`      | `equalizeMedians`       | `vsn::justvsn` on raw |
| `v3_quantile` | `quantile`              | log2 + `limma::normalizeBetweenArrays(method="quantile")` |

## Per-method documentation

Each method block lives inside `run_msstats.R` or `run_nonmsstats.R` and is
wrapped in two-phase timing (`tic()` / `toc()` from `run_step_common.R`).
Inputs are precomputed once with the helpers in
`../benchmark_experiments_functions.R` (the project-wide
`prepare_data_for_<method>()` functions). G1/G2 group annotation comes from
`CSF_annotation.csv`; the TP/TN labels come from `CSF_protein_swap_list.csv`.

### MSstats+ — `run_msstats.R`

- **Preprocess:** `MSstatsConvert::SpectronauttoMSstatsFormat` with
  `intensity="PeakArea"`, `excludedFromQuantificationFilter=TRUE`,
  `filter_with_Qvalue=TRUE`, **anomaly scoring on** (`calculateAnomalyScores=TRUE`
  with FG MS1/MS2 shape quality scores + EG ΔRT, `removeMissingFeatures=.75`,
  `runOrder=run_order`, `numberOfCores=12`).
- **Summarization:** `MSstats::dataProcess` with `featureSubset="topN"`,
  `n_top_feature=100`, `MBimpute=TRUE`, `summaryMethod="linear"`,
  `normalization=<msstats_norm>`.
- **Fit + contrast:** `MSstats::groupComparison` with
  `comparison = Condition2 − Condition1`.
- **Output:** `<tag>/<variant>/MSstats+{,_preswap}/{MSstats+_input.csv,
  MSstats+_summarized.rda, MSstats+_model.csv, MSstats+_timing.csv}`.

### MSstats — `run_msstats.R`

- **Preprocess:** same as MSstats+ but with anomaly scoring OFF.
- **Summarization:** `dataProcess` with `summaryMethod="TMP"` instead of `"linear"`;
  everything else identical.
- **Fit + contrast:** `groupComparison` on the TMP-summarized data.
- **Output:** `<tag>/<variant>/MSstats{,_preswap}/...`.

### msqrob2 (hurdle) — `run_nonmsstats.R`

- **Preprocess:** `prepare_data_for_limma` → `prolfqua::LFQData` → normalize at
  precursor scale → `prolfqua::LFQDataToSummarizedExperiment` →
  `QFeatures::QFeatures` → `aggregateFeatures(fun=MsCoreUtils::robustSummary)`
  (the QFeatures vignette default — per-protein robust rollup).
- **Fit + contrast:** `msqrob2::msqrobHurdle(formula=~ 0 + group_)`
  (intensity = `MASS::rlm`; count = `glm` for MNAR proteins) +
  `msqrob2::hypothesisTestHurdle` with contrast
  `group_Condition2 − group_Condition1`. Intensity estimates preferred;
  proteins missing the intensity fit fall back to the count model.
- **No imputation** — the hurdle's count component handles MNAR natively.
- **Output:** `<tag>/<variant>/msqrob2{,_preswap}/{msqrob2_model.csv,
  msqrob2_timing.csv, msqrob_obj.rda}`.

### MaxLFQ + limma — `run_nonmsstats.R`

- **Preprocess:** `prepare_data_for_limma` → `iq::preprocess` (long → wide,
  log2, no median normalization) → `iq::fast_MaxLFQ` (protein-level log2
  estimate). For v2_vsn the log2 is undone (`2^x`) and replaced with
  `vsn::justvsn`. For v3_quantile, `limma::normalizeBetweenArrays(method="quantile")`
  is applied to the log2 protein matrix.
- **Fit + contrast:** `limma::lmFit` → `contrasts.fit(Condition2 − Condition1)` →
  `eBayes`.
- **Output:** `<tag>/<variant>/limma{,_preswap}/{limma_model.csv,
  limma_timing.csv, qc-plots.pdf}`.

### limpa — `run_nonmsstats.R`

- **Preprocess:** `prepare_data_for_limpa` (filtered wide precursor matrix
  on the linear scale) → optional vsn / log2+quantile normalization.
- **Detection probability curve + protein quantification:** `limpa::dpc` →
  `limpa::dpcQuant` (one log-intensity per protein per sample, accounting
  for MNAR via DPC).
- **Fit + contrast:** `limpa::dpcDE` → `contrasts.fit(Condition2 − Condition1)` →
  `eBayes`.
- **Note:** `dpc()` requires the log-scale variance pattern that vsn destroys;
  v2_vsn limpa is wrapped in `tryCatch` and emits an empty model on failure.
- **Output:** `<tag>/<variant>/limpa{,_preswap}/{limpa_model.csv,
  limpa_timing.csv}`.

### DEqMS — `run_nonmsstats.R`

- **Preprocess:** `prepare_data_for_deqms` (long log2 precursor table) →
  `summarize_deqms_no_ref_col` (row-median of log2 precursors per protein →
  protein matrix). For v2_vsn the log2 is undone (`2^x`) before
  `vsn::justvsn`. For v3_quantile, `normalizeBetweenArrays(method="quantile")`
  is applied to the log2 protein matrix.
- **Fit + contrast:** `limma::lmFit` → `contrasts.fit(Condition2 − Condition1)` →
  `eBayes` → `DEqMS::spectraCounteBayes` (peptide-count-weighted
  variance moderation).
- **Output:** `<tag>/<variant>/DEqMS{,_preswap}/{deqms_model.csv,
  DEqMS_timing.csv}`.

### prolfqua — `run_nonmsstats.R` → `run_prolfqua_step.R`

- **Preprocess:** `prepare_data_for_limma` → `prolfqua::LFQData` (hierarchy
  protein → precursor → fragment, including `F.FrgLossType` so duplicate
  precursor IDs don't collapse) → `intensity_array(log)` →
  `Aggregator("medpolish")` → `intensity_array(exp)` (matches the
  `prolfquapp::aggregate` canonical workflow). Protein-level normalization:
  `log2 + robscale` (V1), `vsn::justvsn` (V2), or `log2 + quantile` (V3).
- **Fit + contrast:** `prolfqua::ContrastsLMImputeFacade$new(~ group_,
  c(Condition2_vs_Condition1 = "group_Condition2 - group_Condition1"))`.
  Per-protein `lm` with LOD imputation + borrowed covariance for proteins
  missing in one group, then moderated.
- **Output:** `<tag>/<variant>/prolfqua{,_preswap}/{prolfqua_model.csv,
  prolfqua_timing.csv}`.

## File map

```
CSF_Spectronaut_swap/
├── README.md                                   this file
├── CSF_annotation.csv                          R.FileName → Condition (G1/G2)
├── CSF_protein_swap_list.csv                   Positive / Negative protein labels
├── 20250130_..._Report.tsv                     modified Spectronaut report (swapped)
├── 20250130_..._Report_swap_ground_truth.tsv   per-pair record from the swap script
├── 20250130_..._Report_swap_true_positives.tsv flat TP list
├── 20250130_..._Report_swap_group_annotation.csv  G1/G2 per run
├── CSF_precursor_swap_visualization.qmd        QC visualisation
├── CSF_precursor_swap_visualization.html
│
├── run_split_benchmark.fish                    orchestrator (this is the entry point)
├── run_msstats.R                               MSstats / MSstats+ fitter
├── run_nonmsstats.R                            msqrob2 / limma / limpa / DEqMS / prolfqua
├── run_prolfqua_step.R                         prolfqua DE block (single source of truth)
├── run_step_common.R                           shared env-var parsing, timing, vsn / quantile helpers
├── CSF_Spectronaut_swap_comparison_table.R     TPR / PPV / timing aggregator (read-only)
├── CSF_Spectronaut_swap_comparison_table.fish  wrapper that iterates over <run> tags
│
├── all_dilutions/                              outputs of the full-dataset run
│   ├── V1_log2/<method>{,_preswap}/<method>_{model,timing}.csv
│   ├── v2_vsn/.../
│   ├── v3_quantile/.../
│   ├── comparison_table.{csv,txt}
│   └── run_<variant>_{msstats,nonmsstats}{,_preswap}.log
│
└── no_high_dilutions/                          run with 1to32 + 1to64 excluded
    └── (same layout as all_dilutions/)
```

The two long-running R sources (`run_msstats.R`, `run_nonmsstats.R`) are
driven by env vars set by the fish runner:

| Variable          | Meaning                                                                  |
|-------------------|--------------------------------------------------------------------------|
| `REPORT_PATH`     | Spectronaut TSV (swapped or unswapped).                                  |
| `VARIANT`         | `V1_log2` / `v2_vsn` / `v3_quantile`.                                    |
| `OUT_SUFFIX`      | `""` (post-swap) or `"_preswap"`.                                        |
| `OUT_TAG`         | Parent run directory; defaults to `all_dilutions`.                       |
| `EXCLUDE_DILUTIONS` | Comma-separated `R.Condition` values to drop (e.g. `"1to32,1to64"`). |
| `NORMALIZATION`   | `none` / `equalizeMedians` / `quantile` / `vsn` — picked per variant.    |
