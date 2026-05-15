# Step-by-step plan — msstatsplus review (RMSV000000701.3-rerun/quant)

Derived from `TODO/WhatWeStillNeed.md` plus clarifying Q&A. Out of scope: `CSF_DIANN`, `K562_Spectronaut`, `K562_DIANN`, `TPCPA` — leave untouched.

---

## 0. Conventions fixed up front

### 0.1 Unified folder path
```
<csffolder>/<dataset>/<normalization>/swap/<modellingPackage>
```
- `swap` is always the literal segment. The `noswap` branch is disabled for now (re-enabled at the very end of the project).
- `<csffolder>` ∈ {`CSF_Spectronaut`, `CSF_Spectronaut_protein_swap`, `CSF_Spectronaut_sample_swap`, `Mix_of_Proteome`}.
- `<dataset>`: see §0.2.
- `<normalization>` ∈ {`log2`, `median`, `quantile`}.
- `<modellingPackage>` ∈ {`MSstats`, `MSstatsPlus`, `MaxLFQ_limma`, `DEqMS`, `msqrob2`, `prolfqua`, `limpa`}.

### 0.2 Datasets per folder

| folder | `all_data` | `good_data` | `small` | `small_good_data` |
|---|---|---|---|---|
| `CSF_Spectronaut` | ✓ | neat + 1/2, balanced reps | — | — |
| `CSF_Spectronaut_protein_swap` | ✓ | **neat only** (TP/TN from protein swap, not dilution) | ✓ | 4 vs 4 from `good_data` |
| `CSF_Spectronaut_sample_swap` | ✓ | neat + 1/2, balanced reps | — | 4 vs 4 from `good_data` |
| `Mix_of_Proteome` | ✓ | — | — | — |

### 0.3 Code organisation
- `quant/R/` — shared R helpers (preprocessing, modelling wrappers, ground-truth utilities, table/figure code).
- `quant/vignettes/` — shared `.qmd`/`.Rmd` (visualisation, diagnostics, `review.qmd`).
- `quant/src/` — Python helpers (TSV swap scripts: existing `swap_spectronaut_report.py` and new `swap_spectronaut_report_samples.py`).
- `quant/results/` — intermediate aggregated CSVs (combined Table 1 across folders, p-value tables, etc.).
- **Exception (load-bearing):** `CSF_Spectronaut/` keeps its own author-faithful R scripts — *no* shared-code refactor of that folder. Only its output directories are reshaped to match §0.1.

### 0.4 Modelling-package list (same set everywhere)
MSstats, MSstats+, MaxLFQ+limma, DEqMS, msqrob2 (faster code), prolfqua, limpa.

### 0.5 Normalization caveat to investigate
Before running, audit whether each package's "median" normalization is median-**centering only** or median-centering **+ scaling** (z-style). Document per-package and align where feasible. Record in `R/README_normalization.md`.

---

## 1. Shared infrastructure — `quant/R/` and `quant/vignettes/`

**1.1 Inventory existing scripts to extract from.** Sources of reusable logic:
- `CSF_Spectronaut_swap/run_msstats.R`, `run_nonmsstats.R`, `run_prolfqua_step.R`, `run_step_common.R`
- `CSF_Spectronaut_swap/CSF_Spectronaut_swap_comparison_table.R`
- `Mix_of_Proteome/Mixture_of_proteomes_processing*.R` and `Mixture_of_proteomes_analysis.R`
- `benchmark_experiments_functions.R` (top-level)

**1.2 Create `quant/R/` with these modules:**
- `preprocess.R` — Spectronaut TSV → tool-specific inputs; reads annotation; filters by dataset (`all_data`/`good_data`/`small`/`small_good_data`).
- `normalize.R` — `log2`, `median`, `quantile`; uniform interface per package; document each tool's actual normalization behaviour (cf. §0.5).
- `dataset_subsets.R` — `<good_data>` and `<small_good_data>` rules; per-folder.
- `models_*.R` — one file per package (`models_msstats.R`, `models_msstatsplus.R`, `models_maxlfq_limma.R`, `models_deqms.R`, `models_msqrob2.R`, `models_prolfqua.R`, `models_limpa.R`). Each exports a single function `run_<pkg>(input_dir, out_dir, normalization)`.
- `ground_truth.R` — TP/TN lookup; protein-swap version (from `CSF_protein_swap_list.csv`) and sample-swap version (from the new ground-truth TSV emitted by §2.2).
- `comparison_table.R` — generalised Table 1 (counts of TP/FP/TN/FN at p-value and FDR thresholds, across packages/normalizations).
- `figures.R` — heatmaps, NA heatmaps, density plots, p-value histograms, effect-size vs variance scatter, FDR curves.

**1.3 Create `quant/vignettes/`:**
- `diagnostics.qmd` — generalised version of `all_dilutions/V1_log_diagnostics.qmd`, parameterised by `<csffolder>` and `<dataset>`.
- `swap_visualization.qmd` — generalised version of `CSF_swap_visualization.qmd`, parameterised by swap type (sample/protein) and folder.
- `review.qmd` — see §7.

**1.4 Top-level driver script** `quant/run_all.fish` that, for each (folder, dataset, normalization, package) combination defined in a manifest CSV, invokes the per-package `run_<pkg>()` and writes outputs into the §0.1 path.

---

## 2. New TSV-level sample-swap script

**2.1** Read [swap_spectronaut_report.py](../swap_spectronaut_report.py) and `CSF_swap_design.md` to mirror its CLI surface and output schema.

**2.2** Write `quant/src/swap_spectronaut_report_samples.py` (create `quant/src/` if absent; this is the home for Python helpers — also relocate `quant/swap_spectronaut_report.py` into `quant/src/` and update any callers/`pyproject.toml` references):
- Inputs: original Spectronaut TSV, annotation CSV, swap fraction (default 0.9 — matches the description "swapped samples for 90% of the proteins"), seed.
- For each protein selected for swap: swap intensities between G1/G2 runs in a way that matches the label-level swap currently done in `CSF_Spectronaut_preswap_variants.R::swap_condition_labels_msstats`.
- Outputs (next to the TSV):
  - `*_Report_sample_swap.tsv` (schema-identical Spectronaut TSV with swapped intensities).
  - `*_Report_sample_swap_ground_truth.tsv` (TP/TN labels per protein).
  - `*_Report_sample_swap_true_positives.tsv` and `*_group_annotation.csv` mirroring the protein-swap outputs.

**2.3** Smoke test: rerun MSstats label-swap on the original TSV and MSstats no-swap on the new TSV; verify protein-level results match within numerical tolerance.

---

## 3. `CSF_Spectronaut` — author-faithful replication, output restructured

**3.1** Do **not** rewrite the R scripts. Keep `CSF_Spectronaut_analysis.R`, `CSF_Spectronaut_processing.R`, `CSF_Spectronaut_preswap_variants.R`, `CSF_Spectronaut_model_nonmsstats_variant.R` as they are — only adapt their output-path arguments.

**3.2** Add a thin wrapper `CSF_Spectronaut/run_replication.fish` that, for each `(dataset ∈ {all_data, good_data}, normalization ∈ {log2, median, quantile}, package)`, writes results into:
```
CSF_Spectronaut/<dataset>/<normalization>/swap/<package>/
```

**3.3** Datasets:
- `all_data`: full annotation (current `all_dilutions`).
- `good_data`: neat + 1/2 only, balanced replicate count.

**3.4** Run all 7 packages × 3 normalizations × 2 datasets. Archive the existing `all_dilutions/V*/` results under `CSF_Spectronaut/_legacy/` before regenerating.

**3.5** Sanity check: compare the new `all_data/log2/swap/MSstats/` numbers against the existing legacy `all_dilutions/V1_log2/.../MSstats/` to confirm we have not regressed.

---

## 4. `CSF_Spectronaut_protein_swap` — rename + restructure

**4.1** Rename folder: `CSF_Spectronaut_swap` → `CSF_Spectronaut_protein_swap`. Update all internal references.

**4.2** Disable the `noswap` branch (currently visible in scripts as the non-`_preswap_` model variants). Comment out / gate behind a flag — re-enable later.

**4.3** Datasets:
- `all_data`: keep (was `all_dilutions`).
- `small`: keep.
- `good_data`: **neat-only** subset (TP/TN come from protein swap, dilution is not the ground-truth axis).
- `small_good_data`: 4 vs 4 sampled from `good_data` with a fixed seed.
- Drop `no_high_dilutions` entirely.

**4.4** Replace per-folder R scripts with calls into `quant/R/`. The Python protein-swap TSV is the input; ground truth is `CSF_protein_swap_list.csv`.

**4.5** Run 7 packages × 3 normalizations × 4 datasets, output under §0.1.

---

## 5. `CSF_Spectronaut_sample_swap` — new folder

**5.1** Create `CSF_Spectronaut_sample_swap/`. Run §2 to produce the swapped Spectronaut TSV + ground-truth files inside it.

**5.2** Datasets:
- `all_data`
- `good_data`: neat + 1/2 balanced.
- `small_good_data`: 4 vs 4 from `good_data`.

**5.3** Run the same `quant/R/` pipeline as protein_swap, swap=sample. 7 packages × 3 normalizations × 3 datasets.

**5.4** Cross-check: for MSstats `all_data/log2/swap`, compare results here vs the label-level swap inside `CSF_Spectronaut/all_data/log2/swap/MSstats/`. They should agree closely; large discrepancies indicate a bug in the new sample-swap script.

---

## 6. `Mix_of_Proteome` — bring under the shared pipeline

**6.1** Ground truth for Mix_of_Proteome comes from species labels, not from a swap. Add a `ground_truth.R::mix_of_proteome_truth()` that reads `idmapping.tsv` / species annotation.

**6.2** Create `Mix_of_Proteome/all_data/<normalization>/swap/<package>/` outputs (the `swap` literal is retained for path uniformity even though no swap is performed; document this).

**6.3** Run 7 packages × 3 normalizations × 1 dataset. Archive existing `V1_log2/`, `DEqMS/`, `limma/`, `MSstats/`, `MSstats+/` outputs under `_legacy/` first.

---

## 7. The review document — `quant/vignettes/review.qmd`

Structure mirrors `WhatWeStillNeed.md` end-section. The qmd reads aggregated CSVs from `quant/results/` (built by `R/comparison_table.R` and `R/figures.R`).

**7.1 Introduction**
- Olga Vitek and her contribution; MSstats history (first published 2014); citations where competing tools outperformed MSstats: msqrob2, limpa, prolfqua.
- Summary of what MSstatsPlus claims.

**7.2 Benchmark Table 1 — biorxiv vs resubmission**
- Reproduce the two Table 1 variants from the manuscript (p-value vs FDR thresholding), with explicit threshold values and provenance.
- Describe how CSF and K562 benchmarks were constructed.

**7.3 Our reproduction (CSF_Spectronaut)**
- Table 1 from our run, p-value and FDR, `all_data` and `good_data`. Note numerical drift due to msqrob2 fast path.
- Contrast with `Mix_of_Proteome` `all_data` Table 1. Expect Mix-of-Proteomes to look uniform across tools; CSF `good_data` to look different — that asymmetry motivates the rest of the review.
- p-value distribution histograms for `good_data` (uses `diagnostics.qmd`). Discuss H0 uniformity expectation and what non-uniformity implies.

**7.4 The Sample Swap benchmark dataset**
- Figure 1 from `CSF_swap_visualization.qmd` (swap schema).
- Visualisations from `swap_visualization.qmd` parameterised on `CSF_Spectronaut_sample_swap`: heatmap (`all_data`), NA heatmap (`all_data`), density plots (`good_data`).
- Per-protein effect size and within-condition variance figures (from same qmd).
- Discussion of variance differences H0 vs H1 and how prior-N drives shrinkage.
- Table 1 for `good_data` and `small_good_data`, contrasting sample-size effects.

**7.5 The Protein Swap benchmark dataset**
- Schematic figure of the protein swap (to be drawn — see §8.4).
- Same visualisation suite on `CSF_Spectronaut_protein_swap` (`all_data` heatmaps; `good_data` density; effect-size/variance scatter).

**7.6 Benchmark results — protein swap**
- Compare `good_data` (= neat only) on protein_swap vs `all_data` on Mix_of_Proteome — expect similar performance profiles across tools.
- Contrast with sample_swap `good_data` to highlight the variance-shrinkage failure mode.
- Show `all_data` and `small` for protein_swap, contrasted with sample_swap counterparts.

**7.7 Conclusions** — placeholder.

---

## 8. Deliverables — Zenodo & GitHub

**8.1 Authors' upload (Zenodo only).** Package the original `RMSV000000701.3-rerun/quant/` as received from the authors. Link to the MSstatsPlus biorxiv preprint. README states purpose: "Snapshot of the authors' submitted scripts; no modifications."

**8.2 Our work — GitHub repo.** Repo contents:
- `CSF_Spectronaut/`, `CSF_Spectronaut_protein_swap/`, `CSF_Spectronaut_sample_swap/`, `Mix_of_Proteome/`
- `R/`, `vignettes/`, `results/`
- `src/` (contains `swap_spectronaut_report.py`, `swap_spectronaut_report_samples.py`), `benchmark_experiments_functions.R`
- Top-level `README.md`, `LICENSE`, `pyproject.toml`
- Exclude: `CSF_DIANN`, `K562_*`, `TPCPA`, `__pycache__`, `_legacy/`

**8.3** Use the existing git repo in `quant/`. Push to GitHub. Connect to Zenodo via the GitHub integration; tag a release; reference the Zenodo DOI from the repo README and from `review.qmd`.

**8.4** Draw the protein-swap schematic (§7.5) — TikZ/SVG; commit under `vignettes/figures/`.

---

## 9. Execution order (suggested)

1. §0.5 — audit median normalization across tools (small, do first; results inform §1.2 `normalize.R`).
2. §1 — build `quant/R/` and `quant/vignettes/` scaffolding.
3. §2 — write and validate sample-swap script.
4. §4 — rename + restructure protein_swap (most code already exists, easiest port).
5. §6 — port Mix_of_Proteome onto shared pipeline.
6. §5 — sample_swap end-to-end run.
7. §3 — re-run CSF_Spectronaut with restructured outputs.
8. §7 — write review.qmd against the now-populated `results/`.
9. §8 — Zenodo + GitHub publication.

## 10. Open items to revisit at the end

- Re-enable the `noswap` branch in `CSF_Spectronaut_protein_swap` (§4.2).
- Decide whether to also include sample_swap `small` (currently only `all_data`, `good_data`, `small_good_data` are planned).
- Confirm 4v4 stratification rule for `small_good_data` (random with fixed seed vs first-4-by-replicate-id).
