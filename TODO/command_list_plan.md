# Plan — Explicit command list per target (Makefile design comes after)

## Context

User direction: drop the auto-generated cells.mk idea (premature
abstraction). First enumerate every concrete shell command — literal
paths, no variables — for each target. Once the command list is on
paper, decide what becomes a Make rule and what does not.

All commands assume `cwd = quant/`. Datasets per folder:

| Folder | Datasets | Norms |
|---|---|---|
| CSF_Spectronaut (replicate authors) | all_data, good_data | log2 only |
| CSF_Spectronaut_sample_swap | all_data, good_data, small_good_data | log2, median, quantile |
| CSF_Spectronaut_protein_swap | all_data, good_data, small, small_good_data | log2, median, quantile |
| Mix_of_Proteome | all_data | log2, median, quantile |

Packages everywhere: `MSstats`, `MSstats+`, `MaxLFQ_limma`, `DEqMS`,
`msqrob2`, `prolfqua`, `limpa`.

Pkg → on-disk filename (from `R/comparison_table.R::.model_filename`):
`MSstats → MSstats_model.csv`, `MSstats+ → MSstats+_model.csv`,
`MaxLFQ_limma → limma_model.csv`, `DEqMS → deqms_model.csv`,
`msqrob2 → msqrob2_model.csv`, `prolfqua → prolfqua_model.csv`,
`limpa → limpa_model.csv`.

---

## Explicit command list

### 0. Manifest

```sh
Rscript R/build_manifest.R
# writes manifest.csv (after build_manifest.R is updated to CSF_Spectronaut=log2-only)
```

### 1. CSF_Spectronaut (authors' replication)

**Data prep**: none — author's raw TSV is the input.
Truth file: `CSF_Spectronaut/CSF_protein_swap_list.csv` (1,820 list,
existing).

**Modelling** — 2 datasets × 1 norm × 7 packages = **14 cells**.
MSstats run individually; the 5 non-MSstats bundled per (dataset, norm).

```sh
# MSstats / MSstats+ (one Rscript per cell)
Rscript R/run_cell.R CSF_Spectronaut all_data  log2 MSstats
Rscript R/run_cell.R CSF_Spectronaut all_data  log2 "MSstats+"
Rscript R/run_cell.R CSF_Spectronaut good_data log2 MSstats
Rscript R/run_cell.R CSF_Spectronaut good_data log2 "MSstats+"

# Non-MSstats bundle (one Rscript per (dataset, norm) — produces 5 model CSVs)
Rscript R/run_nonmsstats_block.R CSF_Spectronaut all_data  log2
Rscript R/run_nonmsstats_block.R CSF_Spectronaut good_data log2
```

**Diagnostics** — 2 cells:

```sh
quarto render vignettes/diagnostics.qmd \
  --output-dir ../CSF_Spectronaut/all_data/log2/swap \
  -P csffolder=CSF_Spectronaut -P dataset=all_data -P normalization=log2 \
  -P truth_kind=sample_swap \
  -P truth_path=CSF_Spectronaut/CSF_protein_swap_list.csv \
  -P base_dir=..

quarto render vignettes/diagnostics.qmd \
  --output-dir ../CSF_Spectronaut/good_data/log2/swap \
  -P csffolder=CSF_Spectronaut -P dataset=good_data -P normalization=log2 \
  -P truth_kind=sample_swap \
  -P truth_path=CSF_Spectronaut/CSF_protein_swap_list.csv \
  -P base_dir=..
```

### 2. CSF_Spectronaut_sample_swap

**Data prep**:

```sh
# (a) Regenerate the sample-swap ground-truth list on the 3,041-protein universe
Rscript R/build_sample_swap_list.R
#   reads  CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv (3,041 proteins)
#   writes CSF_Spectronaut_sample_swap/CSF_protein_swap_list.csv (~304 Pos / ~2737 Neg)

# (b) Rebuild the swapped TSV from the canonical Spectronaut report
#     (only needed if Negative protein set changed — which it just did)
python src/swap_spectronaut_report_samples.py \
  --report CSF_Spectronaut/'20250130_163144_CSF dilutions Jan 2025 no normalization_Report.tsv' \
  --annotation CSF_Spectronaut/CSF_annotation.csv \
  --negatives CSF_Spectronaut_sample_swap/CSF_protein_swap_list.csv \
  --out-dir CSF_Spectronaut_sample_swap
#   writes:
#     CSF_Spectronaut_sample_swap/...Report_sample_swap.tsv
#     ...Report_sample_swap_ground_truth.tsv
#     ...Report_sample_swap_group_annotation.csv
#     ...Report_sample_swap_true_positives.tsv
#
# (verify the actual argparse interface of swap_spectronaut_report_samples.py
# before running — may need a --negatives or --truth-list argument added)
```

**Modelling** — 3 datasets × 3 norms × 7 packages = **63 cells**.

```sh
# MSstats / MSstats+
for d in all_data good_data small_good_data; do
  for n in log2 median quantile; do
    Rscript R/run_cell.R CSF_Spectronaut_sample_swap $d $n MSstats
    Rscript R/run_cell.R CSF_Spectronaut_sample_swap $d $n "MSstats+"
  done
done

# Non-MSstats bundles
for d in all_data good_data small_good_data; do
  for n in log2 median quantile; do
    Rscript R/run_nonmsstats_block.R CSF_Spectronaut_sample_swap $d $n
  done
done
```

**Diagnostics** — 9 cells (3 datasets × 3 norms):

```sh
for d in all_data good_data small_good_data; do
  for n in log2 median quantile; do
    quarto render vignettes/diagnostics.qmd \
      --output-dir ../CSF_Spectronaut_sample_swap/$d/$n/swap \
      -P csffolder=CSF_Spectronaut_sample_swap -P dataset=$d -P normalization=$n \
      -P truth_kind=sample_swap \
      -P truth_path=CSF_Spectronaut_sample_swap/CSF_protein_swap_list.csv \
      -P base_dir=..
  done
done
```

### 3. CSF_Spectronaut_protein_swap

**Data prep** (already done — kept for completeness):

```sh
# (a) Generate the protein-swap TSV + 3,041-protein ground truth list
python src/swap_spectronaut_report.py \
  --report CSF_Spectronaut/'20250130_163144_CSF dilutions Jan 2025 no normalization_Report.tsv' \
  --annotation CSF_Spectronaut/CSF_annotation.csv \
  --out-dir CSF_Spectronaut_protein_swap \
  --swap-fraction 0.05 --target-log2fc-min 1.4 --target-log2fc-max 1.8 --seed 42
#   writes the protein_swap report TSV + group_annotation.csv + ground_truth.tsv +
#   CSF_protein_swap_list.csv (3041 rows, 224 Pos / 2817 Neg)
```

**Modelling** — 4 datasets × 3 norms × 7 packages = **84 cells**.

```sh
for d in all_data good_data small small_good_data; do
  for n in log2 median quantile; do
    Rscript R/run_cell.R CSF_Spectronaut_protein_swap $d $n MSstats
    Rscript R/run_cell.R CSF_Spectronaut_protein_swap $d $n "MSstats+"
    Rscript R/run_nonmsstats_block.R CSF_Spectronaut_protein_swap $d $n
  done
done
```

**Diagnostics** — 12 cells:

```sh
for d in all_data good_data small small_good_data; do
  for n in log2 median quantile; do
    quarto render vignettes/diagnostics.qmd \
      --output-dir ../CSF_Spectronaut_protein_swap/$d/$n/swap \
      -P csffolder=CSF_Spectronaut_protein_swap -P dataset=$d -P normalization=$n \
      -P truth_kind=protein_swap \
      -P truth_path=CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv \
      -P base_dir=..
  done
done
```

### 4. Mix_of_Proteome

**Data prep**: none — input TSV is the species mix from authors.
Truth file: `Mix_of_Proteome/idmapping_*.tsv` (species → Pos/Neg).

**Modelling** — 1 dataset × 3 norms × 7 packages = **21 cells**.

```sh
for n in log2 median quantile; do
  Rscript R/run_cell.R Mix_of_Proteome all_data $n MSstats
  Rscript R/run_cell.R Mix_of_Proteome all_data $n "MSstats+"
  Rscript R/run_nonmsstats_block.R Mix_of_Proteome all_data $n
done
```

**Diagnostics** — 3 cells:

```sh
for n in log2 median quantile; do
  quarto render vignettes/diagnostics.qmd \
    --output-dir ../Mix_of_Proteome/all_data/$n/swap \
    -P csffolder=Mix_of_Proteome -P dataset=all_data -P normalization=$n \
    -P truth_kind=mix_of_proteome \
    -P truth_path=Mix_of_Proteome/idmapping_*.tsv \
    -P base_dir=..
done
```

### 5. Per-folder swap visualizations

(Heatmap + density figures. Currently embedded inside `review.qmd` — the
parametrized `swap_visualization.qmd` would be an alternative if we
want stand-alone HTMLs per folder × dataset. Not strictly needed if
review.qmd renders them. Decide before adding rules.)

Candidate commands:

```sh
# CSF_Spectronaut_sample_swap (all_data heatmap, good_data variance)
quarto render vignettes/swap_visualization.qmd \
  --output-dir ../CSF_Spectronaut_sample_swap/all_data \
  -P csffolder=CSF_Spectronaut_sample_swap -P dataset=all_data \
  -P normalization=log2 -P base_dir=..

quarto render vignettes/swap_visualization.qmd \
  --output-dir ../CSF_Spectronaut_sample_swap/good_data \
  -P csffolder=CSF_Spectronaut_sample_swap -P dataset=good_data \
  -P normalization=log2 -P base_dir=..

# … same pattern for CSF_Spectronaut_protein_swap (all_data, good_data)
```

### 6. Review

```sh
cd vignettes && quarto render review.qmd --to html
cd vignettes && quarto render review.qmd --to pdf
```

---

## Totals

- Modelling cells: **14 + 63 + 84 + 21 = 182 cells**.
  Of those, MSstats/MSstats+ are 2/cell-block × 26 (dataset, norm)
  blocks = 52 individual Rscript calls; non-MSstats are 26 bundle
  calls. **78 R invocations total** for a full rebuild.
- Diagnostics renders: **2 + 9 + 12 + 3 = 26 HTMLs**.
- Swap visualizations (if enabled): ~4 HTMLs.
- Review: 2 renders (HTML + PDF).

## Concrete edits required before any of the above runs

1. `R/build_manifest.R` — CSF_Spectronaut log2-only.
2. `vignettes/diagnostics.qmd` — default `base_dir: ".."`.
3. `vignettes/swap_visualization.qmd` — default `base_dir: ".."`.
4. `R/run_nonmsstats_block.R` — non-zero exit on any package failure.
5. New `R/build_sample_swap_list.R` — 3,041-universe 10/90 sampler.
6. Verify `src/swap_spectronaut_report_samples.py` accepts a
   `--negatives` (or equivalent) argument so the sample-swap TSV can
   be rebuilt from the new ground-truth list. If not, add it.
7. `vignettes/review.qmd` — §6 heatmap/NA chunks → `all_data`; rewrite
   §6 counts paragraph (lines 820-823) and stale comment (937-939).

## Once command list is approved

Decide which slices become Make rules. Suggested first-pass Makefile
(no auto-generation, just literal targets):

```make
# Top-level phony targets
.PHONY: all manifest cells diagnostics review clean \
        sample_swap_list \
        cells-csf cells-sample-swap cells-protein-swap cells-mix \
        diag-csf diag-sample-swap diag-protein-swap diag-mix

all: cells diagnostics review

manifest:
	Rscript R/build_manifest.R

sample_swap_list:
	Rscript R/build_sample_swap_list.R

cells: cells-csf cells-sample-swap cells-protein-swap cells-mix

cells-csf:
	# (the 6 commands from §1)

cells-sample-swap: sample_swap_list
	# (the 18 commands from §2, expressed as a loop in the recipe)

# ...etc, one explicit recipe per folder.
```

Per-cell dependency tracking can stay out of v1: each `cells-<folder>`
recipe re-runs the loop wholesale and relies on the run scripts to
detect already-current outputs (skip logic in `run_cell.R` /
`run_nonmsstats_block.R`).

## Verification

- Read every command in §§1-6 and confirm paths are exact.
- Spot-run one command per section to validate it executes.
- Then design the Makefile rules.

## Out of scope

- Auto-generated `cells.mk`.
- §8 Zenodo / GitHub publication.
