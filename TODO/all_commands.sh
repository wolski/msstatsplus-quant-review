#!/bin/sh
# All commands to run end-to-end. cwd = quant/.
# Every command spelled out — no loops, no variables. Paths relative to quant/.

# =============================================================================
# 0. CANONICAL RAW INPUTS (one-time; relative symlinks, originals preserved)
# =============================================================================
# Only the two raw experiment sources need symlinks. Everything else under
# CSF_Spectronaut_protein_swap/ and CSF_Spectronaut_sample_swap/ is generated
# by the Python swap scripts in §2 / §3 below — those scripts write
# Report.tsv and annotation.csv directly.
ln -s "20250130_163144_CSF dilutions Jan 2025 no normalization_Report.tsv"        CSF_Spectronaut/Report.tsv
ln -s  CSF_annotation.csv                                                          CSF_Spectronaut/annotation.csv
ln -s "20250422_140629_OKE_April22_2025_Olsen_astral_benchmark_mix_200ng_NE_report_nodecoy.tsv" Mix_of_Proteome/Report.tsv
ln -s  Mix_of_Proteome_annotation_contrast.csv                                     Mix_of_Proteome/annotation.csv


# =============================================================================
# 1. CSF_Spectronaut  (authors' replication; log2 only)
# =============================================================================

# --- Data prep: build the two subset dirs (no swap) ---
python src/build_subsets.py --report CSF_Spectronaut/Report.tsv --annotation CSF_Spectronaut/annotation.csv --out-dir CSF_Spectronaut --subsets all_data good_data --good-rule label_good

# --- Modelling MSstats / MSstats+ (4 commands) ---
Rscript R/run_cell.R CSF_Spectronaut/all_data  log2 MSstats
Rscript R/run_cell.R CSF_Spectronaut/all_data  log2 "MSstats+"
Rscript R/run_cell.R CSF_Spectronaut/good_data log2 MSstats
Rscript R/run_cell.R CSF_Spectronaut/good_data log2 "MSstats+"

# --- Modelling non-MSstats bundles (2 commands, 10 cells) ---
Rscript R/run_nonmsstats_block.R CSF_Spectronaut/all_data  log2
Rscript R/run_nonmsstats_block.R CSF_Spectronaut/good_data log2

# --- Diagnostics (2 HTMLs) ---
quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut/all_data/log2/swap  -P subset_dir=CSF_Spectronaut/all_data  -P normalization=log2 -P truth_kind=sample_swap -P truth_path=CSF_Spectronaut/CSF_protein_swap_list.csv -P base_dir=..
quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut/good_data/log2/swap -P subset_dir=CSF_Spectronaut/good_data -P normalization=log2 -P truth_kind=sample_swap -P truth_path=CSF_Spectronaut/CSF_protein_swap_list.csv -P base_dir=..


# =============================================================================
# 2. CSF_Spectronaut_sample_swap  (NEW: 3041-protein universe, 10/90)
# =============================================================================

# --- Data prep: regenerate swap list, run TSV-level swap, build subsets ---
Rscript R/build_sample_swap_list.R
python src/swap_spectronaut_report_samples.py --report CSF_Spectronaut/Report.tsv --annotation CSF_Spectronaut/annotation.csv --negatives CSF_Spectronaut_sample_swap/CSF_protein_swap_list.csv --out-dir CSF_Spectronaut_sample_swap
python src/build_subsets.py --report CSF_Spectronaut_sample_swap/Report.tsv --annotation CSF_Spectronaut_sample_swap/annotation.csv --out-dir CSF_Spectronaut_sample_swap --subsets all_data good_data small_good_data --good-rule label_good

# --- Modelling MSstats / MSstats+ (18 commands) ---
Rscript R/run_cell.R CSF_Spectronaut_sample_swap/all_data        log2     MSstats
Rscript R/run_cell.R CSF_Spectronaut_sample_swap/all_data        log2     "MSstats+"
Rscript R/run_cell.R CSF_Spectronaut_sample_swap/all_data        median   MSstats
Rscript R/run_cell.R CSF_Spectronaut_sample_swap/all_data        median   "MSstats+"
Rscript R/run_cell.R CSF_Spectronaut_sample_swap/all_data        quantile MSstats
Rscript R/run_cell.R CSF_Spectronaut_sample_swap/all_data        quantile "MSstats+"
Rscript R/run_cell.R CSF_Spectronaut_sample_swap/good_data       log2     MSstats
Rscript R/run_cell.R CSF_Spectronaut_sample_swap/good_data       log2     "MSstats+"
Rscript R/run_cell.R CSF_Spectronaut_sample_swap/good_data       median   MSstats
Rscript R/run_cell.R CSF_Spectronaut_sample_swap/good_data       median   "MSstats+"
Rscript R/run_cell.R CSF_Spectronaut_sample_swap/good_data       quantile MSstats
Rscript R/run_cell.R CSF_Spectronaut_sample_swap/good_data       quantile "MSstats+"
Rscript R/run_cell.R CSF_Spectronaut_sample_swap/small_good_data log2     MSstats
Rscript R/run_cell.R CSF_Spectronaut_sample_swap/small_good_data log2     "MSstats+"
Rscript R/run_cell.R CSF_Spectronaut_sample_swap/small_good_data median   MSstats
Rscript R/run_cell.R CSF_Spectronaut_sample_swap/small_good_data median   "MSstats+"
Rscript R/run_cell.R CSF_Spectronaut_sample_swap/small_good_data quantile MSstats
Rscript R/run_cell.R CSF_Spectronaut_sample_swap/small_good_data quantile "MSstats+"

# --- Modelling non-MSstats bundles (9 commands, 45 cells) ---
Rscript R/run_nonmsstats_block.R CSF_Spectronaut_sample_swap/all_data        log2
Rscript R/run_nonmsstats_block.R CSF_Spectronaut_sample_swap/all_data        median
Rscript R/run_nonmsstats_block.R CSF_Spectronaut_sample_swap/all_data        quantile
Rscript R/run_nonmsstats_block.R CSF_Spectronaut_sample_swap/good_data       log2
Rscript R/run_nonmsstats_block.R CSF_Spectronaut_sample_swap/good_data       median
Rscript R/run_nonmsstats_block.R CSF_Spectronaut_sample_swap/good_data       quantile
Rscript R/run_nonmsstats_block.R CSF_Spectronaut_sample_swap/small_good_data log2
Rscript R/run_nonmsstats_block.R CSF_Spectronaut_sample_swap/small_good_data median
Rscript R/run_nonmsstats_block.R CSF_Spectronaut_sample_swap/small_good_data quantile

# --- Diagnostics (9 HTMLs) ---
quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_sample_swap/all_data/log2/swap            -P subset_dir=CSF_Spectronaut_sample_swap/all_data        -P normalization=log2     -P truth_kind=sample_swap -P truth_path=CSF_Spectronaut_sample_swap/CSF_protein_swap_list.csv -P base_dir=..
quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_sample_swap/all_data/median/swap          -P subset_dir=CSF_Spectronaut_sample_swap/all_data        -P normalization=median   -P truth_kind=sample_swap -P truth_path=CSF_Spectronaut_sample_swap/CSF_protein_swap_list.csv -P base_dir=..
quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_sample_swap/all_data/quantile/swap        -P subset_dir=CSF_Spectronaut_sample_swap/all_data        -P normalization=quantile -P truth_kind=sample_swap -P truth_path=CSF_Spectronaut_sample_swap/CSF_protein_swap_list.csv -P base_dir=..
quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_sample_swap/good_data/log2/swap           -P subset_dir=CSF_Spectronaut_sample_swap/good_data       -P normalization=log2     -P truth_kind=sample_swap -P truth_path=CSF_Spectronaut_sample_swap/CSF_protein_swap_list.csv -P base_dir=..
quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_sample_swap/good_data/median/swap         -P subset_dir=CSF_Spectronaut_sample_swap/good_data       -P normalization=median   -P truth_kind=sample_swap -P truth_path=CSF_Spectronaut_sample_swap/CSF_protein_swap_list.csv -P base_dir=..
quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_sample_swap/good_data/quantile/swap       -P subset_dir=CSF_Spectronaut_sample_swap/good_data       -P normalization=quantile -P truth_kind=sample_swap -P truth_path=CSF_Spectronaut_sample_swap/CSF_protein_swap_list.csv -P base_dir=..
quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_sample_swap/small_good_data/log2/swap     -P subset_dir=CSF_Spectronaut_sample_swap/small_good_data -P normalization=log2     -P truth_kind=sample_swap -P truth_path=CSF_Spectronaut_sample_swap/CSF_protein_swap_list.csv -P base_dir=..
quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_sample_swap/small_good_data/median/swap   -P subset_dir=CSF_Spectronaut_sample_swap/small_good_data -P normalization=median   -P truth_kind=sample_swap -P truth_path=CSF_Spectronaut_sample_swap/CSF_protein_swap_list.csv -P base_dir=..
quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_sample_swap/small_good_data/quantile/swap -P subset_dir=CSF_Spectronaut_sample_swap/small_good_data -P normalization=quantile -P truth_kind=sample_swap -P truth_path=CSF_Spectronaut_sample_swap/CSF_protein_swap_list.csv -P base_dir=..


# =============================================================================
# 3. CSF_Spectronaut_protein_swap
# =============================================================================

# --- Data prep: protein-swap TSV, then build subsets ---
python src/swap_spectronaut_report.py --report CSF_Spectronaut/Report.tsv --annotation CSF_Spectronaut/annotation.csv --out-dir CSF_Spectronaut_protein_swap --swap-fraction 0.05 --target-log2fc-min 1.4 --target-log2fc-max 1.8 --seed 42
python src/build_subsets.py --report CSF_Spectronaut_protein_swap/Report.tsv --annotation CSF_Spectronaut_protein_swap/annotation.csv --out-dir CSF_Spectronaut_protein_swap --subsets all_data good_data small small_good_data --good-rule neat_only

# --- Modelling MSstats / MSstats+ (24 commands) ---
Rscript R/run_cell.R CSF_Spectronaut_protein_swap/all_data        log2     MSstats
Rscript R/run_cell.R CSF_Spectronaut_protein_swap/all_data        log2     "MSstats+"
Rscript R/run_cell.R CSF_Spectronaut_protein_swap/all_data        median   MSstats
Rscript R/run_cell.R CSF_Spectronaut_protein_swap/all_data        median   "MSstats+"
Rscript R/run_cell.R CSF_Spectronaut_protein_swap/all_data        quantile MSstats
Rscript R/run_cell.R CSF_Spectronaut_protein_swap/all_data        quantile "MSstats+"
Rscript R/run_cell.R CSF_Spectronaut_protein_swap/good_data       log2     MSstats
Rscript R/run_cell.R CSF_Spectronaut_protein_swap/good_data       log2     "MSstats+"
Rscript R/run_cell.R CSF_Spectronaut_protein_swap/good_data       median   MSstats
Rscript R/run_cell.R CSF_Spectronaut_protein_swap/good_data       median   "MSstats+"
Rscript R/run_cell.R CSF_Spectronaut_protein_swap/good_data       quantile MSstats
Rscript R/run_cell.R CSF_Spectronaut_protein_swap/good_data       quantile "MSstats+"
Rscript R/run_cell.R CSF_Spectronaut_protein_swap/small           log2     MSstats
Rscript R/run_cell.R CSF_Spectronaut_protein_swap/small           log2     "MSstats+"
Rscript R/run_cell.R CSF_Spectronaut_protein_swap/small           median   MSstats
Rscript R/run_cell.R CSF_Spectronaut_protein_swap/small           median   "MSstats+"
Rscript R/run_cell.R CSF_Spectronaut_protein_swap/small           quantile MSstats
Rscript R/run_cell.R CSF_Spectronaut_protein_swap/small           quantile "MSstats+"
Rscript R/run_cell.R CSF_Spectronaut_protein_swap/small_good_data log2     MSstats
Rscript R/run_cell.R CSF_Spectronaut_protein_swap/small_good_data log2     "MSstats+"
Rscript R/run_cell.R CSF_Spectronaut_protein_swap/small_good_data median   MSstats
Rscript R/run_cell.R CSF_Spectronaut_protein_swap/small_good_data median   "MSstats+"
Rscript R/run_cell.R CSF_Spectronaut_protein_swap/small_good_data quantile MSstats
Rscript R/run_cell.R CSF_Spectronaut_protein_swap/small_good_data quantile "MSstats+"

# --- Modelling non-MSstats bundles (12 commands, 60 cells) ---
Rscript R/run_nonmsstats_block.R CSF_Spectronaut_protein_swap/all_data        log2
Rscript R/run_nonmsstats_block.R CSF_Spectronaut_protein_swap/all_data        median
Rscript R/run_nonmsstats_block.R CSF_Spectronaut_protein_swap/all_data        quantile
Rscript R/run_nonmsstats_block.R CSF_Spectronaut_protein_swap/good_data       log2
Rscript R/run_nonmsstats_block.R CSF_Spectronaut_protein_swap/good_data       median
Rscript R/run_nonmsstats_block.R CSF_Spectronaut_protein_swap/good_data       quantile
Rscript R/run_nonmsstats_block.R CSF_Spectronaut_protein_swap/small           log2
Rscript R/run_nonmsstats_block.R CSF_Spectronaut_protein_swap/small           median
Rscript R/run_nonmsstats_block.R CSF_Spectronaut_protein_swap/small           quantile
Rscript R/run_nonmsstats_block.R CSF_Spectronaut_protein_swap/small_good_data log2
Rscript R/run_nonmsstats_block.R CSF_Spectronaut_protein_swap/small_good_data median
Rscript R/run_nonmsstats_block.R CSF_Spectronaut_protein_swap/small_good_data quantile

# --- Diagnostics (12 HTMLs) ---
quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_protein_swap/all_data/log2/swap            -P subset_dir=CSF_Spectronaut_protein_swap/all_data        -P normalization=log2     -P truth_kind=protein_swap -P truth_path=CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv -P base_dir=..
quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_protein_swap/all_data/median/swap          -P subset_dir=CSF_Spectronaut_protein_swap/all_data        -P normalization=median   -P truth_kind=protein_swap -P truth_path=CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv -P base_dir=..
quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_protein_swap/all_data/quantile/swap        -P subset_dir=CSF_Spectronaut_protein_swap/all_data        -P normalization=quantile -P truth_kind=protein_swap -P truth_path=CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv -P base_dir=..
quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_protein_swap/good_data/log2/swap           -P subset_dir=CSF_Spectronaut_protein_swap/good_data       -P normalization=log2     -P truth_kind=protein_swap -P truth_path=CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv -P base_dir=..
quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_protein_swap/good_data/median/swap         -P subset_dir=CSF_Spectronaut_protein_swap/good_data       -P normalization=median   -P truth_kind=protein_swap -P truth_path=CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv -P base_dir=..
quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_protein_swap/good_data/quantile/swap       -P subset_dir=CSF_Spectronaut_protein_swap/good_data       -P normalization=quantile -P truth_kind=protein_swap -P truth_path=CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv -P base_dir=..
quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_protein_swap/small/log2/swap               -P subset_dir=CSF_Spectronaut_protein_swap/small           -P normalization=log2     -P truth_kind=protein_swap -P truth_path=CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv -P base_dir=..
quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_protein_swap/small/median/swap             -P subset_dir=CSF_Spectronaut_protein_swap/small           -P normalization=median   -P truth_kind=protein_swap -P truth_path=CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv -P base_dir=..
quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_protein_swap/small/quantile/swap           -P subset_dir=CSF_Spectronaut_protein_swap/small           -P normalization=quantile -P truth_kind=protein_swap -P truth_path=CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv -P base_dir=..
quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_protein_swap/small_good_data/log2/swap     -P subset_dir=CSF_Spectronaut_protein_swap/small_good_data -P normalization=log2     -P truth_kind=protein_swap -P truth_path=CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv -P base_dir=..
quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_protein_swap/small_good_data/median/swap   -P subset_dir=CSF_Spectronaut_protein_swap/small_good_data -P normalization=median   -P truth_kind=protein_swap -P truth_path=CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv -P base_dir=..
quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_protein_swap/small_good_data/quantile/swap -P subset_dir=CSF_Spectronaut_protein_swap/small_good_data -P normalization=quantile -P truth_kind=protein_swap -P truth_path=CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv -P base_dir=..


# =============================================================================
# 4. Mix_of_Proteome
# =============================================================================

# --- Data prep ---
python src/build_subsets.py --report Mix_of_Proteome/Report.tsv --annotation Mix_of_Proteome/annotation.csv --out-dir Mix_of_Proteome --subsets all_data

# --- Modelling MSstats / MSstats+ (6 commands) ---
Rscript R/run_cell.R Mix_of_Proteome/all_data log2     MSstats
Rscript R/run_cell.R Mix_of_Proteome/all_data log2     "MSstats+"
Rscript R/run_cell.R Mix_of_Proteome/all_data median   MSstats
Rscript R/run_cell.R Mix_of_Proteome/all_data median   "MSstats+"
Rscript R/run_cell.R Mix_of_Proteome/all_data quantile MSstats
Rscript R/run_cell.R Mix_of_Proteome/all_data quantile "MSstats+"

# --- Modelling non-MSstats bundles (3 commands, 15 cells) ---
Rscript R/run_nonmsstats_block.R Mix_of_Proteome/all_data log2
Rscript R/run_nonmsstats_block.R Mix_of_Proteome/all_data median
Rscript R/run_nonmsstats_block.R Mix_of_Proteome/all_data quantile

# --- Diagnostics (3 HTMLs) ---
# NOTE: replace idmapping_TODO.tsv with the actual filename in Mix_of_Proteome/
quarto render vignettes/diagnostics.qmd --output-dir ../Mix_of_Proteome/all_data/log2/swap     -P subset_dir=Mix_of_Proteome/all_data -P normalization=log2     -P truth_kind=mix_of_proteome -P truth_path=Mix_of_Proteome/idmapping_TODO.tsv -P base_dir=..
quarto render vignettes/diagnostics.qmd --output-dir ../Mix_of_Proteome/all_data/median/swap   -P subset_dir=Mix_of_Proteome/all_data -P normalization=median   -P truth_kind=mix_of_proteome -P truth_path=Mix_of_Proteome/idmapping_TODO.tsv -P base_dir=..
quarto render vignettes/diagnostics.qmd --output-dir ../Mix_of_Proteome/all_data/quantile/swap -P subset_dir=Mix_of_Proteome/all_data -P normalization=quantile -P truth_kind=mix_of_proteome -P truth_path=Mix_of_Proteome/idmapping_TODO.tsv -P base_dir=..


# =============================================================================
# 5. Per-folder swap visualizations  (OPTIONAL — review.qmd already embeds these)
# =============================================================================
quarto render vignettes/swap_visualization.qmd --output-dir ../CSF_Spectronaut_sample_swap/all_data   -P subset_dir=CSF_Spectronaut_sample_swap/all_data  -P normalization=log2 -P base_dir=..
quarto render vignettes/swap_visualization.qmd --output-dir ../CSF_Spectronaut_sample_swap/good_data  -P subset_dir=CSF_Spectronaut_sample_swap/good_data -P normalization=log2 -P base_dir=..
quarto render vignettes/swap_visualization.qmd --output-dir ../CSF_Spectronaut_protein_swap/all_data  -P subset_dir=CSF_Spectronaut_protein_swap/all_data  -P normalization=log2 -P base_dir=..
quarto render vignettes/swap_visualization.qmd --output-dir ../CSF_Spectronaut_protein_swap/good_data -P subset_dir=CSF_Spectronaut_protein_swap/good_data -P normalization=log2 -P base_dir=..


# =============================================================================
# 6. Review
# =============================================================================
( cd vignettes && quarto render review.qmd --to html )
( cd vignettes && quarto render review.qmd --to pdf )
