# CSF_Spectronaut_protein_swap/Makefile — protein-swap benchmark folder.
# 4 subsets × 3 norms × 6 default packages = 72 model fits across 12 cell-blocks.
# Each cell-block has two stamps: .stamp.msstats and .stamp.bundle.
# Downstream depends on .stamp.bundle so MSstats+ failures (Bioc quantile bug)
# do not block diagnostics or review.
# cwd assumed to be quant/. Invoke from quant/ as:
#   make -f CSF_Spectronaut_protein_swap/Makefile <target>

include mk/common.mk

.PHONY: prep-protein-swap cells-protein-swap diag-protein-swap \
        cells-ps-all_data-log2     cells-ps-all_data-median     cells-ps-all_data-quantile \
        cells-ps-good_data-log2    cells-ps-good_data-median    cells-ps-good_data-quantile \
        cells-ps-small-log2        cells-ps-small-median        cells-ps-small-quantile \
        cells-ps-small_good-log2   cells-ps-small_good-median   cells-ps-small_good-quantile \
        diag-ps-all_data-log2      diag-ps-all_data-median      diag-ps-all_data-quantile \
        diag-ps-good_data-log2     diag-ps-good_data-median     diag-ps-good_data-quantile \
        diag-ps-small-log2         diag-ps-small-median         diag-ps-small-quantile \
        diag-ps-small_good-log2    diag-ps-small_good-median    diag-ps-small_good-quantile \
        clean-ps-models clean-ps-subsets clean-ps-prep


# --- Prep --------------------------------------------------------------------
prep-protein-swap: CSF_Spectronaut_protein_swap/.prep.stamp

CSF_Spectronaut_protein_swap/.prep.stamp: CSF_Spectronaut/Report.tsv CSF_Spectronaut/annotation.csv src/swap_spectronaut_report.py src/build_subsets.py
	$(PYTHON) src/swap_spectronaut_report.py --report CSF_Spectronaut/Report.tsv --annotation CSF_Spectronaut/annotation.csv --out-dir CSF_Spectronaut_protein_swap --swap-fraction 0.05 --target-log2fc-min 1.4 --target-log2fc-max 1.8 --seed 42
	$(PYTHON) src/build_subsets.py --report CSF_Spectronaut_protein_swap/Report.tsv --annotation CSF_Spectronaut_protein_swap/annotation.csv --out-dir CSF_Spectronaut_protein_swap --subsets all_data good_data small small_good_data --good-rule neat_only
	touch $@


# --- Cells -------------------------------------------------------------------
cells-protein-swap: \
    cells-ps-all_data-log2     cells-ps-all_data-median     cells-ps-all_data-quantile \
    cells-ps-good_data-log2    cells-ps-good_data-median    cells-ps-good_data-quantile \
    cells-ps-small-log2        cells-ps-small-median        cells-ps-small-quantile \
    cells-ps-small_good-log2   cells-ps-small_good-median   cells-ps-small_good-quantile

cells-ps-all_data-log2:       CSF_Spectronaut_protein_swap/all_data/log2/swap/.stamp.msstats        CSF_Spectronaut_protein_swap/all_data/log2/swap/.stamp.bundle
cells-ps-all_data-median:     CSF_Spectronaut_protein_swap/all_data/median/swap/.stamp.msstats      CSF_Spectronaut_protein_swap/all_data/median/swap/.stamp.bundle
cells-ps-all_data-quantile:   CSF_Spectronaut_protein_swap/all_data/quantile/swap/.stamp.msstats    CSF_Spectronaut_protein_swap/all_data/quantile/swap/.stamp.bundle
cells-ps-good_data-log2:      CSF_Spectronaut_protein_swap/good_data/log2/swap/.stamp.msstats       CSF_Spectronaut_protein_swap/good_data/log2/swap/.stamp.bundle
cells-ps-good_data-median:    CSF_Spectronaut_protein_swap/good_data/median/swap/.stamp.msstats     CSF_Spectronaut_protein_swap/good_data/median/swap/.stamp.bundle
cells-ps-good_data-quantile:  CSF_Spectronaut_protein_swap/good_data/quantile/swap/.stamp.msstats   CSF_Spectronaut_protein_swap/good_data/quantile/swap/.stamp.bundle
cells-ps-small-log2:          CSF_Spectronaut_protein_swap/small/log2/swap/.stamp.msstats           CSF_Spectronaut_protein_swap/small/log2/swap/.stamp.bundle
cells-ps-small-median:        CSF_Spectronaut_protein_swap/small/median/swap/.stamp.msstats         CSF_Spectronaut_protein_swap/small/median/swap/.stamp.bundle
cells-ps-small-quantile:      CSF_Spectronaut_protein_swap/small/quantile/swap/.stamp.msstats       CSF_Spectronaut_protein_swap/small/quantile/swap/.stamp.bundle
cells-ps-small_good-log2:     CSF_Spectronaut_protein_swap/small_good_data/log2/swap/.stamp.msstats     CSF_Spectronaut_protein_swap/small_good_data/log2/swap/.stamp.bundle
cells-ps-small_good-median:   CSF_Spectronaut_protein_swap/small_good_data/median/swap/.stamp.msstats   CSF_Spectronaut_protein_swap/small_good_data/median/swap/.stamp.bundle
cells-ps-small_good-quantile: CSF_Spectronaut_protein_swap/small_good_data/quantile/swap/.stamp.msstats CSF_Spectronaut_protein_swap/small_good_data/quantile/swap/.stamp.bundle

# MSstats stamps (MSstats + MSstats+) ----------------------------------------
CSF_Spectronaut_protein_swap/all_data/log2/swap/.stamp.msstats: CSF_Spectronaut_protein_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_cell.R CSF_Spectronaut_protein_swap/all_data log2 MSstats
	Rscript R/run_cell.R CSF_Spectronaut_protein_swap/all_data log2 "MSstats+"
	touch $@

CSF_Spectronaut_protein_swap/all_data/median/swap/.stamp.msstats: CSF_Spectronaut_protein_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_cell.R CSF_Spectronaut_protein_swap/all_data median MSstats
	Rscript R/run_cell.R CSF_Spectronaut_protein_swap/all_data median "MSstats+"
	touch $@

CSF_Spectronaut_protein_swap/all_data/quantile/swap/.stamp.msstats: CSF_Spectronaut_protein_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_cell.R CSF_Spectronaut_protein_swap/all_data quantile MSstats
	Rscript R/run_cell.R CSF_Spectronaut_protein_swap/all_data quantile "MSstats+"
	touch $@

CSF_Spectronaut_protein_swap/good_data/log2/swap/.stamp.msstats: CSF_Spectronaut_protein_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_cell.R CSF_Spectronaut_protein_swap/good_data log2 MSstats
	Rscript R/run_cell.R CSF_Spectronaut_protein_swap/good_data log2 "MSstats+"
	touch $@

CSF_Spectronaut_protein_swap/good_data/median/swap/.stamp.msstats: CSF_Spectronaut_protein_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_cell.R CSF_Spectronaut_protein_swap/good_data median MSstats
	Rscript R/run_cell.R CSF_Spectronaut_protein_swap/good_data median "MSstats+"
	touch $@

CSF_Spectronaut_protein_swap/good_data/quantile/swap/.stamp.msstats: CSF_Spectronaut_protein_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_cell.R CSF_Spectronaut_protein_swap/good_data quantile MSstats
	Rscript R/run_cell.R CSF_Spectronaut_protein_swap/good_data quantile "MSstats+"
	touch $@

CSF_Spectronaut_protein_swap/small/log2/swap/.stamp.msstats: CSF_Spectronaut_protein_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_cell.R CSF_Spectronaut_protein_swap/small log2 MSstats
	Rscript R/run_cell.R CSF_Spectronaut_protein_swap/small log2 "MSstats+"
	touch $@

CSF_Spectronaut_protein_swap/small/median/swap/.stamp.msstats: CSF_Spectronaut_protein_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_cell.R CSF_Spectronaut_protein_swap/small median MSstats
	Rscript R/run_cell.R CSF_Spectronaut_protein_swap/small median "MSstats+"
	touch $@

CSF_Spectronaut_protein_swap/small/quantile/swap/.stamp.msstats: CSF_Spectronaut_protein_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_cell.R CSF_Spectronaut_protein_swap/small quantile MSstats
	Rscript R/run_cell.R CSF_Spectronaut_protein_swap/small quantile "MSstats+"
	touch $@

CSF_Spectronaut_protein_swap/small_good_data/log2/swap/.stamp.msstats: CSF_Spectronaut_protein_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_cell.R CSF_Spectronaut_protein_swap/small_good_data log2 MSstats
	Rscript R/run_cell.R CSF_Spectronaut_protein_swap/small_good_data log2 "MSstats+"
	touch $@

CSF_Spectronaut_protein_swap/small_good_data/median/swap/.stamp.msstats: CSF_Spectronaut_protein_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_cell.R CSF_Spectronaut_protein_swap/small_good_data median MSstats
	Rscript R/run_cell.R CSF_Spectronaut_protein_swap/small_good_data median "MSstats+"
	touch $@

CSF_Spectronaut_protein_swap/small_good_data/quantile/swap/.stamp.msstats: CSF_Spectronaut_protein_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_cell.R CSF_Spectronaut_protein_swap/small_good_data quantile MSstats
	Rscript R/run_cell.R CSF_Spectronaut_protein_swap/small_good_data quantile "MSstats+"
	touch $@


# Bundle stamps (non-MSstats packages: MaxLFQ_limma, DEqMS, prolfqua, limpa, msqrob2)
CSF_Spectronaut_protein_swap/all_data/log2/swap/.stamp.bundle: CSF_Spectronaut_protein_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_nonmsstats_block.R CSF_Spectronaut_protein_swap/all_data log2
	touch $@

CSF_Spectronaut_protein_swap/all_data/median/swap/.stamp.bundle: CSF_Spectronaut_protein_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_nonmsstats_block.R CSF_Spectronaut_protein_swap/all_data median
	touch $@

CSF_Spectronaut_protein_swap/all_data/quantile/swap/.stamp.bundle: CSF_Spectronaut_protein_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_nonmsstats_block.R CSF_Spectronaut_protein_swap/all_data quantile
	touch $@

CSF_Spectronaut_protein_swap/good_data/log2/swap/.stamp.bundle: CSF_Spectronaut_protein_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_nonmsstats_block.R CSF_Spectronaut_protein_swap/good_data log2
	touch $@

CSF_Spectronaut_protein_swap/good_data/median/swap/.stamp.bundle: CSF_Spectronaut_protein_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_nonmsstats_block.R CSF_Spectronaut_protein_swap/good_data median
	touch $@

CSF_Spectronaut_protein_swap/good_data/quantile/swap/.stamp.bundle: CSF_Spectronaut_protein_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_nonmsstats_block.R CSF_Spectronaut_protein_swap/good_data quantile
	touch $@

CSF_Spectronaut_protein_swap/small/log2/swap/.stamp.bundle: CSF_Spectronaut_protein_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_nonmsstats_block.R CSF_Spectronaut_protein_swap/small log2
	touch $@

CSF_Spectronaut_protein_swap/small/median/swap/.stamp.bundle: CSF_Spectronaut_protein_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_nonmsstats_block.R CSF_Spectronaut_protein_swap/small median
	touch $@

CSF_Spectronaut_protein_swap/small/quantile/swap/.stamp.bundle: CSF_Spectronaut_protein_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_nonmsstats_block.R CSF_Spectronaut_protein_swap/small quantile
	touch $@

CSF_Spectronaut_protein_swap/small_good_data/log2/swap/.stamp.bundle: CSF_Spectronaut_protein_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_nonmsstats_block.R CSF_Spectronaut_protein_swap/small_good_data log2
	touch $@

CSF_Spectronaut_protein_swap/small_good_data/median/swap/.stamp.bundle: CSF_Spectronaut_protein_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_nonmsstats_block.R CSF_Spectronaut_protein_swap/small_good_data median
	touch $@

CSF_Spectronaut_protein_swap/small_good_data/quantile/swap/.stamp.bundle: CSF_Spectronaut_protein_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_nonmsstats_block.R CSF_Spectronaut_protein_swap/small_good_data quantile
	touch $@


# --- Diagnostics (depend on .stamp.bundle; MSstats failures don't block) ----
diag-protein-swap: \
    diag-ps-all_data-log2     diag-ps-all_data-median     diag-ps-all_data-quantile \
    diag-ps-good_data-log2    diag-ps-good_data-median    diag-ps-good_data-quantile \
    diag-ps-small-log2        diag-ps-small-median        diag-ps-small-quantile \
    diag-ps-small_good-log2   diag-ps-small_good-median   diag-ps-small_good-quantile

diag-ps-all_data-log2:        CSF_Spectronaut_protein_swap/all_data/log2/swap/diagnostics.html
diag-ps-all_data-median:      CSF_Spectronaut_protein_swap/all_data/median/swap/diagnostics.html
diag-ps-all_data-quantile:    CSF_Spectronaut_protein_swap/all_data/quantile/swap/diagnostics.html
diag-ps-good_data-log2:       CSF_Spectronaut_protein_swap/good_data/log2/swap/diagnostics.html
diag-ps-good_data-median:     CSF_Spectronaut_protein_swap/good_data/median/swap/diagnostics.html
diag-ps-good_data-quantile:   CSF_Spectronaut_protein_swap/good_data/quantile/swap/diagnostics.html
diag-ps-small-log2:           CSF_Spectronaut_protein_swap/small/log2/swap/diagnostics.html
diag-ps-small-median:         CSF_Spectronaut_protein_swap/small/median/swap/diagnostics.html
diag-ps-small-quantile:       CSF_Spectronaut_protein_swap/small/quantile/swap/diagnostics.html
diag-ps-small_good-log2:      CSF_Spectronaut_protein_swap/small_good_data/log2/swap/diagnostics.html
diag-ps-small_good-median:    CSF_Spectronaut_protein_swap/small_good_data/median/swap/diagnostics.html
diag-ps-small_good-quantile:  CSF_Spectronaut_protein_swap/small_good_data/quantile/swap/diagnostics.html

CSF_Spectronaut_protein_swap/all_data/log2/swap/diagnostics.html: CSF_Spectronaut_protein_swap/all_data/log2/swap/.stamp.bundle vignettes/diagnostics.qmd $(VIGNETTE_HELPERS)
	quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_protein_swap/all_data/log2/swap            -P subset_dir=CSF_Spectronaut_protein_swap/all_data        -P normalization=log2     -P truth_kind=protein_swap -P truth_path=CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv -P base_dir=..

CSF_Spectronaut_protein_swap/all_data/median/swap/diagnostics.html: CSF_Spectronaut_protein_swap/all_data/median/swap/.stamp.bundle vignettes/diagnostics.qmd $(VIGNETTE_HELPERS)
	quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_protein_swap/all_data/median/swap          -P subset_dir=CSF_Spectronaut_protein_swap/all_data        -P normalization=median   -P truth_kind=protein_swap -P truth_path=CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv -P base_dir=..

CSF_Spectronaut_protein_swap/all_data/quantile/swap/diagnostics.html: CSF_Spectronaut_protein_swap/all_data/quantile/swap/.stamp.bundle vignettes/diagnostics.qmd $(VIGNETTE_HELPERS)
	quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_protein_swap/all_data/quantile/swap        -P subset_dir=CSF_Spectronaut_protein_swap/all_data        -P normalization=quantile -P truth_kind=protein_swap -P truth_path=CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv -P base_dir=..

CSF_Spectronaut_protein_swap/good_data/log2/swap/diagnostics.html: CSF_Spectronaut_protein_swap/good_data/log2/swap/.stamp.bundle vignettes/diagnostics.qmd $(VIGNETTE_HELPERS)
	quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_protein_swap/good_data/log2/swap           -P subset_dir=CSF_Spectronaut_protein_swap/good_data       -P normalization=log2     -P truth_kind=protein_swap -P truth_path=CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv -P base_dir=..

CSF_Spectronaut_protein_swap/good_data/median/swap/diagnostics.html: CSF_Spectronaut_protein_swap/good_data/median/swap/.stamp.bundle vignettes/diagnostics.qmd $(VIGNETTE_HELPERS)
	quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_protein_swap/good_data/median/swap         -P subset_dir=CSF_Spectronaut_protein_swap/good_data       -P normalization=median   -P truth_kind=protein_swap -P truth_path=CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv -P base_dir=..

CSF_Spectronaut_protein_swap/good_data/quantile/swap/diagnostics.html: CSF_Spectronaut_protein_swap/good_data/quantile/swap/.stamp.bundle vignettes/diagnostics.qmd $(VIGNETTE_HELPERS)
	quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_protein_swap/good_data/quantile/swap       -P subset_dir=CSF_Spectronaut_protein_swap/good_data       -P normalization=quantile -P truth_kind=protein_swap -P truth_path=CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv -P base_dir=..

CSF_Spectronaut_protein_swap/small/log2/swap/diagnostics.html: CSF_Spectronaut_protein_swap/small/log2/swap/.stamp.bundle vignettes/diagnostics.qmd $(VIGNETTE_HELPERS)
	quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_protein_swap/small/log2/swap               -P subset_dir=CSF_Spectronaut_protein_swap/small           -P normalization=log2     -P truth_kind=protein_swap -P truth_path=CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv -P base_dir=..

CSF_Spectronaut_protein_swap/small/median/swap/diagnostics.html: CSF_Spectronaut_protein_swap/small/median/swap/.stamp.bundle vignettes/diagnostics.qmd $(VIGNETTE_HELPERS)
	quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_protein_swap/small/median/swap             -P subset_dir=CSF_Spectronaut_protein_swap/small           -P normalization=median   -P truth_kind=protein_swap -P truth_path=CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv -P base_dir=..

CSF_Spectronaut_protein_swap/small/quantile/swap/diagnostics.html: CSF_Spectronaut_protein_swap/small/quantile/swap/.stamp.bundle vignettes/diagnostics.qmd $(VIGNETTE_HELPERS)
	quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_protein_swap/small/quantile/swap           -P subset_dir=CSF_Spectronaut_protein_swap/small           -P normalization=quantile -P truth_kind=protein_swap -P truth_path=CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv -P base_dir=..

CSF_Spectronaut_protein_swap/small_good_data/log2/swap/diagnostics.html: CSF_Spectronaut_protein_swap/small_good_data/log2/swap/.stamp.bundle vignettes/diagnostics.qmd $(VIGNETTE_HELPERS)
	quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_protein_swap/small_good_data/log2/swap     -P subset_dir=CSF_Spectronaut_protein_swap/small_good_data -P normalization=log2     -P truth_kind=protein_swap -P truth_path=CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv -P base_dir=..

CSF_Spectronaut_protein_swap/small_good_data/median/swap/diagnostics.html: CSF_Spectronaut_protein_swap/small_good_data/median/swap/.stamp.bundle vignettes/diagnostics.qmd $(VIGNETTE_HELPERS)
	quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_protein_swap/small_good_data/median/swap   -P subset_dir=CSF_Spectronaut_protein_swap/small_good_data -P normalization=median   -P truth_kind=protein_swap -P truth_path=CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv -P base_dir=..

CSF_Spectronaut_protein_swap/small_good_data/quantile/swap/diagnostics.html: CSF_Spectronaut_protein_swap/small_good_data/quantile/swap/.stamp.bundle vignettes/diagnostics.qmd $(VIGNETTE_HELPERS)
	quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_protein_swap/small_good_data/quantile/swap -P subset_dir=CSF_Spectronaut_protein_swap/small_good_data -P normalization=quantile -P truth_kind=protein_swap -P truth_path=CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv -P base_dir=..


# --- Clean -------------------------------------------------------------------
clean-ps-models:
	rm -rf CSF_Spectronaut_protein_swap/all_data/log2/swap
	rm -rf CSF_Spectronaut_protein_swap/all_data/median/swap
	rm -rf CSF_Spectronaut_protein_swap/all_data/quantile/swap
	rm -rf CSF_Spectronaut_protein_swap/good_data/log2/swap
	rm -rf CSF_Spectronaut_protein_swap/good_data/median/swap
	rm -rf CSF_Spectronaut_protein_swap/good_data/quantile/swap
	rm -rf CSF_Spectronaut_protein_swap/small/log2/swap
	rm -rf CSF_Spectronaut_protein_swap/small/median/swap
	rm -rf CSF_Spectronaut_protein_swap/small/quantile/swap
	rm -rf CSF_Spectronaut_protein_swap/small_good_data/log2/swap
	rm -rf CSF_Spectronaut_protein_swap/small_good_data/median/swap
	rm -rf CSF_Spectronaut_protein_swap/small_good_data/quantile/swap

clean-ps-subsets: clean-ps-models
	rm -rf CSF_Spectronaut_protein_swap/all_data
	rm -rf CSF_Spectronaut_protein_swap/good_data
	rm -rf CSF_Spectronaut_protein_swap/small
	rm -rf CSF_Spectronaut_protein_swap/small_good_data
	rm -f  CSF_Spectronaut_protein_swap/.prep.stamp

clean-ps-prep: clean-ps-subsets
	rm -f CSF_Spectronaut_protein_swap/Report.tsv
	rm -f CSF_Spectronaut_protein_swap/annotation.csv
	rm -f CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv
	rm -f CSF_Spectronaut_protein_swap/*_swap_ground_truth.tsv
	rm -f CSF_Spectronaut_protein_swap/*_swap_group_annotation.csv
	rm -f CSF_Spectronaut_protein_swap/*_swap_true_positives.tsv
