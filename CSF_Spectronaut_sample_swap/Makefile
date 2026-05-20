# CSF_Spectronaut_sample_swap/Makefile — sample-swap benchmark folder.
# 3 subsets × 3 norms × 7 packages = 63 model fits across 9 cell-blocks.
# Reads the protein universe directly from CSF_Spectronaut/Report.tsv;
# no dependency on the protein-swap folder.
# cwd assumed to be quant/. Invoke from quant/ as:
#   make -f CSF_Spectronaut_sample_swap/Makefile <target>
# or include from the top-level Makefile.

include mk/common.mk

.PHONY: prep-sample-swap cells-sample-swap diag-sample-swap \
        cells-ss-all_data-log2      cells-ss-all_data-median      cells-ss-all_data-quantile \
        cells-ss-good_data-log2     cells-ss-good_data-median     cells-ss-good_data-quantile \
        cells-ss-small_good-log2    cells-ss-small_good-median    cells-ss-small_good-quantile \
        diag-ss-all_data-log2       diag-ss-all_data-median       diag-ss-all_data-quantile \
        diag-ss-good_data-log2      diag-ss-good_data-median      diag-ss-good_data-quantile \
        diag-ss-small_good-log2     diag-ss-small_good-median     diag-ss-small_good-quantile \
        clean-ss-models clean-ss-subsets clean-ss-prep


# --- Prep --------------------------------------------------------------------
prep-sample-swap: CSF_Spectronaut_sample_swap/.prep.stamp

CSF_Spectronaut_sample_swap/.prep.stamp: CSF_Spectronaut/Report.tsv CSF_Spectronaut/annotation.csv R/build_sample_swap_list.R src/swap_spectronaut_report_samples.py src/build_subsets.py
	Rscript R/build_sample_swap_list.R
	$(PYTHON) src/swap_spectronaut_report_samples.py --report CSF_Spectronaut/Report.tsv --annotation CSF_Spectronaut/annotation.csv --protein-swap-list CSF_Spectronaut_sample_swap/CSF_protein_swap_list.csv --out-dir CSF_Spectronaut_sample_swap
	$(PYTHON) src/build_subsets.py --report CSF_Spectronaut_sample_swap/Report.tsv --annotation CSF_Spectronaut_sample_swap/annotation.csv --out-dir CSF_Spectronaut_sample_swap --subsets all_data good_data small_good_data --good-rule label_good
	touch $@


# --- Cells -------------------------------------------------------------------
cells-sample-swap: \
    cells-ss-all_data-log2     cells-ss-all_data-median     cells-ss-all_data-quantile \
    cells-ss-good_data-log2    cells-ss-good_data-median    cells-ss-good_data-quantile \
    cells-ss-small_good-log2   cells-ss-small_good-median   cells-ss-small_good-quantile

cells-ss-all_data-log2:       CSF_Spectronaut_sample_swap/all_data/log2/swap/.stamp
cells-ss-all_data-median:     CSF_Spectronaut_sample_swap/all_data/median/swap/.stamp
cells-ss-all_data-quantile:   CSF_Spectronaut_sample_swap/all_data/quantile/swap/.stamp
cells-ss-good_data-log2:      CSF_Spectronaut_sample_swap/good_data/log2/swap/.stamp
cells-ss-good_data-median:    CSF_Spectronaut_sample_swap/good_data/median/swap/.stamp
cells-ss-good_data-quantile:  CSF_Spectronaut_sample_swap/good_data/quantile/swap/.stamp
cells-ss-small_good-log2:     CSF_Spectronaut_sample_swap/small_good_data/log2/swap/.stamp
cells-ss-small_good-median:   CSF_Spectronaut_sample_swap/small_good_data/median/swap/.stamp
cells-ss-small_good-quantile: CSF_Spectronaut_sample_swap/small_good_data/quantile/swap/.stamp

CSF_Spectronaut_sample_swap/all_data/log2/swap/.stamp: CSF_Spectronaut_sample_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_cell.R CSF_Spectronaut_sample_swap/all_data log2 MSstats
	Rscript R/run_cell.R CSF_Spectronaut_sample_swap/all_data log2 "MSstats+"
	Rscript R/run_nonmsstats_block.R CSF_Spectronaut_sample_swap/all_data log2
	touch $@

CSF_Spectronaut_sample_swap/all_data/median/swap/.stamp: CSF_Spectronaut_sample_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_cell.R CSF_Spectronaut_sample_swap/all_data median MSstats
	Rscript R/run_cell.R CSF_Spectronaut_sample_swap/all_data median "MSstats+"
	Rscript R/run_nonmsstats_block.R CSF_Spectronaut_sample_swap/all_data median
	touch $@

CSF_Spectronaut_sample_swap/all_data/quantile/swap/.stamp: CSF_Spectronaut_sample_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_cell.R CSF_Spectronaut_sample_swap/all_data quantile MSstats
	Rscript R/run_cell.R CSF_Spectronaut_sample_swap/all_data quantile "MSstats+"
	Rscript R/run_nonmsstats_block.R CSF_Spectronaut_sample_swap/all_data quantile
	touch $@

CSF_Spectronaut_sample_swap/good_data/log2/swap/.stamp: CSF_Spectronaut_sample_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_cell.R CSF_Spectronaut_sample_swap/good_data log2 MSstats
	Rscript R/run_cell.R CSF_Spectronaut_sample_swap/good_data log2 "MSstats+"
	Rscript R/run_nonmsstats_block.R CSF_Spectronaut_sample_swap/good_data log2
	touch $@

CSF_Spectronaut_sample_swap/good_data/median/swap/.stamp: CSF_Spectronaut_sample_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_cell.R CSF_Spectronaut_sample_swap/good_data median MSstats
	Rscript R/run_cell.R CSF_Spectronaut_sample_swap/good_data median "MSstats+"
	Rscript R/run_nonmsstats_block.R CSF_Spectronaut_sample_swap/good_data median
	touch $@

CSF_Spectronaut_sample_swap/good_data/quantile/swap/.stamp: CSF_Spectronaut_sample_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_cell.R CSF_Spectronaut_sample_swap/good_data quantile MSstats
	Rscript R/run_cell.R CSF_Spectronaut_sample_swap/good_data quantile "MSstats+"
	Rscript R/run_nonmsstats_block.R CSF_Spectronaut_sample_swap/good_data quantile
	touch $@

CSF_Spectronaut_sample_swap/small_good_data/log2/swap/.stamp: CSF_Spectronaut_sample_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_cell.R CSF_Spectronaut_sample_swap/small_good_data log2 MSstats
	Rscript R/run_cell.R CSF_Spectronaut_sample_swap/small_good_data log2 "MSstats+"
	Rscript R/run_nonmsstats_block.R CSF_Spectronaut_sample_swap/small_good_data log2
	touch $@

CSF_Spectronaut_sample_swap/small_good_data/median/swap/.stamp: CSF_Spectronaut_sample_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_cell.R CSF_Spectronaut_sample_swap/small_good_data median MSstats
	Rscript R/run_cell.R CSF_Spectronaut_sample_swap/small_good_data median "MSstats+"
	Rscript R/run_nonmsstats_block.R CSF_Spectronaut_sample_swap/small_good_data median
	touch $@

CSF_Spectronaut_sample_swap/small_good_data/quantile/swap/.stamp: CSF_Spectronaut_sample_swap/.prep.stamp $(SCRIPTS)
	Rscript R/run_cell.R CSF_Spectronaut_sample_swap/small_good_data quantile MSstats
	Rscript R/run_cell.R CSF_Spectronaut_sample_swap/small_good_data quantile "MSstats+"
	Rscript R/run_nonmsstats_block.R CSF_Spectronaut_sample_swap/small_good_data quantile
	touch $@


# --- Diagnostics -------------------------------------------------------------
diag-sample-swap: \
    diag-ss-all_data-log2     diag-ss-all_data-median     diag-ss-all_data-quantile \
    diag-ss-good_data-log2    diag-ss-good_data-median    diag-ss-good_data-quantile \
    diag-ss-small_good-log2   diag-ss-small_good-median   diag-ss-small_good-quantile

diag-ss-all_data-log2:        CSF_Spectronaut_sample_swap/all_data/log2/swap/diagnostics.html
diag-ss-all_data-median:      CSF_Spectronaut_sample_swap/all_data/median/swap/diagnostics.html
diag-ss-all_data-quantile:    CSF_Spectronaut_sample_swap/all_data/quantile/swap/diagnostics.html
diag-ss-good_data-log2:       CSF_Spectronaut_sample_swap/good_data/log2/swap/diagnostics.html
diag-ss-good_data-median:     CSF_Spectronaut_sample_swap/good_data/median/swap/diagnostics.html
diag-ss-good_data-quantile:   CSF_Spectronaut_sample_swap/good_data/quantile/swap/diagnostics.html
diag-ss-small_good-log2:      CSF_Spectronaut_sample_swap/small_good_data/log2/swap/diagnostics.html
diag-ss-small_good-median:    CSF_Spectronaut_sample_swap/small_good_data/median/swap/diagnostics.html
diag-ss-small_good-quantile:  CSF_Spectronaut_sample_swap/small_good_data/quantile/swap/diagnostics.html

CSF_Spectronaut_sample_swap/all_data/log2/swap/diagnostics.html: CSF_Spectronaut_sample_swap/all_data/log2/swap/.stamp vignettes/diagnostics.qmd $(VIGNETTE_HELPERS)
	quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_sample_swap/all_data/log2/swap            -P subset_dir=CSF_Spectronaut_sample_swap/all_data        -P normalization=log2     -P truth_kind=sample_swap -P truth_path=CSF_Spectronaut_sample_swap/CSF_protein_swap_list.csv -P base_dir=..

CSF_Spectronaut_sample_swap/all_data/median/swap/diagnostics.html: CSF_Spectronaut_sample_swap/all_data/median/swap/.stamp vignettes/diagnostics.qmd $(VIGNETTE_HELPERS)
	quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_sample_swap/all_data/median/swap          -P subset_dir=CSF_Spectronaut_sample_swap/all_data        -P normalization=median   -P truth_kind=sample_swap -P truth_path=CSF_Spectronaut_sample_swap/CSF_protein_swap_list.csv -P base_dir=..

CSF_Spectronaut_sample_swap/all_data/quantile/swap/diagnostics.html: CSF_Spectronaut_sample_swap/all_data/quantile/swap/.stamp vignettes/diagnostics.qmd $(VIGNETTE_HELPERS)
	quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_sample_swap/all_data/quantile/swap        -P subset_dir=CSF_Spectronaut_sample_swap/all_data        -P normalization=quantile -P truth_kind=sample_swap -P truth_path=CSF_Spectronaut_sample_swap/CSF_protein_swap_list.csv -P base_dir=..

CSF_Spectronaut_sample_swap/good_data/log2/swap/diagnostics.html: CSF_Spectronaut_sample_swap/good_data/log2/swap/.stamp vignettes/diagnostics.qmd $(VIGNETTE_HELPERS)
	quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_sample_swap/good_data/log2/swap           -P subset_dir=CSF_Spectronaut_sample_swap/good_data       -P normalization=log2     -P truth_kind=sample_swap -P truth_path=CSF_Spectronaut_sample_swap/CSF_protein_swap_list.csv -P base_dir=..

CSF_Spectronaut_sample_swap/good_data/median/swap/diagnostics.html: CSF_Spectronaut_sample_swap/good_data/median/swap/.stamp vignettes/diagnostics.qmd $(VIGNETTE_HELPERS)
	quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_sample_swap/good_data/median/swap         -P subset_dir=CSF_Spectronaut_sample_swap/good_data       -P normalization=median   -P truth_kind=sample_swap -P truth_path=CSF_Spectronaut_sample_swap/CSF_protein_swap_list.csv -P base_dir=..

CSF_Spectronaut_sample_swap/good_data/quantile/swap/diagnostics.html: CSF_Spectronaut_sample_swap/good_data/quantile/swap/.stamp vignettes/diagnostics.qmd $(VIGNETTE_HELPERS)
	quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_sample_swap/good_data/quantile/swap       -P subset_dir=CSF_Spectronaut_sample_swap/good_data       -P normalization=quantile -P truth_kind=sample_swap -P truth_path=CSF_Spectronaut_sample_swap/CSF_protein_swap_list.csv -P base_dir=..

CSF_Spectronaut_sample_swap/small_good_data/log2/swap/diagnostics.html: CSF_Spectronaut_sample_swap/small_good_data/log2/swap/.stamp vignettes/diagnostics.qmd $(VIGNETTE_HELPERS)
	quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_sample_swap/small_good_data/log2/swap     -P subset_dir=CSF_Spectronaut_sample_swap/small_good_data -P normalization=log2     -P truth_kind=sample_swap -P truth_path=CSF_Spectronaut_sample_swap/CSF_protein_swap_list.csv -P base_dir=..

CSF_Spectronaut_sample_swap/small_good_data/median/swap/diagnostics.html: CSF_Spectronaut_sample_swap/small_good_data/median/swap/.stamp vignettes/diagnostics.qmd $(VIGNETTE_HELPERS)
	quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_sample_swap/small_good_data/median/swap   -P subset_dir=CSF_Spectronaut_sample_swap/small_good_data -P normalization=median   -P truth_kind=sample_swap -P truth_path=CSF_Spectronaut_sample_swap/CSF_protein_swap_list.csv -P base_dir=..

CSF_Spectronaut_sample_swap/small_good_data/quantile/swap/diagnostics.html: CSF_Spectronaut_sample_swap/small_good_data/quantile/swap/.stamp vignettes/diagnostics.qmd $(VIGNETTE_HELPERS)
	quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut_sample_swap/small_good_data/quantile/swap -P subset_dir=CSF_Spectronaut_sample_swap/small_good_data -P normalization=quantile -P truth_kind=sample_swap -P truth_path=CSF_Spectronaut_sample_swap/CSF_protein_swap_list.csv -P base_dir=..


# --- Clean -------------------------------------------------------------------
clean-ss-models:
	rm -rf CSF_Spectronaut_sample_swap/all_data/log2/swap
	rm -rf CSF_Spectronaut_sample_swap/all_data/median/swap
	rm -rf CSF_Spectronaut_sample_swap/all_data/quantile/swap
	rm -rf CSF_Spectronaut_sample_swap/good_data/log2/swap
	rm -rf CSF_Spectronaut_sample_swap/good_data/median/swap
	rm -rf CSF_Spectronaut_sample_swap/good_data/quantile/swap
	rm -rf CSF_Spectronaut_sample_swap/small_good_data/log2/swap
	rm -rf CSF_Spectronaut_sample_swap/small_good_data/median/swap
	rm -rf CSF_Spectronaut_sample_swap/small_good_data/quantile/swap

clean-ss-subsets: clean-ss-models
	rm -rf CSF_Spectronaut_sample_swap/all_data
	rm -rf CSF_Spectronaut_sample_swap/good_data
	rm -rf CSF_Spectronaut_sample_swap/small_good_data
	rm -f  CSF_Spectronaut_sample_swap/.prep.stamp

clean-ss-prep: clean-ss-subsets
	rm -f CSF_Spectronaut_sample_swap/Report.tsv
	rm -f CSF_Spectronaut_sample_swap/annotation.csv
	rm -f CSF_Spectronaut_sample_swap/CSF_protein_swap_list.csv
	rm -f CSF_Spectronaut_sample_swap/*_sample_swap_ground_truth.tsv
	rm -f CSF_Spectronaut_sample_swap/*_sample_swap_group_annotation.csv
	rm -f CSF_Spectronaut_sample_swap/*_sample_swap_true_positives.tsv
