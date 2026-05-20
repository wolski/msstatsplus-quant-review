# Mix_of_Proteome/Makefile — species-mix benchmark folder.
# 1 subset × 3 norms × 7 packages = 21 model fits across 3 cell-blocks.
# cwd assumed to be quant/. Invoke from quant/ as:
#   make -f Mix_of_Proteome/Makefile <target>
# or include from the top-level Makefile.

include mk/common.mk

.PHONY: prep-mix cells-mix diag-mix symlinks-mix \
        cells-mix-all_data-log2 cells-mix-all_data-median cells-mix-all_data-quantile \
        diag-mix-all_data-log2  diag-mix-all_data-median  diag-mix-all_data-quantile \
        clean-mix-models clean-mix-subsets clean-mix-prep


# --- Symlinks ----------------------------------------------------------------
# Raw input symlinks live in this folder; folder Makefile owns the rules so
# standalone `make -f Mix_of_Proteome/Makefile prep-mix` works.
symlinks-mix: Mix_of_Proteome/Report.tsv Mix_of_Proteome/annotation.csv

Mix_of_Proteome/Report.tsv:
	ln -s "20250422_140629_OKE_April22_2025_Olsen_astral_benchmark_mix_200ng_NE_report_nodecoy.tsv" $@

Mix_of_Proteome/annotation.csv:
	ln -s Mix_of_Proteome_annotation_contrast.csv $@


# --- Prep --------------------------------------------------------------------
prep-mix: Mix_of_Proteome/.prep.stamp

Mix_of_Proteome/.prep.stamp: Mix_of_Proteome/Report.tsv Mix_of_Proteome/annotation.csv src/build_subsets.py
	$(PYTHON) src/build_subsets.py --report Mix_of_Proteome/Report.tsv --annotation Mix_of_Proteome/annotation.csv --out-dir Mix_of_Proteome --subsets all_data
	touch $@


# --- Cells -------------------------------------------------------------------
cells-mix: cells-mix-all_data-log2 cells-mix-all_data-median cells-mix-all_data-quantile

cells-mix-all_data-log2:     Mix_of_Proteome/all_data/log2/swap/.stamp
cells-mix-all_data-median:   Mix_of_Proteome/all_data/median/swap/.stamp
cells-mix-all_data-quantile: Mix_of_Proteome/all_data/quantile/swap/.stamp

Mix_of_Proteome/all_data/log2/swap/.stamp: Mix_of_Proteome/.prep.stamp $(SCRIPTS)
	Rscript R/run_cell.R Mix_of_Proteome/all_data log2 MSstats
	Rscript R/run_cell.R Mix_of_Proteome/all_data log2 "MSstats+"
	Rscript R/run_nonmsstats_block.R Mix_of_Proteome/all_data log2
	touch $@

Mix_of_Proteome/all_data/median/swap/.stamp: Mix_of_Proteome/.prep.stamp $(SCRIPTS)
	Rscript R/run_cell.R Mix_of_Proteome/all_data median MSstats
	Rscript R/run_cell.R Mix_of_Proteome/all_data median "MSstats+"
	Rscript R/run_nonmsstats_block.R Mix_of_Proteome/all_data median
	touch $@

Mix_of_Proteome/all_data/quantile/swap/.stamp: Mix_of_Proteome/.prep.stamp $(SCRIPTS)
	Rscript R/run_cell.R Mix_of_Proteome/all_data quantile MSstats
	Rscript R/run_cell.R Mix_of_Proteome/all_data quantile "MSstats+"
	Rscript R/run_nonmsstats_block.R Mix_of_Proteome/all_data quantile
	touch $@


# --- Diagnostics -------------------------------------------------------------
diag-mix: diag-mix-all_data-log2 diag-mix-all_data-median diag-mix-all_data-quantile

diag-mix-all_data-log2:     Mix_of_Proteome/all_data/log2/swap/diagnostics.html
diag-mix-all_data-median:   Mix_of_Proteome/all_data/median/swap/diagnostics.html
diag-mix-all_data-quantile: Mix_of_Proteome/all_data/quantile/swap/diagnostics.html

Mix_of_Proteome/all_data/log2/swap/diagnostics.html: Mix_of_Proteome/all_data/log2/swap/.stamp vignettes/diagnostics.qmd $(VIGNETTE_HELPERS)
	quarto render vignettes/diagnostics.qmd --output-dir ../Mix_of_Proteome/all_data/log2/swap     -P subset_dir=Mix_of_Proteome/all_data -P normalization=log2     -P truth_kind=mix_of_proteome -P truth_path=Mix_of_Proteome/idmapping.tsv -P base_dir=..

Mix_of_Proteome/all_data/median/swap/diagnostics.html: Mix_of_Proteome/all_data/median/swap/.stamp vignettes/diagnostics.qmd $(VIGNETTE_HELPERS)
	quarto render vignettes/diagnostics.qmd --output-dir ../Mix_of_Proteome/all_data/median/swap   -P subset_dir=Mix_of_Proteome/all_data -P normalization=median   -P truth_kind=mix_of_proteome -P truth_path=Mix_of_Proteome/idmapping.tsv -P base_dir=..

Mix_of_Proteome/all_data/quantile/swap/diagnostics.html: Mix_of_Proteome/all_data/quantile/swap/.stamp vignettes/diagnostics.qmd $(VIGNETTE_HELPERS)
	quarto render vignettes/diagnostics.qmd --output-dir ../Mix_of_Proteome/all_data/quantile/swap -P subset_dir=Mix_of_Proteome/all_data -P normalization=quantile -P truth_kind=mix_of_proteome -P truth_path=Mix_of_Proteome/idmapping.tsv -P base_dir=..


# --- Clean -------------------------------------------------------------------
clean-mix-models:
	rm -rf Mix_of_Proteome/all_data/log2/swap
	rm -rf Mix_of_Proteome/all_data/median/swap
	rm -rf Mix_of_Proteome/all_data/quantile/swap

clean-mix-subsets: clean-mix-models
	rm -rf Mix_of_Proteome/all_data
	rm -f  Mix_of_Proteome/.prep.stamp

clean-mix-prep: clean-mix-subsets
	@true
