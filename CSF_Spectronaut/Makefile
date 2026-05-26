# CSF_Spectronaut/Makefile — authors' replication folder. log2 only.
# cwd assumed to be quant/. Invoke from quant/ as:
#   make -f CSF_Spectronaut/Makefile <target>
# or include from the top-level Makefile.
#
# Each cell-block has TWO stamps:
#   .stamp.msstats — MSstats + MSstats+ both done
#   .stamp.bundle  — non-MSstats bundle done
# Downstream (diagnostics, review) depends on .stamp.bundle so a failure
# in MSstats+ does not block the rest of the pipeline.

include mk/common.mk

.PHONY: prep-csf cells-csf diag-csf \
        cells-csf-all_data-log2 cells-csf-good_data-log2 \
        diag-csf-all_data-log2  diag-csf-good_data-log2 \
        clean-csf-models clean-csf-subsets clean-csf-prep


# --- Symlinks ----------------------------------------------------------------
# Raw-input symlink rules live in mk/common.mk so the protein-swap and
# sample-swap Makefiles can use them standalone.


# --- Prep --------------------------------------------------------------------
prep-csf: CSF_Spectronaut/.prep.stamp

CSF_Spectronaut/.prep.stamp: CSF_Spectronaut/Report.tsv CSF_Spectronaut/annotation.csv src/build_subsets.py
	$(PYTHON) src/build_subsets.py --report CSF_Spectronaut/Report.tsv --annotation CSF_Spectronaut/annotation.csv --out-dir CSF_Spectronaut --subsets all_data --subsets good_data --good-rule label_good
	touch $@


# --- Cells (authors' scripts; swap done post-summarisation inside the R code) -
cells-csf: cells-csf-all_data-log2 cells-csf-good_data-log2

cells-csf-all_data-log2: \
    CSF_Spectronaut/all_data/log2/swap/.stamp.msstats \
    CSF_Spectronaut/all_data/log2/swap/.stamp.bundle

cells-csf-good_data-log2: \
    CSF_Spectronaut/good_data/log2/swap/.stamp.msstats \
    CSF_Spectronaut/good_data/log2/swap/.stamp.bundle

CSF_Spectronaut/all_data/log2/swap/.stamp.msstats: CSF_Spectronaut/.prep.stamp $(AUTHORS_SCRIPTS)
	cd CSF_Spectronaut && Rscript run_msstats.R all_data
	touch $@

CSF_Spectronaut/all_data/log2/swap/.stamp.bundle: CSF_Spectronaut/.prep.stamp $(AUTHORS_SCRIPTS)
	cd CSF_Spectronaut && RUN_MSQROB2=1 Rscript run_nonmsstats.R all_data
	touch $@

CSF_Spectronaut/good_data/log2/swap/.stamp.msstats: CSF_Spectronaut/.prep.stamp $(AUTHORS_SCRIPTS)
	cd CSF_Spectronaut && Rscript run_msstats.R good_data
	touch $@

CSF_Spectronaut/good_data/log2/swap/.stamp.bundle: CSF_Spectronaut/.prep.stamp $(AUTHORS_SCRIPTS)
	cd CSF_Spectronaut && RUN_MSQROB2=1 Rscript run_nonmsstats.R good_data
	touch $@


# --- Diagnostics -------------------------------------------------------------
diag-csf: diag-csf-all_data-log2 diag-csf-good_data-log2

diag-csf-all_data-log2:  CSF_Spectronaut/all_data/log2/swap/diagnostics.html
diag-csf-good_data-log2: CSF_Spectronaut/good_data/log2/swap/diagnostics.html

CSF_Spectronaut/all_data/log2/swap/diagnostics.html: CSF_Spectronaut/all_data/log2/swap/.stamp.bundle vignettes/diagnostics.qmd $(VIGNETTE_HELPERS)
	quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut/all_data/log2/swap  -P subset_dir=CSF_Spectronaut/all_data  -P normalization=log2 -P truth_kind=sample_swap -P truth_path=CSF_Spectronaut/CSF_protein_swap_list.csv -P base_dir=..

CSF_Spectronaut/good_data/log2/swap/diagnostics.html: CSF_Spectronaut/good_data/log2/swap/.stamp.bundle vignettes/diagnostics.qmd $(VIGNETTE_HELPERS)
	quarto render vignettes/diagnostics.qmd --output-dir ../CSF_Spectronaut/good_data/log2/swap -P subset_dir=CSF_Spectronaut/good_data -P normalization=log2 -P truth_kind=sample_swap -P truth_path=CSF_Spectronaut/CSF_protein_swap_list.csv -P base_dir=..


# --- Clean -------------------------------------------------------------------
# rm -rf swap/ removes both .stamp.msstats and .stamp.bundle along with model dirs.
clean-csf-models:
	rm -rf CSF_Spectronaut/all_data/log2/swap
	rm -rf CSF_Spectronaut/good_data/log2/swap
	rm -f  CSF_Spectronaut/MSstats_*.log
	rm -f  CSF_Spectronaut/Rplots.pdf

clean-csf-subsets: clean-csf-models
	rm -rf CSF_Spectronaut/all_data
	rm -rf CSF_Spectronaut/good_data
	rm -f  CSF_Spectronaut/.prep.stamp

clean-csf-prep: clean-csf-subsets
	@true
