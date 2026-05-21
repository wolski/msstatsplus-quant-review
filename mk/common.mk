# mk/common.mk — shared variables for quant/ Makefile and per-folder Makefiles.
# Included by quant/Makefile and each folder's Makefile (guarded so multi-include
# is harmless). cwd is always quant/ when these variables are dereferenced.

ifndef COMMON_MK_INCLUDED
COMMON_MK_INCLUDED := 1

# Python interpreter. Prefer the project's venv (has polars + deps installed).
# Override on the command line with `make PYTHON=/some/other/python ...`.
PYTHON ?= $(if $(wildcard .venv/bin/python3),.venv/bin/python3,python3)

# Shared R model code. Editing any of these invalidates every cell stamp
# across all folders (touch propagates via Make's mtime check).
SCRIPTS = R/run_cell.R R/run_nonmsstats_block.R \
          R/models_msstats.R R/models_maxlfq_limma.R R/models_deqms.R \
          R/models_prolfqua.R R/models_limpa.R R/models_msqrob2.R \
          R/preprocess.R R/paths.R

# Authors' scripts — referenced only by CSF_Spectronaut/Makefile cell rules.
AUTHORS_SCRIPTS = CSF_Spectronaut/run_msstats.R \
                  CSF_Spectronaut/run_nonmsstats.R \
                  CSF_Spectronaut/run_msqrob2_step.R \
                  CSF_Spectronaut/run_step_common.R \
                  CSF_Spectronaut/run_prolfqua_step.R

# Helpers sourced by vignettes/*.qmd at render time. Editing any of these
# re-renders diagnostics + review, but does not invalidate cell stamps.
VIGNETTE_HELPERS = R/figures.R R/comparison_table.R R/ground_truth.R R/paths.R

# Shared raw-input symlink rules. These live here (not in CSF_Spectronaut/Makefile)
# because both swap folders' prep also consumes CSF_Spectronaut/Report.tsv +
# annotation.csv. Putting the rules in common.mk ensures they are visible
# whether you invoke the top-level Makefile or a folder Makefile standalone.
# Mix_of_Proteome's symlinks stay folder-local — only that folder uses them.

.PHONY: symlinks-csf
symlinks-csf: CSF_Spectronaut/Report.tsv CSF_Spectronaut/annotation.csv

CSF_Spectronaut/Report.tsv:
	ln -s "20250130_163144_CSF dilutions Jan 2025 no normalization_Report.tsv" $@

CSF_Spectronaut/annotation.csv:
	ln -s CSF_annotation.csv $@

endif
