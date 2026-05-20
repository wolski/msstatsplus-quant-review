# quant/Makefile — top-level orchestrator.
#
# Includes the four per-folder Makefiles (each editable in isolation) plus
# shared variables from mk/common.mk. Single dependency graph; `make -jN`
# parallelism spans folders without any recursive-make magic.
#
# Invoke from quant/:
#   make all                           # symlinks + prep + cells + diagnostics + review
#   make symlinks
#   make prep                          # all four folders' prep
#   make cells                         # all 26 cell-blocks
#   make diagnostics                   # all 26 diagnostics HTMLs
#   make review                        # review.html + review.pdf
#   make review-best-effort            # render review against partial outputs
#   make cells-csf                     # one folder's cells (and similarly for ps/ss/mix)
#   make clean-models | clean-subsets | clean-prep
#
# Parallelism: each cell-block runs MSstats + MSstats+ + bundle sequentially
# (MSstats internally uses 8 cores). Distinct cell-blocks run in parallel.
#   -j 2  : two cell-blocks at a time (~16 cores). Laptop default.
#   -j 4  : ~32 cores. Workstation.

include mk/common.mk

include CSF_Spectronaut/Makefile
include CSF_Spectronaut_protein_swap/Makefile
include CSF_Spectronaut_sample_swap/Makefile
include Mix_of_Proteome/Makefile


.PHONY: all symlinks prep cells diagnostics review \
        review-html review-pdf review-best-effort \
        clean clean-models clean-subsets clean-prep


all: symlinks prep cells diagnostics review


# =============================================================================
# §0. SYMLINKS — delegated to folder Makefiles (so standalone invocation works).
# =============================================================================
symlinks: symlinks-csf symlinks-mix


# =============================================================================
# §1-3. Aggregators (folder Makefiles supply the per-folder rules).
# =============================================================================
prep:        prep-csf prep-protein-swap prep-sample-swap prep-mix
cells:       cells-csf cells-protein-swap cells-sample-swap cells-mix
diagnostics: diag-csf  diag-protein-swap  diag-sample-swap  diag-mix


# =============================================================================
# §4. REVIEW
# =============================================================================
review: review-html review-pdf

# Strict review: depends on every cell stamp + the qmd + helpers.
# KNOWN ISSUE: some cells fail upstream (MSstats+ × quantile across multiple
# subsets, DEqMS × small × log2). If a stamp is permanently missing,
# `make review` will block. Workarounds:
#   (a) `make -k review` — keep-going; review.qmd guards each chunk with
#       file.exists() and skips missing models.
#   (b) `make review-best-effort` — render against whatever is on disk
#       without requiring stamps.
REVIEW_DEPS = vignettes/review.qmd $(VIGNETTE_HELPERS) \
    CSF_Spectronaut/all_data/log2/swap/.stamp                            CSF_Spectronaut/good_data/log2/swap/.stamp \
    CSF_Spectronaut_sample_swap/all_data/log2/swap/.stamp                CSF_Spectronaut_sample_swap/all_data/median/swap/.stamp                CSF_Spectronaut_sample_swap/all_data/quantile/swap/.stamp \
    CSF_Spectronaut_sample_swap/good_data/log2/swap/.stamp               CSF_Spectronaut_sample_swap/good_data/median/swap/.stamp               CSF_Spectronaut_sample_swap/good_data/quantile/swap/.stamp \
    CSF_Spectronaut_sample_swap/small_good_data/log2/swap/.stamp         CSF_Spectronaut_sample_swap/small_good_data/median/swap/.stamp         CSF_Spectronaut_sample_swap/small_good_data/quantile/swap/.stamp \
    CSF_Spectronaut_protein_swap/all_data/log2/swap/.stamp               CSF_Spectronaut_protein_swap/all_data/median/swap/.stamp               CSF_Spectronaut_protein_swap/all_data/quantile/swap/.stamp \
    CSF_Spectronaut_protein_swap/good_data/log2/swap/.stamp              CSF_Spectronaut_protein_swap/good_data/median/swap/.stamp              CSF_Spectronaut_protein_swap/good_data/quantile/swap/.stamp \
    CSF_Spectronaut_protein_swap/small/log2/swap/.stamp                  CSF_Spectronaut_protein_swap/small/median/swap/.stamp                  CSF_Spectronaut_protein_swap/small/quantile/swap/.stamp \
    CSF_Spectronaut_protein_swap/small_good_data/log2/swap/.stamp        CSF_Spectronaut_protein_swap/small_good_data/median/swap/.stamp        CSF_Spectronaut_protein_swap/small_good_data/quantile/swap/.stamp \
    Mix_of_Proteome/all_data/log2/swap/.stamp                            Mix_of_Proteome/all_data/median/swap/.stamp                            Mix_of_Proteome/all_data/quantile/swap/.stamp

review-html: $(REVIEW_DEPS)
	cd vignettes && quarto render review.qmd --to html

review-pdf: $(REVIEW_DEPS)
	cd vignettes && quarto render review.qmd --to pdf

# review-best-effort: render review.qmd against whatever cell outputs happen
# to exist on disk. Bypasses the cell stamps (so missing/broken cells don't
# block) but still requires the four .prep.stamps because review.qmd
# unconditionally reads truth files (CSF_protein_swap_list.csv etc.) that
# prep produces. review.qmd's per-cell chunks already guard with
# file.exists() and skip models that aren't on disk.
review-best-effort: vignettes/review.qmd $(VIGNETTE_HELPERS) \
                    CSF_Spectronaut/.prep.stamp \
                    CSF_Spectronaut_protein_swap/.prep.stamp \
                    CSF_Spectronaut_sample_swap/.prep.stamp \
                    Mix_of_Proteome/.prep.stamp
	cd vignettes && quarto render review.qmd --to html
	cd vignettes && quarto render review.qmd --to pdf


# =============================================================================
# §5. CLEAN — layered. `make clean` is a NO-OP (deliberate). Pick the level:
#
#   clean-models   - remove every <subset>/<norm>/swap/ tree (model CSVs +
#                    diagnostics HTML + .stamp).
#   clean-subsets  - also remove every <subset>/ dir (subset Report.tsv +
#                    annotation.csv) and the .prep.stamp files.
#   clean-prep     - also remove swap-script outputs at folder root for the
#                    two swap folders. Keeps the raw Spectronaut TSV symlinks
#                    in CSF_Spectronaut/ and Mix_of_Proteome/.
#
# Per-folder clean targets (clean-csf-*, clean-ps-*, clean-ss-*, clean-mix-*)
# are defined in each folder Makefile.
# =============================================================================
clean:
	@echo "No-op. Use:"
	@echo "  make clean-models   - remove model CSVs + diagnostics HTMLs + stamps"
	@echo "  make clean-subsets  - also remove subset Report.tsv + annotation.csv + prep stamps"
	@echo "  make clean-prep     - also remove swap-script outputs at folder root"
	@echo "Original Spectronaut TSVs and annotations are NEVER touched."

clean-models:  clean-csf-models  clean-ps-models  clean-ss-models  clean-mix-models clean-root-debris
clean-subsets: clean-csf-subsets clean-ps-subsets clean-ss-subsets clean-mix-subsets
clean-prep:    clean-csf-prep    clean-ps-prep    clean-ss-prep    clean-mix-prep

# MSstats writes log files to cwd when invoked. Cells for protein_swap,
# sample_swap, and Mix_of_Proteome are run from quant/ root, so their debris
# lands here.
.PHONY: clean-root-debris
clean-root-debris:
	rm -f MSstats_*.log
	rm -f Rplots.pdf
