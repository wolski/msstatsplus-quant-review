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
# Mix_of_Proteome temporarily excluded — its 3v3 contrast is too small for
# DEqMS's spectraCounteBayes; needs multi-condition annotation + 2-line
# tweaks to R/models_msqrob2.R and R/models_msstats.R. Re-enable by
# un-commenting this line and removing the `mix` references in the
# aggregator targets below.
# include Mix_of_Proteome/Makefile


.PHONY: all symlinks prep cells diagnostics review \
        review-html review-pdf review-best-effort \
        render-vignettes gh-pages \
        clean clean-models clean-subsets clean-prep \
        all_log2 all_median all_quantile \
        all_log2_swap \
        cells-log2 cells-median cells-quantile cells-log2-swap \
        diag-log2  diag-median  diag-quantile diag-log2-swap


all: symlinks prep cells diagnostics review


# =============================================================================
# Per-normalization shortcuts. `make all_<norm>` runs symlinks + prep + only
# the <norm> cells + their diagnostics + a best-effort review render.
#
# CSF_Spectronaut is log2-only (authors' replication), so it appears only in
# all_log2. protein_swap and sample_swap appear in all three norms.
# =============================================================================
cells-log2: \
    cells-csf-all_data-log2 cells-csf-good_data-log2 \
    cells-ps-all_data-log2  cells-ps-good_data-log2  cells-ps-small-log2  cells-ps-small_good-log2 \
    cells-ss-all_data-log2  cells-ss-good_data-log2  cells-ss-small_good-log2

cells-median: \
    cells-ps-all_data-median  cells-ps-good_data-median  cells-ps-small-median  cells-ps-small_good-median \
    cells-ss-all_data-median  cells-ss-good_data-median  cells-ss-small_good-median

cells-quantile: \
    cells-ps-all_data-quantile  cells-ps-good_data-quantile  cells-ps-small-quantile  cells-ps-small_good-quantile \
    cells-ss-all_data-quantile  cells-ss-good_data-quantile  cells-ss-small_good-quantile

diag-log2: \
    diag-csf-all_data-log2 diag-csf-good_data-log2 \
    diag-ps-all_data-log2  diag-ps-good_data-log2  diag-ps-small-log2  diag-ps-small_good-log2 \
    diag-ss-all_data-log2  diag-ss-good_data-log2  diag-ss-small_good-log2

diag-median: \
    diag-ps-all_data-median  diag-ps-good_data-median  diag-ps-small-median  diag-ps-small_good-median \
    diag-ss-all_data-median  diag-ss-good_data-median  diag-ss-small_good-median

diag-quantile: \
    diag-ps-all_data-quantile  diag-ps-good_data-quantile  diag-ps-small-quantile  diag-ps-small_good-quantile \
    diag-ss-all_data-quantile  diag-ss-good_data-quantile  diag-ss-small_good-quantile

all_log2:     symlinks prep cells-log2     diag-log2     review-best-effort
all_median:   symlinks prep cells-median   diag-median   review-best-effort
all_quantile: symlinks prep cells-quantile diag-quantile review-best-effort

# Swap-folders only (drops CSF_Spectronaut authors' replication).
# CSF is log2-only, so all_median / all_quantile already exclude it — only the
# log2 case needs a swap-only variant.
cells-log2-swap: \
    cells-ps-all_data-log2  cells-ps-good_data-log2  cells-ps-small-log2  cells-ps-small_good-log2 \
    cells-ss-all_data-log2  cells-ss-good_data-log2  cells-ss-small_good-log2

diag-log2-swap: \
    diag-ps-all_data-log2   diag-ps-good_data-log2   diag-ps-small-log2   diag-ps-small_good-log2 \
    diag-ss-all_data-log2   diag-ss-good_data-log2   diag-ss-small_good-log2

all_log2_swap: symlinks-csf prep-protein-swap prep-sample-swap cells-log2-swap diag-log2-swap review-best-effort


# =============================================================================
# §0. SYMLINKS — delegated to folder Makefiles (so standalone invocation works).
# =============================================================================
symlinks: symlinks-csf
# symlinks-mix     # Mix_of_Proteome temporarily excluded


# =============================================================================
# §1-3. Aggregators (folder Makefiles supply the per-folder rules).
# =============================================================================
prep:        prep-csf prep-protein-swap prep-sample-swap     # prep-mix excluded
cells:       cells-csf cells-protein-swap cells-sample-swap  # cells-mix excluded
diagnostics: diag-csf  diag-protein-swap  diag-sample-swap   # diag-mix excluded


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
    CSF_Spectronaut/all_data/log2/swap/.stamp.bundle \
    CSF_Spectronaut/good_data/log2/swap/.stamp.bundle \
    CSF_Spectronaut_sample_swap/all_data/log2/swap/.stamp.bundle \
    CSF_Spectronaut_sample_swap/all_data/median/swap/.stamp.bundle \
    CSF_Spectronaut_sample_swap/all_data/quantile/swap/.stamp.bundle \
    CSF_Spectronaut_sample_swap/good_data/log2/swap/.stamp.bundle \
    CSF_Spectronaut_sample_swap/good_data/median/swap/.stamp.bundle \
    CSF_Spectronaut_sample_swap/good_data/quantile/swap/.stamp.bundle \
    CSF_Spectronaut_sample_swap/small_good_data/log2/swap/.stamp.bundle \
    CSF_Spectronaut_sample_swap/small_good_data/median/swap/.stamp.bundle \
    CSF_Spectronaut_sample_swap/small_good_data/quantile/swap/.stamp.bundle \
    CSF_Spectronaut_protein_swap/all_data/log2/swap/.stamp.bundle \
    CSF_Spectronaut_protein_swap/all_data/median/swap/.stamp.bundle \
    CSF_Spectronaut_protein_swap/all_data/quantile/swap/.stamp.bundle \
    CSF_Spectronaut_protein_swap/good_data/log2/swap/.stamp.bundle \
    CSF_Spectronaut_protein_swap/good_data/median/swap/.stamp.bundle \
    CSF_Spectronaut_protein_swap/good_data/quantile/swap/.stamp.bundle \
    CSF_Spectronaut_protein_swap/small/log2/swap/.stamp.bundle \
    CSF_Spectronaut_protein_swap/small/median/swap/.stamp.bundle \
    CSF_Spectronaut_protein_swap/small/quantile/swap/.stamp.bundle \
    CSF_Spectronaut_protein_swap/small_good_data/log2/swap/.stamp.bundle \
    CSF_Spectronaut_protein_swap/small_good_data/median/swap/.stamp.bundle \
    CSF_Spectronaut_protein_swap/small_good_data/quantile/swap/.stamp.bundle
    # Mix_of_Proteome stamps excluded for now

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
                    # Mix_of_Proteome/.prep.stamp     # excluded for now
	cd vignettes && quarto render review.qmd --to html
	cd vignettes && quarto render review.qmd --to pdf


# =============================================================================
# GitHub Pages — render and staging are deliberately separate. Chain:
#   make render-vignettes && make gh-pages
# `gh-pages` does not depend on `render-vignettes`; the publish script
# aborts if any HTML it needs is missing in vignettes/, rather than
# re-rendering silently. See TODO/TODO_ghpages.md.
# =============================================================================
render-vignettes:
	bash src/render_vignettes.sh

gh-pages:
	bash src/publish_gh_pages.sh


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

clean-models:  clean-csf-models  clean-ps-models  clean-ss-models  clean-root-debris   # clean-mix-models excluded
clean-subsets: clean-csf-subsets clean-ps-subsets clean-ss-subsets                      # clean-mix-subsets excluded
clean-prep:    clean-csf-prep    clean-ps-prep    clean-ss-prep                         # clean-mix-prep excluded

# MSstats writes log files to cwd when invoked. Cells for protein_swap,
# sample_swap, and Mix_of_Proteome are run from quant/ root, so their debris
# lands here.
.PHONY: clean-root-debris
clean-root-debris:
	rm -f MSstats_*.log
	rm -f Rplots.pdf
