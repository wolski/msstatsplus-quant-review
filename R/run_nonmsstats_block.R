## run_nonmsstats_block.R - execute the fast non-MSstats packages for one
## (subset_dir, normalization) cell-block in a single R session, sharing
## the Report.tsv + annotation fread/merge cost across packages.
##
## Usage:
##   Rscript R/run_nonmsstats_block.R <subset_dir> <normalization> [pkg1,pkg2,...]
##
## Defaults to the non-MSstats packages: MaxLFQ_limma, DEqMS, prolfqua, limpa,
## msqrob2. Override by passing an explicit comma-separated package list as
## the optional third argument. Each package writes to
## <subset_dir>/<normalization>/swap/<pkg>/.
## Per-package failures are caught (tryCatch); other packages still run.
## The script exits non-zero if any requested package failed.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 3) {
  stop("Usage: Rscript R/run_nonmsstats_block.R <subset_dir> <normalization> [pkg1,pkg2,...]")
}
subset_dir    <- args[[1]]
normalization <- args[[2]]
pkgs <- if (length(args) == 3) {
  strsplit(args[[3]], ",", fixed = TRUE)[[1]]
} else {
  c("MaxLFQ_limma", "DEqMS", "prolfqua", "limpa", "msqrob2")
}

suppressPackageStartupMessages(library(data.table))
source("R/paths.R")
source("R/preprocess.R")
source("R/models_maxlfq_limma.R")
source("R/models_deqms.R")
source("R/models_prolfqua.R")
source("R/models_limpa.R")
if ("msqrob2" %in% pkgs) {
  source("R/models_msqrob2.R")
}

report_path     <- file.path(subset_dir, "Report.tsv")
annotation_path <- file.path(subset_dir, "annotation.csv")
if (!file.exists(report_path))     stop("Not found: ", report_path)
if (!file.exists(annotation_path)) stop("Not found: ", annotation_path)

t_load <- proc.time()[3]
raw_input  <- fread(report_path, sep = "\t")
annotation <- fread(annotation_path)
raw_input  <- raw_input[tolower(R.Condition) != "blank"]
annotation <- annotation[tolower(Condition) != "blank"]
raw_input  <- raw_input[R.FileName %in% annotation$R.FileName]
merged_input <- merge(raw_input, annotation, by = "R.FileName", all.x = FALSE)
merged_input[, Order := as.integer(stringr::str_split_i(R.FileName, "Seq", 2))]
cat(sprintf("[block] load+merge: %.1f s; %d rows; %d runs\n",
            proc.time()[3] - t_load, nrow(merged_input), nrow(annotation)))

runners <- list(
  MaxLFQ_limma = function(c) run_maxlfq_limma(merged_input, annotation, normalization, c),
  DEqMS        = function(c) run_deqms(merged_input,        annotation, normalization, c),
  msqrob2      = function(c) run_msqrob2(merged_input,      annotation, normalization, c),
  prolfqua     = function(c) run_prolfqua(merged_input,     annotation, normalization, c),
  limpa        = function(c) run_limpa(merged_input,        annotation, normalization, c)
)

n_fail <- 0L
for (pkg in pkgs) {
  if (!pkg %in% names(runners)) {
    cat(sprintf("[skip] unknown package: %s\n", pkg)); next
  }
  cell <- out_dir(subset_dir, normalization, pkg)
  # Always re-run when invoked. Make decides skip-vs-rebuild at the stamp
  # level via mtime against $(SCRIPTS); if we're here, the stamp is stale
  # so the model must be regenerated regardless of any existing CSV.
  t0 <- proc.time()[3]
  cat(sprintf("[run ] %s\n", pkg))
  res <- tryCatch(runners[[pkg]](cell),
                  error = function(e) {
                    cat(sprintf("[fail] %s: %s\n", pkg, conditionMessage(e)))
                    NULL
                  })
  # Treat both an outright error (NULL) AND a "soft failure" (the runner
  # returned a result but model is NULL — limpa's tryCatch path) as failures.
  failed <- is.null(res) ||
            (is.list(res) && "model" %in% names(res) && is.null(res$model))
  if (failed) {
    cat(sprintf("[fail] %s: no model produced\n", pkg))
    n_fail <- n_fail + 1L
  } else {
    cat(sprintf("[ ok ] %s in %.1fs\n", pkg, proc.time()[3] - t0))
  }
}

if (n_fail > 0L) {
  stop(sprintf("[block] %d package(s) failed", n_fail))
}
