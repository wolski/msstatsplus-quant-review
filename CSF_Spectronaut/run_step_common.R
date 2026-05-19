## Common setup for the CSF_Spectronaut authors' replication runners.
## Reads pre-prepared subset data from CSF_Spectronaut/<subset>/Report.tsv +
## annotation.csv (written by src/build_subsets.py). Subset name comes from
## the first positional commandArg.
##
## Provides (for downstream sourcers):
##   subset_name, variant, normalization,
##   raw_input, annotation, merged_input, run_order,
##   protein_swap_list, true_positives, all_proteins, no_swap,
##   out_dir(method), label_proteins(df), tic(), toc(t0), write_timing(...),
##   vsn_normalize_matrix(), quantile_normalize_log2_matrix(), median_normalize_log2_matrix().

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
})
source("../benchmark_experiments_functions.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript <runner>.R <subset>   (subset in: all_data, good_data)")
}
subset_name <- args[[1]]
stopifnot(subset_name %in% c("all_data", "good_data"))

# Fixed pipeline parameters for this authors' replication folder.
variant       <- "V1_log2"
normalization <- "none"

subset_dir  <- subset_name  # cwd is CSF_Spectronaut/
report_path <- file.path(subset_dir, "Report.tsv")
ann_path    <- file.path(subset_dir, "annotation.csv")
if (!file.exists(report_path)) stop("Not found: ", report_path,
                                    " (run `make prep-csf` first)")
if (!file.exists(ann_path))    stop("Not found: ", ann_path,
                                    " (run `make prep-csf` first)")

cat(sprintf("[step] subset=%s report=%s\n", subset_name, report_path))

raw_input  <- fread(report_path, sep = "\t")
annotation <- fread(ann_path)
raw_input  <- raw_input[tolower(raw_input$R.Condition) != "blank", ]
annotation <- annotation[tolower(annotation$Condition) != "blank", ]

annotation$Run <- annotation$R.FileName
run_order <- unique(annotation[, .(Run, Order)])

merged_input <- merge(raw_input, annotation, by = "R.FileName",
                      all.x = TRUE, all.y = FALSE)

protein_swap_list <- fread("CSF_protein_swap_list.csv")
true_positives <- protein_swap_list[Label == "Positive", Protein]
all_proteins   <- protein_swap_list$Protein
no_swap        <- character(0)

# Canonical output path: <subset>/log2/swap/<method>/  (relative to CSF_Spectronaut/)
out_dir <- function(method) {
  d <- file.path(subset_dir, "log2", "swap", method)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

label_proteins <- function(df, protein_col = "Protein") {
  df$Label <- ifelse(df[[protein_col]] %in% true_positives,
                     "Positive", "Negative")
  df
}

# Per-method wall-clock timing (preprocess + model phases). Each runner
# calls tic()/toc() and then write_timing() to emit <method>_timing.csv.
tic <- function() proc.time()[3]
toc <- function(t0) as.numeric(proc.time()[3] - t0)

write_timing <- function(method, method_dir,
                         preprocess_seconds, model_seconds) {
  fwrite(
    data.frame(method = method,
               preprocess_seconds = preprocess_seconds,
               model_seconds      = model_seconds),
    file = file.path(method_dir, paste0(method, "_timing.csv"))
  )
  cat(sprintf("[time] %-18s subset=%s pre=%6.1f s  fit=%6.1f s\n",
              method, subset_name, preprocess_seconds, model_seconds))
  invisible(NULL)
}

# Normalization helpers retained for compatibility with run_nonmsstats.R
# (which references them via the `normalization` constant; they are not
# applied at variant=V1_log2 / normalization=none).
vsn_normalize_matrix <- function(m) {
  m <- as.matrix(m); storage.mode(m) <- "double"; dn <- dimnames(m)
  m[!is.finite(m)] <- NA_real_
  out <- vsn::justvsn(m); dimnames(out) <- dn; out
}
quantile_normalize_log2_matrix <- function(m) {
  m <- as.matrix(m); storage.mode(m) <- "double"; dn <- dimnames(m)
  m[!is.finite(m)] <- NA_real_
  out <- limma::normalizeBetweenArrays(m, method = "quantile")
  dimnames(out) <- dn; out
}
median_normalize_log2_matrix <- function(m) {
  m <- as.matrix(m); storage.mode(m) <- "double"; dn <- dimnames(m)
  m[!is.finite(m)] <- NA_real_
  cm <- matrixStats::colMedians(m, na.rm = TRUE)
  out <- sweep(m, 2, cm, FUN = "-"); dimnames(out) <- dn; out
}
apply_vsn      <- FALSE
apply_quantile <- FALSE
apply_median   <- FALSE
