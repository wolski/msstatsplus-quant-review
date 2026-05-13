## Common setup for standalone single-method runners.
## Provides: report_path, variant, out_suffix, apply_vsn, vsn_normalize_matrix,
##            merged_input, annotation, run_order, true_positives,
##            all_proteins, no_swap, out_dir(), label_proteins().

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
})
source("../benchmark_experiments_functions.R")

report_path = Sys.getenv(
  "REPORT_PATH",
  unset = "20250130_163144_CSF dilutions Jan 2025 no normalization_Report.tsv"
)
variant    = Sys.getenv("VARIANT",    unset = "V1_log2")
out_suffix = Sys.getenv("OUT_SUFFIX", unset = "")
# OUT_TAG is now the parent run directory (e.g. "all_dilutions" or
# "no_high_dilutions"), so outputs land in <OUT_TAG>/<variant><OUT_SUFFIX>/.
# Defaults to "all_dilutions" when not specified.
out_tag    = Sys.getenv("OUT_TAG",    unset = "all_dilutions")
exclude_dilutions = Sys.getenv("EXCLUDE_DILUTIONS", unset = "")
exclude_dilutions = if (nchar(exclude_dilutions) > 0) {
  trimws(strsplit(exclude_dilutions, ",", fixed = TRUE)[[1]])
} else {
  character(0)
}
# NORMALIZATION:
#   "none" (default)        : log2 only, no inter-sample normalization
#   "equalizeMedians"       : MSstats default (used in run_msstats.R)
#   "quantile"              : MSstats quantile normalization
#   "vsn"                   : vsn::justvsn (used in run_nonmsstats.R)
# Each script interprets only the values that make sense for it.
normalization = Sys.getenv("NORMALIZATION", unset = "none")
stopifnot(normalization %in% c("none", "equalizeMedians", "quantile", "vsn"))
stopifnot(variant %in% c("V1_log2", "v2_vsn", "v3_quantile"))
apply_vsn = (normalization == "vsn")

cat(sprintf("[step] report=%s variant=%s suffix='%s' tag='%s' norm=%s exclude=%s\n",
            report_path, variant, out_suffix, out_tag, normalization,
            if (length(exclude_dilutions) > 0)
              paste(exclude_dilutions, collapse = ",") else "(none)"))

vsn_normalize_matrix = function(m) {
  m = as.matrix(m)
  storage.mode(m) = "double"
  dn = dimnames(m)
  # vsn::justvsn rejects NaN; coerce all non-finite to NA before fitting.
  m[!is.finite(m)] = NA_real_
  out = vsn::justvsn(m)
  dimnames(out) = dn
  out
}

# Quantile normalization on a log2-scale matrix. Uses limma's
# normalizeBetweenArrays which handles NAs (it computes per-column ranks
# among available values and interpolates onto a common reference). Input
# is assumed to be log2; output stays on log2.
quantile_normalize_log2_matrix = function(m) {
  m = as.matrix(m)
  storage.mode(m) = "double"
  dn = dimnames(m)
  m[!is.finite(m)] = NA_real_
  out = limma::normalizeBetweenArrays(m, method = "quantile")
  dimnames(out) = dn
  out
}

raw_input  = fread(report_path, sep = "\t")
annotation = fread("CSF_annotation.csv")
raw_input  = raw_input[tolower(raw_input$R.Condition) != "blank", ]
annotation = annotation[tolower(annotation$Condition) != "blank", ]
if (length(exclude_dilutions) > 0) {
  drop_runs = raw_input[R.Condition %in% exclude_dilutions, unique(R.FileName)]
  raw_input  = raw_input[!R.Condition %in% exclude_dilutions, ]
  annotation = annotation[!R.FileName %in% drop_runs, ]
  cat(sprintf("[step] dropped %d runs in dilutions {%s}\n",
              length(drop_runs), paste(exclude_dilutions, collapse = ", ")))
}
annotation$Run = annotation$R.FileName
run_order = unique(annotation[, .(Run, Order)])

merged_input = merge(raw_input, annotation, by = "R.FileName",
                     all.x = TRUE, all.y = FALSE)

protein_swap_list = fread("CSF_protein_swap_list.csv")
true_positives = protein_swap_list[Label == "Positive", Protein]
all_proteins   = protein_swap_list$Protein
no_swap        = character(0)

variant_dir = file.path(out_tag, variant)
out_dir = function(method) {
  d = file.path(variant_dir, paste0(method, out_suffix))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

label_proteins = function(df, protein_col = "Protein") {
  df$Label = ifelse(df[[protein_col]] %in% true_positives,
                     "Positive", "Negative")
  df
}

# Per-method wall-clock timing, split into two phases:
#   preprocess = data prep, summarisation, normalisation
#   model      = parameter estimation + contrast / hypothesis test
# Each method calls tic() at the start of each phase; write_timing writes
# both into <method>_timing.csv.
tic = function() proc.time()[3]
toc = function(t0) as.numeric(proc.time()[3] - t0)

write_timing = function(method, method_dir,
                         preprocess_seconds, model_seconds) {
  fwrite(
    data.frame(method = method,
               preprocess_seconds = preprocess_seconds,
               model_seconds      = model_seconds),
    file = file.path(method_dir, paste0(method, "_timing.csv"))
  )
  cat(sprintf("[time] %-18s %s%s pre=%6.1f s  fit=%6.1f s\n",
              method, variant, out_suffix,
              preprocess_seconds, model_seconds))
  invisible(NULL)
}
