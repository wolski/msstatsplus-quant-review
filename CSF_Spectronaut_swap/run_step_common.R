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
out_tag    = Sys.getenv("OUT_TAG",    unset = "")
exclude_dilutions = Sys.getenv("EXCLUDE_DILUTIONS", unset = "")
exclude_dilutions = if (nchar(exclude_dilutions) > 0) {
  trimws(strsplit(exclude_dilutions, ",", fixed = TRUE)[[1]])
} else {
  character(0)
}
stopifnot(variant %in% c("V1_log2", "v2_vsn"))
apply_vsn = (variant == "v2_vsn")

cat(sprintf("[step] report=%s variant=%s suffix='%s' tag='%s' exclude=%s\n",
            report_path, variant, out_suffix, out_tag,
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

raw_input  = fread(report_path, sep = "\t")
annotation = fread("CSF_annotation.csv")
raw_input  = raw_input[tolower(raw_input$R.Condition) != "blank", ]
annotation = annotation[tolower(annotation$Condition) != "blank", ]
if (length(exclude_dilutions) > 0) {
  raw_input  = raw_input[!R.Condition %in% exclude_dilutions, ]
  annotation = annotation[!R.Condition %in% exclude_dilutions, ]
}
annotation$Run = annotation$R.FileName
run_order = unique(annotation[, .(Run, Order)])

merged_input = merge(raw_input, annotation, by = "R.FileName",
                     all.x = TRUE, all.y = FALSE)

protein_swap_list = fread("CSF_protein_swap_list.csv")
true_positives = protein_swap_list[Label == "Positive", Protein]
all_proteins   = protein_swap_list$Protein
no_swap        = character(0)

variant_dir = paste0(variant, out_tag)
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
