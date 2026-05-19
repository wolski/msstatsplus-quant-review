## run_cell.R - execute one (subset_dir, normalization, package) cell.
##
## Usage:
##   Rscript R/run_cell.R <subset_dir> <normalization> <package>
##
## <subset_dir>  is a path like "CSF_Spectronaut_sample_swap/good_data"
##               (relative to quant/) containing Report.tsv + annotation.csv
##               pre-filtered to the desired subset (by src/build_subsets.py).
## <normalization> is one of log2, median, quantile.
## <package>     is one of MSstats, MSstats+, MaxLFQ_limma, DEqMS, msqrob2,
##               prolfqua, limpa.
##
## Output:       <subset_dir>/<normalization>/swap/<package>/<package>_model.csv

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript R/run_cell.R <subset_dir> <normalization> <package>")
}
subset_dir    <- args[[1]]
normalization <- args[[2]]
package       <- args[[3]]

suppressPackageStartupMessages(library(data.table))
source("R/paths.R")
source("R/preprocess.R")

report_path     <- file.path(subset_dir, "Report.tsv")
annotation_path <- file.path(subset_dir, "annotation.csv")
if (!file.exists(report_path))     stop("Not found: ", report_path)
if (!file.exists(annotation_path)) stop("Not found: ", annotation_path)

raw_input  <- fread(report_path, sep = "\t")
annotation <- fread(annotation_path)
raw_input  <- raw_input[tolower(R.Condition) != "blank"]
annotation <- annotation[tolower(Condition) != "blank"]
raw_input  <- raw_input[R.FileName %in% annotation$R.FileName]
merged_input <- merge(raw_input, annotation, by = "R.FileName", all.x = FALSE)
merged_input[, Order := as.integer(stringr::str_split_i(R.FileName, "Seq", 2))]

cell <- out_dir(subset_dir, normalization, package)

if (package %in% c("MSstats", "MSstats+")) {
  source("R/models_msstats.R")
  run_msstats(raw_input, annotation, normalization, cell,
              plus = identical(package, "MSstats+"))
} else if (package == "MaxLFQ_limma") {
  source("R/models_maxlfq_limma.R")
  run_maxlfq_limma(merged_input, annotation, normalization, cell)
} else if (package == "DEqMS") {
  source("R/models_deqms.R");  run_deqms(merged_input, annotation, normalization, cell)
} else if (package == "msqrob2") {
  source("R/models_msqrob2.R"); run_msqrob2(merged_input, annotation, normalization, cell)
} else if (package == "prolfqua") {
  source("R/models_prolfqua.R"); run_prolfqua(merged_input, annotation, normalization, cell)
} else if (package == "limpa") {
  source("R/models_limpa.R");   run_limpa(merged_input, annotation, normalization, cell)
} else {
  stop("Unknown package: ", package)
}
