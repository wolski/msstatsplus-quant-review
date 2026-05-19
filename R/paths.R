## paths.R - canonical output path:  <subset_dir>/<normalization>/swap/<package>/
## subset_dir is a folder/subset pair like "CSF_Spectronaut_sample_swap/good_data"
## (relative to the cwd, which is the quant/ root).

out_dir <- function(subset_dir, normalization, package, base_dir = ".") {
  stopifnot(
    is.character(subset_dir), length(subset_dir) == 1L,
    is.character(normalization), length(normalization) == 1L,
    normalization %in% c("log2", "median", "quantile"),
    is.character(package), length(package) == 1L
  )
  d <- file.path(base_dir, subset_dir, normalization, "swap", package)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
