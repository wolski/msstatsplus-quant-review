## normalize.R - canonical normalization for the shared pipeline. See
## R/README_normalization.md (§0.5 audit) for semantics. All operations are
## centering-only or distribution-matching; no per-sample scaling.

## NA-safe per-column median centering on a log2-scale matrix.
median_normalize_log2_matrix <- function(m) {
  m <- as.matrix(m)
  storage.mode(m) <- "double"
  dn <- dimnames(m)
  m[!is.finite(m)] <- NA_real_
  cm <- matrixStats::colMedians(m, na.rm = TRUE)
  out <- sweep(m, 2, cm, FUN = "-")
  dimnames(out) <- dn
  out
}

## Quantile normalization via limma::normalizeBetweenArrays.
quantile_normalize_log2_matrix <- function(m) {
  m <- as.matrix(m)
  storage.mode(m) <- "double"
  dn <- dimnames(m)
  m[!is.finite(m)] <- NA_real_
  out <- limma::normalizeBetweenArrays(m, method = "quantile")
  dimnames(out) <- dn
  out
}

## Dispatch by label for non-MSstats packages. `log2` is identity here
## because the input is assumed already on log2 scale.
apply_normalization <- function(m, label) {
  stopifnot(label %in% c("log2", "median", "quantile"))
  switch(label,
    log2 = m,
    median = median_normalize_log2_matrix(m),
    quantile = quantile_normalize_log2_matrix(m)
  )
}

## Map our shared label onto the argument MSstats::dataProcess expects.
msstats_normalization_arg <- function(label) {
  stopifnot(label %in% c("log2", "median", "quantile"))
  switch(label,
    log2 = FALSE,
    median = "equalizeMedians",
    quantile = "quantile"
  )
}
