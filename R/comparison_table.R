## comparison_table.R - Table 1 (TPR/PPV) generalised from
## CSF_Spectronaut_swap_comparison_table.R. Works on the canonical
## 6-column model schema written by run_<pkg>(): Protein, logFC, SE, DF,
## pvalue, adj.pvalue.

suppressPackageStartupMessages(library(data.table))

# label_proteins() comes from R/ground_truth.R. Callers (review.qmd,
# diagnostics.qmd, run_cell.R helpers) source ground_truth.R themselves
# before this file, so we don't re-source it here — that previous
# unconditional source() broke whenever cwd != quant/ (e.g. when this
# file was loaded from vignettes/).

.canonical_packages <- c("MSstats+", "MSstats", "limpa", "MaxLFQ_limma",
                          "msqrob2", "DEqMS", "prolfqua")

.model_filename <- list(
  "MSstats+"     = "MSstats+_model.csv",
  "MSstats"      = "MSstats_model.csv",
  "MaxLFQ_limma" = "limma_model.csv",
  "limpa"        = "limpa_model.csv",
  "msqrob2"      = "msqrob2_model.csv",
  "DEqMS"        = "deqms_model.csv",
  "prolfqua"     = "prolfqua_model.csv"
)

tpr_ppv <- function(p, fc, label, threshold = 0.05) {
  finite_fc <- is.finite(fc)
  hit <- !is.na(p) & p < threshold & finite_fc
  n_disc <- sum(hit, na.rm = TRUE)
  false_disc <- sum(hit & label == "Negative", na.rm = TRUE)
  list(
    TPR = sum(hit & label == "Positive", na.rm = TRUE) /
          sum(finite_fc & label == "Positive", na.rm = TRUE),
    PPV = if (n_disc > 0) 1 - false_disc / n_disc else NA_real_
  )
}

## build_table1: for one (csffolder, dataset, normalization) cell, iterate
## over the 7 packages and compute TPR/PPV at p<0.05 and adj.pvalue<0.05.
## `cell_dir` is the parent of the per-package dirs (i.e.,
## <csffolder>/<dataset>/<normalization>/swap/).
build_table1 <- function(cell_dir, truth_df,
                          threshold_p = 0.05, threshold_fdr = 0.05) {
  rows <- lapply(.canonical_packages, function(pkg) {
    f <- file.path(cell_dir, pkg, .model_filename[[pkg]])
    if (!file.exists(f)) {
      return(data.table(Package = pkg,
                        TPR_pval = NA_real_, PPV_pval = NA_real_,
                        TPR_FDR  = NA_real_, PPV_FDR  = NA_real_))
    }
    m <- fread(f)
    m <- as.data.frame(label_proteins(as.data.frame(m), truth_df))
    # MSstats/MSstats+ write log2FC; the other packages write logFC. Normalise.
    fc <- if ("logFC" %in% names(m)) m$logFC else m$log2FC
    p_metrics <- tpr_ppv(m$pvalue, fc, m$Label, threshold_p)
    f_metrics <- tpr_ppv(m$adj.pvalue, fc, m$Label, threshold_fdr)
    data.table(Package = pkg,
               TPR_pval = p_metrics$TPR, PPV_pval = p_metrics$PPV,
               TPR_FDR  = f_metrics$TPR, PPV_FDR  = f_metrics$PPV)
  })
  rbindlist(rows)
}

write_table1 <- function(tbl, out_prefix) {
  fwrite(tbl, file = paste0(out_prefix, ".csv"))
  md <- c(
    paste("|", paste(colnames(tbl), collapse = " | "), "|"),
    paste("|", paste(rep("---", ncol(tbl)), collapse = " | "), "|"),
    apply(tbl, 1, function(row) paste("|", paste(row, collapse = " | "), "|"))
  )
  writeLines(md, paste0(out_prefix, ".md"))
  invisible(NULL)
}
