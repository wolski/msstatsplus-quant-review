## figures.R - shared figure helpers for diagnostics.qmd, swap_visualization.qmd
## and review.qmd. Each function returns the ggplot object AND writes a PNG to
## out_path (NULL out_path = don't write).
##
## Minimal viable set for §1: p-value histogram, FDR curve, density. Heatmap
## and effect-vs-variance figures are wrapped around prolfqua's plotters.

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

.write_png <- function(p, out_path, width = 7, height = 5) {
  if (!is.null(out_path)) {
    ggplot2::ggsave(out_path, plot = p, width = width, height = height, dpi = 150)
  }
  invisible(p)
}

fig_pvalue_histogram <- function(model_df, out_path = NULL,
                                  title = "p-value distribution") {
  p <- ggplot(model_df, aes(x = pvalue)) +
    geom_histogram(binwidth = 0.025, boundary = 0, fill = "#3b6fb6", color = "white") +
    facet_wrap(~ Package, scales = "free_y") +
    labs(x = "p-value", y = "count", title = title) +
    theme_bw()
  .write_png(p, out_path)
}

fig_fdr_curve <- function(model_df, out_path = NULL) {
  stopifnot(all(c("adj.pvalue", "Label", "Package") %in% names(model_df)))
  dt <- as.data.table(model_df)[!is.na(adj.pvalue) & Label %in% c("Positive","Negative")]
  setorder(dt, Package, adj.pvalue)
  dt[, rank := seq_len(.N), by = Package]
  dt[, n_disc := rank]
  dt[, n_fd := cumsum(Label == "Negative"), by = Package]
  dt[, observed_FDR := n_fd / pmax(n_disc, 1)]
  p <- ggplot(dt, aes(x = adj.pvalue, y = observed_FDR, color = Package)) +
    geom_line() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    coord_cartesian(xlim = c(0, 0.25), ylim = c(0, 0.5)) +
    labs(x = "claimed adj.pvalue", y = "observed FDR") +
    theme_bw()
  .write_png(p, out_path)
}

fig_density <- function(matrix_log2, out_path = NULL,
                         title = "per-sample density (log2 protein matrix)") {
  m <- as.matrix(matrix_log2)
  dt <- data.table(value = as.numeric(m),
                   sample = rep(colnames(m), each = nrow(m)))
  p <- ggplot(dt[is.finite(value)], aes(x = value, group = sample, color = sample)) +
    geom_density(alpha = 0.4) +
    labs(x = "log2 intensity", y = "density", title = title) +
    theme_bw() + theme(legend.position = "none")
  .write_png(p, out_path)
}

## fig_heatmap and fig_na_heatmap delegate to prolfqua when an LFQData object
## is available. For raw matrices we fall back to pheatmap.
fig_heatmap <- function(matrix_log2, annotation = NULL, out_path = NULL) {
  stopifnot(requireNamespace("pheatmap", quietly = TRUE))
  m <- as.matrix(matrix_log2)
  m <- m[rowSums(is.finite(m)) >= ncol(m) * 0.5, , drop = FALSE]
  m <- m[order(-matrixStats::rowSds(m, na.rm = TRUE))[seq_len(min(500, nrow(m)))], , drop = FALSE]
  anno_col <- if (!is.null(annotation)) {
    df <- as.data.frame(annotation[match(colnames(m), annotation$R.FileName), ])
    rownames(df) <- df$R.FileName
    df["Condition"]
  } else NULL
  if (!is.null(out_path)) {
    pheatmap::pheatmap(m, annotation_col = anno_col, filename = out_path,
                       width = 8, height = 10, fontsize_row = 4, show_rownames = FALSE)
  } else {
    pheatmap::pheatmap(m, annotation_col = anno_col, fontsize_row = 4,
                       show_rownames = FALSE)
  }
}

fig_na_heatmap <- function(matrix_log2, annotation = NULL, out_path = NULL) {
  m <- !is.finite(as.matrix(matrix_log2))
  storage.mode(m) <- "integer"
  fig_heatmap(m, annotation = annotation, out_path = out_path)
}

## ---------------------------------------------------------------------------
## Protein-swap per-condition SD density — shared canonical implementation.
## Computed directly from the Spectronaut Report.tsv (NOT from the MSstats
## summary RDA), so reviewers see exactly the values the swap script writes
## without any imputation / TMP summarization. Used by both review.qmd
## (`protein-swap-sd`) and review_supplement.qmd
## (`protein-swap-effect-variance`) so the two figures are byte-identical.
##
## Inputs:
##   report_path     - good_data-subset Spectronaut report (e.g.
##                     CSF_Spectronaut_protein_swap/good_data/Report.tsv)
##   annotation_path - matching annotation.csv with R.FileName + Condition
##                     (Condition1 / Condition2 = G1 / G2 of the swap design)
##   swap_list_path  - CSF_protein_swap_list.csv with Protein + Label
##                     (Positive / Negative)
##
## Returns a ggplot. The same input feeds the supplement's
## `condition_sd_stats_for_review` + `plot_effect_variance` helpers, so
## differences would be purely cosmetic.
## ---------------------------------------------------------------------------

protein_swap_sd_density_from_report <- function(report_path,
                                                  annotation_path,
                                                  swap_list_path,
                                                  title_prefix = NULL) {
  stopifnot(file.exists(report_path),
            file.exists(annotation_path),
            file.exists(swap_list_path))

  rep <- data.table::fread(
    report_path,
    select = c("R.FileName", "PG.ProteinGroups",
               "EG.PrecursorId", "FG.Quantity")
  )
  rep <- rep[nzchar(PG.ProteinGroups) &
               is.finite(FG.Quantity) & FG.Quantity > 0]
  rep <- unique(rep, by = c("R.FileName", "PG.ProteinGroups",
                              "EG.PrecursorId"))
  pld <- rep[, .(LogIntensity = mean(log2(FG.Quantity), na.rm = TRUE)),
                by = .(R.FileName, Protein = PG.ProteinGroups)]

  ann <- data.table::fread(annotation_path)[, .(R.FileName, Condition)]
  pld <- merge(pld, ann, by = "R.FileName")

  swap_list <- data.table::fread(swap_list_path)[, .(Protein, Label)]
  pld <- merge(pld, swap_list, by = "Protein")

  per_cond <- pld[
    ,
    .(sd_log = stats::sd(LogIntensity, na.rm = TRUE), n = .N),
    by = .(Protein, Label, Condition)
  ][is.finite(sd_log) & n >= 2]

  pal <- c(Positive = "#2ca02c", Negative = "#7f7f7f")
  title <- if (is.null(title_prefix)) {
    "Within-condition SD (per-condition, direct from Report.tsv)"
  } else {
    paste(title_prefix, "within-condition SD")
  }
  ggplot2::ggplot(per_cond, ggplot2::aes(x = sd_log, fill = Label)) +
    ggplot2::geom_density(alpha = 0.55) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::labs(
      title = title,
      x = "sd(log2 intensity) within one condition",
      y = "density",
      fill = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top")
}
