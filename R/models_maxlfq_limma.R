## models_maxlfq_limma.R - MaxLFQ summarization + limma DE.

source("R/normalize.R")
source("R/preprocess.R")
source("R/timing.R")

run_maxlfq_limma <- function(merged_input, annotation, normalization, out_path) {
  suppressPackageStartupMessages({
    library(iq); library(limma); library(data.table)
  })

  t_pre <- tic()
  limma_long <- prepare_data_for_limma(merged_input)
  maxlfq_input <- iq::preprocess(
    limma_long,
    primary_id = "PG.ProteinGroups",
    secondary_id = c("Feature"),
    sample_id = "R.FileName",
    intensity_col = "F.PeakArea",
    median_normalization = FALSE,
    log2_intensity_cutoff = 0,
    pdf_out = file.path(out_path, "qc-plots.pdf"),
    pdf_width = 12, pdf_height = 8,
    intensity_col_sep = NULL,
    intensity_col_id = NULL,
    na_string = "0"
  )
  maxlfq_summarized <- iq::fast_MaxLFQ(maxlfq_input)$estimate
  maxlfq_summarized <- apply_normalization(maxlfq_summarized, normalization)

  class <- annotation$Condition[
    match(colnames(maxlfq_summarized), annotation$R.FileName)
  ] |> as.factor()
  design <- model.matrix(~ 0 + class)
  pre_s <- toc(t_pre)

  t_mod <- tic()
  fit1 <- limma::lmFit(maxlfq_summarized, design = design)
  cont <- limma::makeContrasts(classCondition2 - classCondition1, levels = design)
  fit2 <- limma::contrasts.fit(fit1, contrasts = cont)
  fit3 <- limma::eBayes(fit2)
  mod_s <- toc(t_mod)

  model <- data.frame(
    Protein    = rownames(fit3$coefficients),
    logFC      = as.numeric(fit3$coefficients),
    SE         = as.numeric(sqrt(fit3$s2.post) * fit3$stdev.unscaled),
    DF         = fit3$df.total,
    pvalue     = as.numeric(fit3$p.value)
  )
  model$adj.pvalue <- p.adjust(model$pvalue, method = "BH")

  data.table::fwrite(model, file = file.path(out_path, "limma_model.csv"))
  write_timing("limma", out_path, pre_s, mod_s)
  list(model = model, preprocess_seconds = pre_s, model_seconds = mod_s)
}
