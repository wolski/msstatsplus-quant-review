## models_deqms.R - DEqMS (limma + spectraCounteBayes for pep-count variance).

source("R/normalize.R")
source("R/preprocess.R")
source("R/timing.R")

.summarize_deqms_no_ref_col <- function(dat, group_col = 1) {
  dat_ratio <- dat
  dat_ratio[, 3:ncol(dat)] <- dat_ratio[, 3:ncol(dat)] -
    matrixStats::rowMedians(as.matrix(dat_ratio[, 3:ncol(dat)]), na.rm = TRUE)
  dat_summary <- plyr::ddply(
    dat_ratio, colnames(dat)[group_col],
    function(x) matrixStats::colMedians(as.matrix(x[, 3:ncol(dat)]), na.rm = TRUE)
  )
  colnames(dat_summary)[2:ncol(dat_summary)] <- colnames(dat)[3:ncol(dat)]
  out <- dat_summary[, -1]
  rownames(out) <- dat_summary[, 1]
  out
}

run_deqms <- function(merged_input, annotation, normalization, out_path) {
  suppressPackageStartupMessages({
    library(limma); library(DEqMS); library(data.table)
  })

  t_pre <- tic()
  shared_input <- prepare_data_for_deqms(merged_input)
  deqms_summarized <- .summarize_deqms_no_ref_col(shared_input, group_col = 1)
  deqms_summarized <- as.data.frame(
    apply_normalization(as.matrix(deqms_summarized), normalization)
  )

  pep_count <- shared_input |>
    dplyr::group_by(PG.ProteinGroups) |>
    dplyr::summarise(count = dplyr::n_distinct(Feature) + 1, .groups = "drop") |>
    as.data.frame()
  rownames(pep_count) <- pep_count$PG.ProteinGroups
  pep_count$PG.ProteinGroups <- NULL

  class <- annotation$Condition[
    match(colnames(deqms_summarized), annotation$R.FileName)
  ] |> as.factor()
  design <- model.matrix(~ 0 + class)
  pre_s <- toc(t_pre)

  t_mod <- tic()
  fit1 <- limma::lmFit(deqms_summarized, design = design)
  cont <- limma::makeContrasts(classCondition2 - classCondition1, levels = design)
  fit2 <- limma::contrasts.fit(fit1, contrasts = cont)
  fit3 <- limma::eBayes(fit2)
  fit3$count <- pep_count[rownames(fit3$coefficients), "count"]
  fit4 <- DEqMS::spectraCounteBayes(fit3)
  mod_s <- toc(t_mod)

  model <- data.frame(
    Protein    = rownames(fit4$coefficients),
    logFC      = as.numeric(fit4$coefficients),
    SE         = as.numeric(sqrt(fit4$s2.post) * fit4$stdev.unscaled),
    DF         = fit4$df.total,
    pvalue     = as.numeric(fit4$p.value)
  )
  model$adj.pvalue <- p.adjust(model$pvalue, method = "BH")

  data.table::fwrite(model, file = file.path(out_path, "deqms_model.csv"))
  data.table::fwrite(
    data.frame(Protein = rownames(pep_count), count = pep_count$count),
    file = file.path(out_path, "deqms_pep_count.csv")
  )
  write_timing("DEqMS", out_path, pre_s, mod_s)
  list(model = model, preprocess_seconds = pre_s, model_seconds = mod_s)
}
