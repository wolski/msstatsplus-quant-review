## models_limpa.R - limpa dpc / dpcQuant / dpcDE pipeline.
## Author-faithful tryCatch around dpc() retained: limpa can fail to converge
## on hard datasets, in which case we write an empty model and continue.

source("R/normalize.R")
source("R/preprocess.R")
source("R/timing.R")

run_limpa <- function(merged_input, annotation, normalization, out_path) {
  suppressPackageStartupMessages({
    library(limpa); library(limma); library(data.table)
  })

  t_pre <- tic()
  limpa_input_wide <- prepare_data_for_limpa(merged_input)
  mapper <- limpa_input_wide[c("PG.ProteinGroups", "Feature")]
  rownames(mapper) <- mapper$Feature

  limpa_dt <- limpa_input_wide
  rownames(limpa_dt) <- limpa_dt$Feature
  limpa_dt <- limpa_dt[, !colnames(limpa_dt) %in% c("PG.ProteinGroups", "Feature")]
  # limpa_input_wide is already on the log2 scale (prepare_data_for_limpa
  # applies log2). For log2 label: pass through. For median/quantile: apply
  # the matching helper directly.
  if (normalization == "median") {
    limpa_dt <- as.data.frame(median_normalize_log2_matrix(as.matrix(limpa_dt)))
  } else if (normalization == "quantile") {
    limpa_dt <- as.data.frame(quantile_normalize_log2_matrix(as.matrix(limpa_dt)))
  }

  targets <- as.data.frame(annotation)[c("R.FileName", "Condition")]
  rownames(targets) <- targets$R.FileName

  limpa_elist <- methods::new(
    "EList", list(E = limpa_dt, genes = mapper, targets = targets)
  )
  pre_s <- toc(t_pre)

  t_mod <- tic()
  dpcfit <- tryCatch(limpa::dpc(limpa_elist),
                      error = function(e) {
                        warning("limpa dpc() failed under normalization=",
                                normalization, ": ", conditionMessage(e))
                        NULL
                      })
  if (is.null(dpcfit)) {
    mod_s <- toc(t_mod)
    write_timing("limpa", out_path, pre_s, mod_s)
    return(list(model = NULL, preprocess_seconds = pre_s, model_seconds = mod_s))
  }
  y.protein <- limpa::dpcQuant(limpa_elist, "PG.ProteinGroups", dpc = dpcfit)
  class <- annotation$Condition[
    match(colnames(y.protein$E), annotation$R.FileName)
  ] |> as.factor()
  design <- model.matrix(~ 0 + class)
  fit <- limpa::dpcDE(y.protein, design, plot = FALSE)
  cont <- limma::makeContrasts(classCondition2 - classCondition1, levels = design)
  fit <- limma::contrasts.fit(fit, contrasts = cont)
  fit <- limma::eBayes(fit)
  mod_s <- toc(t_mod)

  model <- data.frame(
    Protein    = rownames(fit$coefficients),
    logFC      = as.numeric(fit$coefficients),
    SE         = as.numeric(sqrt(fit$s2.post) * fit$stdev.unscaled),
    DF         = fit$df.total,
    pvalue     = as.numeric(fit$p.value)
  )
  model$adj.pvalue <- p.adjust(model$pvalue, method = "BH")

  data.table::fwrite(model, file = file.path(out_path, "limpa_model.csv"))
  write_timing("limpa", out_path, pre_s, mod_s)
  list(model = model, preprocess_seconds = pre_s, model_seconds = mod_s)
}
