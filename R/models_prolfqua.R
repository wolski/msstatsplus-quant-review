## models_prolfqua.R - prolfqua medpolish + ContrastsLMImputeFacade.

source("R/normalize.R")
source("R/preprocess.R")
source("R/timing.R")

run_prolfqua <- function(merged_input, annotation, normalization, out_path) {
  suppressPackageStartupMessages({
    library(prolfqua); library(data.table)
  })

  t_pre <- tic()
  prolfqua_input <- as.data.frame(prepare_data_for_limma(merged_input))
  needed <- c("PG.ProteinGroups", "EG.PrecursorId", "Feature", "R.FileName",
              "Condition", "F.PeakArea")
  stopifnot(all(needed %in% colnames(prolfqua_input)))
  prolfqua_input <- prolfqua_input[is.finite(prolfqua_input$F.PeakArea) &
                                     prolfqua_input$F.PeakArea > 0, , drop = FALSE]

  config <- prolfqua::AnalysisConfiguration$new()
  config$file_name <- "R.FileName"
  config$factors["group_"] <- "Condition"
  config$hierarchy[["protein_Id"]]   <- "PG.ProteinGroups"
  config$hierarchy[["precursor_Id"]] <- "EG.PrecursorId"
  config$hierarchy[["fragment_Id"]]  <- c("Feature", "F.FrgLossType")
  config$hierarchy_depth <- 1
  config$set_response("F.PeakArea")
  adata <- prolfqua::setup_analysis(prolfqua_input, config)
  lfqdata <- prolfqua::LFQData$new(adata, config)

  tr_log <- lfqdata$get_Transformer()$intensity_array(log)
  agg <- tr_log$lfq$get_Aggregator("medpolish"); agg$aggregate()
  lfq_protein_log <- agg$lfq_agg
  tr_inv <- lfq_protein_log$get_Transformer()$intensity_array(exp, force = TRUE)
  lfq_protein <- tr_inv$lfq
  lfq_protein$is_transformed(FALSE)

  tr_norm <- lfq_protein$get_Transformer()
  tr_norm$log2()
  if (normalization == "quantile") {
    tr_norm$intensity_matrix(.func = quantile_normalize_log2_matrix,
                              force = TRUE)
  } else if (normalization == "median") {
    tr_norm$intensity_matrix(.func = median_normalize_log2_matrix,
                              force = TRUE)
  }
  lfq_protein <- tr_norm$lfq
  pre_s <- as.numeric(toc(t_pre))

  t_mod <- tic()
  contr_spec <- c("Condition2_vs_Condition1" =
                    "group_Condition2 - group_Condition1")
  fa <- prolfqua::ContrastsLMImputeFacade$new(lfq_protein, "~ group_",
                                               contr_spec)
  res <- fa$get_contrasts()
  mod_s <- as.numeric(toc(t_mod))

  model <- data.frame(
    Protein    = res$protein_Id,
    logFC      = res$diff,
    SE         = res$std.error,
    DF         = res$df,
    pvalue     = res$p.value,
    adj.pvalue = res$FDR,
    stringsAsFactors = FALSE
  )
  data.table::fwrite(model, file = file.path(out_path, "prolfqua_model.csv"))
  write_timing("prolfqua", out_path, pre_s, mod_s)
  list(model = model, preprocess_seconds = pre_s, model_seconds = mod_s)
}
