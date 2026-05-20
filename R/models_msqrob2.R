## models_msqrob2.R - msqrob2 hurdle on protein-aggregated QFeatures.
##
## Normalization is now strictly label-driven:
##   log2     -> no normalization (raw log2 peptides aggregated to protein)
##   median   -> peptide-level center.median (matches msqrob2 vignette)
##   quantile -> protein-level quantile (limma::normalizeBetweenArrays)
## Previously the peptide-level center.median ran unconditionally, which made
## log2 results look like median-centered results.

source("R/normalize.R")
source("R/preprocess.R")
source("R/timing.R")

run_msqrob2 <- function(merged_input, annotation, normalization, out_path) {
  suppressPackageStartupMessages({
    library(QFeatures); library(msqrob2); library(prolfqua); library(data.table)
  })

  t_pre <- tic()
  ms2_input <- as.data.frame(prepare_data_for_limma(merged_input))
  ms2_input <- ms2_input[is.finite(ms2_input$F.PeakArea) &
                           ms2_input$F.PeakArea > 0, , drop = FALSE]

  config <- prolfqua::AnalysisConfiguration$new()
  config$file_name <- "R.FileName"
  config$factors["group_"] <- "Condition"
  config$hierarchy[["protein_Id"]]   <- "PG.ProteinGroups"
  config$hierarchy[["precursor_Id"]] <- "EG.PrecursorId"
  config$hierarchy[["fragment_Id"]]  <- c("Feature", "F.FrgLossType")
  config$hierarchy_depth <- 1
  config$set_response("F.PeakArea")
  adata <- prolfqua::setup_analysis(ms2_input, config)
  lfq_pep <- prolfqua::LFQData$new(adata, config)

  tr_pep <- lfq_pep$get_Transformer(); tr_pep$log2()
  lfq_pep_log <- tr_pep$lfq

  se <- prolfqua::LFQDataToSummarizedExperiment(lfqdata = lfq_pep_log)
  pe <- QFeatures::QFeatures(list(peptide = se),
                              colData = SummarizedExperiment::colData(se))

  # Peptide-level normalization applied ONLY for the "median" label, matching
  # the msqrob2 vignette. log2 / quantile aggregate raw log2 peptides.
  if (normalization == "median") {
    pe <- QFeatures::normalize(pe, i = "peptide", method = "center.median",
                                name = "peptide_norm")
    agg_from <- "peptide_norm"
  } else {
    agg_from <- "peptide"
  }
  pe <- QFeatures::aggregateFeatures(
    pe, i = agg_from, fcol = "protein_Id", name = "protein"
  )

  # Protein-level normalization. log2 = identity; median already done at
  # peptide level (so identity here); quantile applied here.
  if (normalization == "quantile") {
    prot_mat <- SummarizedExperiment::assay(pe[["protein"]])
    prot_mat <- apply_normalization(prot_mat, normalization)
    SummarizedExperiment::assay(pe[["protein"]]) <- prot_mat
  }
  pre_s <- toc(t_pre)

  t_mod <- tic()
  prlm <- msqrob2::msqrobHurdle(pe, i = "protein",
                                 formula = ~ 0 + group_, overwrite = TRUE)
  L <- msqrob2::makeContrast(
    "group_Condition2 - group_Condition1=0",
    parameterNames = c("group_Condition1", "group_Condition2")
  )
  prlm <- msqrob2::hypothesisTestHurdle(prlm, i = "protein", L,
                                         overwrite = TRUE)
  mod_s <- toc(t_mod)
  save(prlm, file = file.path(out_path, "msqrob_obj.rda"))

  xx <- SummarizedExperiment::rowData(prlm[["protein"]])
  hurdle_cols <- grep("^hurdle_", names(xx), value = TRUE)
  stopifnot(length(hurdle_cols) == 1)
  hdf <- as.data.frame(xx[[hurdle_cols[1]]])
  hdf$Protein <- rownames(xx)

  use_intensity <- !is.na(hdf$logFC)
  model <- data.frame(
    Protein = c(hdf$Protein[use_intensity], hdf$Protein[!use_intensity]),
    logFC   = c(hdf$logFC[use_intensity],   hdf$logOR[!use_intensity]),
    SE      = c(hdf$logFCse[use_intensity], hdf$logORse[!use_intensity]),
    DF      = c(hdf$logFCdf[use_intensity], hdf$logORdf[!use_intensity]),
    pvalue  = c(hdf$logFCpval[use_intensity], hdf$logORpval[!use_intensity]),
    source  = c(rep("intensity", sum(use_intensity)),
                rep("count",     sum(!use_intensity))),
    stringsAsFactors = FALSE
  )
  model$adj.pvalue <- p.adjust(model$pvalue, method = "BH")

  data.table::fwrite(model, file = file.path(out_path, "msqrob2_model.csv"))
  write_timing("msqrob2", out_path, pre_s, mod_s)
  list(model = model, preprocess_seconds = pre_s, model_seconds = mod_s)
}
