## models_msstats.R - MSstats and MSstats+ runner.
##
## plus = TRUE  -> MSstats+ (anomaly scoring enabled, summaryMethod = "linear")
## plus = FALSE -> MSstats   (no anomaly scoring,        summaryMethod = "TMP")
##
## Both branches respect the canonical normalization labels:
##   "log2"     -> dataProcess(normalization = FALSE)
##   "median"   -> dataProcess(normalization = "equalizeMedians")
##   "quantile" -> dataProcess(normalization = "quantile")

source("R/normalize.R")
source("R/timing.R")

run_msstats <- function(raw_input, annotation, normalization, out_path,
                         plus = FALSE, n_cores = 4L) {
  suppressPackageStartupMessages({
    library(MSstats); library(MSstatsConvert); library(data.table)
  })
  msstats_norm <- msstats_normalization_arg(normalization)
  method_name <- if (plus) "MSstats+" else "MSstats"

  comparison <- matrix(c(-1, 1), nrow = 1)
  rownames(comparison) <- "Condition2-Condition1"
  colnames(comparison) <- c("Condition1", "Condition2")

  # SpectronauttoMSstatsFormat requires a `Run` column in the annotation;
  # the upstream CSF/annotation files name it `R.FileName`.
  if (!"Run" %in% colnames(annotation)) {
    annotation <- as.data.frame(annotation)
    annotation$Run <- annotation$R.FileName
  }

  t_pre <- tic()
  if (plus) {
    # Anomaly-score features need a (Run, Order) lookup table.
    run_order <- unique(as.data.table(annotation)[, .(Run, Order)])
    msstats_input <- MSstatsConvert::SpectronauttoMSstatsFormat(
      raw_input, annotation,
      intensity = "PeakArea",
      excludedFromQuantificationFilter = TRUE,
      filter_with_Qvalue = TRUE,
      calculateAnomalyScores = TRUE,
      anomalyModelFeatures = c("FGShapeQualityScore(MS2)",
                                "FGShapeQualityScore(MS1)",
                                "EGDeltaRT"),
      anomalyModelFeatureTemporal = c("mean_decrease",
                                       "mean_decrease",
                                       "dispersion_increase"),
      removeMissingFeatures = .75,
      runOrder = run_order,
      max_depth = "auto",
      numberOfCores = n_cores
    )
  } else {
    msstats_input <- MSstatsConvert::SpectronauttoMSstatsFormat(
      raw_input, annotation,
      intensity = "PeakArea",
      excludedFromQuantificationFilter = TRUE,
      filter_with_Qvalue = TRUE
    )
    msstats_input <- as.data.frame(msstats_input) |>
      dplyr::filter(Condition != "Blank")
  }
  data.table::fwrite(msstats_input,
                     file = file.path(out_path, paste0(method_name, "_input.csv")))

  summarized <- MSstats::dataProcess(
    msstats_input,
    normalization = msstats_norm,
    featureSubset = "topN",
    n_top_feature = 100,
    MBimpute = TRUE,
    summaryMethod = if (plus) "linear" else "TMP",
    numberOfCores = n_cores
  )
  save(summarized,
       file = file.path(out_path, paste0(method_name, "_summarized.rda")))

  pld <- summarized$ProteinLevelData
  if (!plus) pld$Variance <- NA
  pld$SUBJECT <- as.numeric(as.factor(paste0(pld$originalRUN, pld$GROUP)))
  summarized$ProteinLevelData <- pld
  pre_s <- toc(t_pre)

  t_mod <- tic()
  cmp <- MSstats::groupComparison(comparison, summarized,
                                   numberOfCores = n_cores)$ComparisonResult
  mod_s <- toc(t_mod)

  model <- data.frame(
    Protein    = cmp$Protein,
    logFC      = cmp$log2FC,
    SE         = cmp$SE,
    DF         = cmp$DF,
    pvalue     = cmp$pvalue,
    adj.pvalue = cmp$adj.pvalue,
    stringsAsFactors = FALSE
  )
  data.table::fwrite(model,
                     file = file.path(out_path,
                                       paste0(method_name, "_model.csv")))
  write_timing(method_name, out_path, pre_s, mod_s)
  list(model = model, preprocess_seconds = pre_s, model_seconds = mod_s)
}
