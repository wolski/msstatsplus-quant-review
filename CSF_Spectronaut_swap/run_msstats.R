## MSstats + MSstats+ runner.
##
## Shares SpectronauttoMSstatsFormat between the two methods. They differ in
## summaryMethod ("linear" vs "TMP") and in whether the MSstats+ anomaly
## scoring is enabled.
##
## Env vars: see run_step_common.R for the full list. NORMALIZATION values
## that this script understands:
##   "none"             -> dataProcess(normalization = FALSE)
##   "equalizeMedians"  -> dataProcess(normalization = "equalizeMedians")
##   "quantile"         -> dataProcess(normalization = "quantile")
## "vsn" is rejected here (use run_nonmsstats.R for vsn-based normalization).
suppressPackageStartupMessages({
  library(MSstats)
  library(MSstatsConvert)
})
source("run_step_common.R")

if (normalization == "vsn") {
  stop("run_msstats.R does not support NORMALIZATION=vsn. ",
       "Use 'equalizeMedians' or 'quantile' to apply MSstats's built-in ",
       "normalization, or 'none' to disable it.")
}
msstats_norm = if (normalization == "none") FALSE else normalization
cat(sprintf("[msstats] dataProcess normalization = %s\n",
            if (isFALSE(msstats_norm)) "FALSE" else sprintf("'%s'", msstats_norm)))

comparison = matrix(c(-1, 1), nrow = 1)
rownames(comparison) = "Condition2-Condition1"
colnames(comparison) = c("Condition1", "Condition2")

## MSstats+ --------------------------------------------------------------------
msstats_input = MSstatsConvert::SpectronauttoMSstatsFormat(
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
  numberOfCores = 12
)
fwrite(msstats_input, file = file.path(out_dir("MSstats+"),
                                          "MSstats+_input.csv"))

summarized = dataProcess(
  msstats_input,
  normalization = msstats_norm,
  featureSubset = "topN",
  n_top_feature = 100,
  MBimpute = TRUE,
  summaryMethod = "linear",
  numberOfCores = 12
)
save(summarized, file = file.path(out_dir("MSstats+"),
                                    "MSstats+_summarized.rda"))

summarized$ProteinLevelData$SUBJECT = as.numeric(as.factor(
  paste0(summarized$ProteinLevelData$originalRUN,
         summarized$ProteinLevelData$GROUP)))

msstatsplus_model = groupComparison(comparison, summarized,
                                     numberOfCores = 12)$ComparisonResult
msstatsplus_model = label_proteins(msstatsplus_model)
fwrite(msstatsplus_model, file = file.path(out_dir("MSstats+"),
                                              "MSstats+_model.csv"))
message("MSstats+ finished")

## MSstats ---------------------------------------------------------------------
## Re-run the Spectronaut->MSstats format conversion *without* the anomaly
## scoring (the only difference at the input level).
base_msstats_input = MSstatsConvert::SpectronauttoMSstatsFormat(
  raw_input, annotation,
  intensity = "PeakArea",
  excludedFromQuantificationFilter = TRUE,
  filter_with_Qvalue = TRUE
)
fwrite(base_msstats_input, file = file.path(out_dir("MSstats"),
                                                "MSstats_input.csv"))

base_msstats_input = as.data.frame(base_msstats_input) |>
  dplyr::filter(Condition != "Blank")

base_msstats_summarized = dataProcess(
  base_msstats_input,
  normalization = msstats_norm,
  featureSubset = "topN",
  n_top_feature = 100,
  MBimpute = TRUE,
  summaryMethod = "TMP",
  numberOfCores = 12
)
save(base_msstats_summarized,
     file = file.path(out_dir("MSstats"), "MSstats_summarized.rda"))

base_msstats_summarized$ProteinLevelData$Variance = NA
base_msstats_summarized$ProteinLevelData$SUBJECT = as.numeric(as.factor(
  paste0(base_msstats_summarized$ProteinLevelData$originalRUN,
         base_msstats_summarized$ProteinLevelData$GROUP)))

msstats_model = groupComparison(comparison, base_msstats_summarized,
                                 numberOfCores = 12)$ComparisonResult
msstats_model = label_proteins(msstats_model)
fwrite(msstats_model, file = file.path(out_dir("MSstats"),
                                          "MSstats_model.csv"))
message("MSstats finished")
