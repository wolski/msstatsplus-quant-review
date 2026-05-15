## MSstats + MSstats+ runner for the original CSF Spectronaut publication
## benchmark. This keeps the original post-dataProcess condition-label swap.

suppressPackageStartupMessages({
  library(MSstats)
  library(MSstatsConvert)
})
source("run_step_common.R")

if (variant != "V1_log2" || normalization != "none") {
  stop("CSF_Spectronaut/run_msstats.R currently supports only ",
       "VARIANT=V1_log2 and NORMALIZATION=none.")
}

true_negatives = protein_swap_list[Label == "Negative", Protein]
comparison = matrix(c(-1, 1), nrow = 1)
rownames(comparison) = "Condition2-Condition1"
colnames(comparison) = c("Condition1", "Condition2")

## MSstats+ --------------------------------------------------------------------
t_pre = tic()
msstats_input = MSstatsConvert::SpectronauttoMSstatsFormat(
  raw_input,
  annotation,
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
fwrite(msstats_input,
       file = file.path(out_dir("MSstats+"), "MSstats+_input.csv"))

summarized = dataProcess(
  msstats_input,
  normalization = FALSE,
  featureSubset = "topN",
  n_top_feature = 100,
  MBimpute = TRUE,
  summaryMethod = "linear",
  numberOfCores = 12
)
save(summarized,
     file = file.path(out_dir("MSstats+"), "MSstats+_summarized.rda"))

weighted_input = summarized$ProteinLevelData
weighted_input$Order = as.integer(stringr::str_split_i(
  weighted_input$originalRUN, "Seq", 2
))
weighted_input = swap_condition_labels(weighted_input,
                                       true_positives,
                                       true_negatives)
summarized$ProteinLevelData = weighted_input
summarized$ProteinLevelData$SUBJECT = as.numeric(as.factor(
  paste0(summarized$ProteinLevelData$originalRUN,
         summarized$ProteinLevelData$GROUP)
))

pre_s = toc(t_pre)
t_mod = tic()
msstatsplus_model = groupComparison(comparison, summarized,
                                    numberOfCores = 12)$ComparisonResult
mod_s = toc(t_mod)

msstatsplus_model = label_proteins(msstatsplus_model)
fwrite(msstatsplus_model,
       file = file.path(out_dir("MSstats+"), "MSstats+_model.csv"))
write_timing("MSstats+", out_dir("MSstats+"), pre_s, mod_s)
message("MSstats+ finished")

## MSstats ---------------------------------------------------------------------
t_pre = tic()
base_msstats_input = MSstatsConvert::SpectronauttoMSstatsFormat(
  raw_input,
  annotation,
  intensity = "PeakArea",
  excludedFromQuantificationFilter = TRUE,
  filter_with_Qvalue = TRUE
)
fwrite(base_msstats_input,
       file = file.path(out_dir("MSstats"), "MSstats_input.csv"))

base_msstats_input = as.data.frame(base_msstats_input) |>
  dplyr::filter(Condition != "Blank")

base_msstats_summarized = dataProcess(
  base_msstats_input,
  normalization = FALSE,
  featureSubset = "topN",
  n_top_feature = 100,
  MBimpute = TRUE,
  summaryMethod = "TMP",
  numberOfCores = 12
)
save(base_msstats_summarized,
     file = file.path(out_dir("MSstats"), "MSstats_summarized.rda"))

weighted_input = base_msstats_summarized$ProteinLevelData
weighted_input$Order = as.integer(stringr::str_split_i(
  weighted_input$originalRUN, "Seq", 2
))
weighted_input = swap_condition_labels(weighted_input,
                                       true_positives,
                                       true_negatives)
weighted_input$Order = NULL
base_msstats_summarized$ProteinLevelData = weighted_input

base_msstats_summarized$ProteinLevelData$Variance = NA
base_msstats_summarized$ProteinLevelData$SUBJECT = as.numeric(as.factor(
  paste0(base_msstats_summarized$ProteinLevelData$originalRUN,
         base_msstats_summarized$ProteinLevelData$GROUP)
))

pre_s = toc(t_pre)
t_mod = tic()
msstats_model = groupComparison(comparison, base_msstats_summarized,
                                numberOfCores = 12)$ComparisonResult
mod_s = toc(t_mod)

msstats_model = label_proteins(msstats_model)
fwrite(msstats_model,
       file = file.path(out_dir("MSstats"), "MSstats_model.csv"))
write_timing("MSstats", out_dir("MSstats"), pre_s, mod_s)
message("MSstats finished")
