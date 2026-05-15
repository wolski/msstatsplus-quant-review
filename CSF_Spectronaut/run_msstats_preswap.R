## MSstats / MSstats+ runner where the condition-label swap is applied before
## dataProcess. This is the pre-swap comparator for the original CSF benchmark.

suppressPackageStartupMessages({
  library(MSstats)
})
source("run_step_common.R")

if (variant != "V1_log2" || normalization != "none") {
  stop("CSF_Spectronaut/run_msstats_preswap.R currently supports only ",
       "VARIANT=V1_log2 and NORMALIZATION=none.")
}

true_negatives = protein_swap_list[Label == "Negative", Protein]
comparison = matrix(c(-1, 1), nrow = 1)
rownames(comparison) = "Condition2-Condition1"
colnames(comparison) = c("Condition1", "Condition2")

run_preswap_variant = function(input_path, out_method, timing_method,
                               out_model, summary_method) {
  if (!file.exists(input_path)) {
    stop("Missing ", input_path, ". Run run_msstats.R first to create the ",
         "MSstats-format input cache.")
  }

  t_pre = tic()
  msstats_input = fread(input_path)
  msstats_input = as.data.frame(msstats_input) |>
    dplyr::filter(Condition != "Blank")

  swapped = swap_condition_labels_msstats(msstats_input,
                                          true_positives,
                                          true_negatives)
  swapped[, BioReplicate := Run]
  swapped = as.data.frame(swapped)

  summarized = dataProcess(
    swapped,
    normalization = FALSE,
    featureSubset = "topN",
    n_top_feature = 100,
    MBimpute = TRUE,
    summaryMethod = summary_method,
    numberOfCores = 12
  )
  save(summarized, file = file.path(out_dir(out_method), "summarized.rda"))

  pre_s = toc(t_pre)
  t_mod = tic()
  model = groupComparison(comparison, summarized,
                          numberOfCores = 12)$ComparisonResult
  mod_s = toc(t_mod)

  model = label_proteins(model)
  fwrite(model, file = file.path(out_dir(out_method), out_model))
  write_timing(timing_method, out_dir(out_method), pre_s, mod_s)
  invisible(model)
}

run_preswap_variant(
  input_path = file.path(out_dir("MSstats+"), "MSstats+_input.csv"),
  out_method = "MSstats+_preswap",
  timing_method = "MSstats+",
  out_model = "MSstats+_preswap_model.csv",
  summary_method = "linear"
)
message("MSstats+_preswap finished")

run_preswap_variant(
  input_path = file.path(out_dir("MSstats"), "MSstats_input.csv"),
  out_method = "MSstats_preswap",
  timing_method = "MSstats",
  out_model = "MSstats_preswap_model.csv",
  summary_method = "TMP"
)
message("MSstats_preswap finished")
