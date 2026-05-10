## MSstats / MSstats+ variants where the condition-label swap is applied to
## the long-format input BEFORE dataProcess, matching the stage at which
## msqrob2/limma/limpa/DEqMS swap. Reuses the cached *_input.csv files written
## by CSF_Spectronaut_processing.R.

library(data.table)
library(tidyverse)
library(MSstats)

source("../benchmark_experiments_functions.R")

data_folder = "."

protein_swap_list = fread(file = "CSF_protein_swap_list.csv")
true_positives = unlist(
  protein_swap_list[protein_swap_list$Label == "Positive", "Protein"])
true_negatives = unlist(
  protein_swap_list[protein_swap_list$Label == "Negative", "Protein"])

comparison = matrix(c(-1, 1), nrow = 1)
row.names(comparison) = "Condition2-Condition1"
colnames(comparison) = c("Condition1", "Condition2")

run_preswap_variant = function(input_path,
                               out_dir,
                               out_model,
                               summary_method) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  msstats_input = fread(input_path)
  msstats_input = as.data.frame(msstats_input) %>% filter(Condition != "Blank")

  swapped = swap_condition_labels_msstats(msstats_input,
                                          true_positives,
                                          true_negatives)
  # After swap a single Condition can contain runs that previously shared a
  # BioReplicate index in their old Condition. Make BioReplicate unique per
  # Run so dataProcess does not collapse them.
  swapped[, BioReplicate := Run]
  swapped = as.data.frame(swapped)

  summarized = dataProcess(swapped,
                           normalization = FALSE,
                           featureSubset = "topN",
                           n_top_feature = 100,
                           MBimpute = TRUE,
                           summaryMethod = summary_method,
                           numberOfCores = 12)

  save(summarized, file = file.path(out_dir, "summarized.rda"))

  model = groupComparison(comparison, summarized, numberOfCores = 12)
  res = model$ComparisonResult
  res$Label = ifelse(res$Protein %in% true_positives, "Positive", "Negative")

  fwrite(res, file = file.path(out_dir, out_model))
  invisible(res)
}

# MSstats+_preswap (linear summarization, anomaly weights already in input)
run_preswap_variant(
  input_path    = file.path("MSstats+", "MSstats+_input.csv"),
  out_dir       = "MSstats+_preswap",
  out_model     = "MSstats+_preswap_model.csv",
  summary_method = "linear"
)
print("MSstats+_preswap finished")

# MSstats_preswap (TMP summarization)
run_preswap_variant(
  input_path    = file.path("MSstats", "MSstats_input.csv"),
  out_dir       = "MSstats_preswap",
  out_model     = "MSstats_preswap_model.csv",
  summary_method = "TMP"
)
print("MSstats_preswap finished")
