required_packages = c("data.table", "tidyverse")

missing_packages = required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

library(data.table)
library(tidyverse)

source("../benchmark_experiments_functions.R")

prepare_feature_peakarea_matrix = function(input,
                                           sample_proteins_sig,
                                           sample_proteins_insig) {
  feature_input = input %>%
    filter(
      F.ExcludedFromQuantification == FALSE &
        !grepl("Blank", R.Condition) &
        F.PeakArea > 1 &
        PG.Qvalue < .01 &
        EG.Qvalue < .01
    )

  sig = feature_input[(feature_input$PG.ProteinGroups %in%
                         sample_proteins_sig),]
  insig = feature_input[(feature_input$PG.ProteinGroups %in%
                           sample_proteins_insig),]

  mixup_run = feature_input %>%
    distinct(R.FileName, Condition, Order) %>%
    filter(Condition != "Blank") %>%
    arrange(Order)

  cond1_idx = which(mixup_run$Condition == "Condition1")
  cond2_idx = which(mixup_run$Condition == "Condition2")

  flip1 = cond1_idx[seq(2, length(cond1_idx), by = 2)]
  flip2 = cond2_idx[seq(2, length(cond2_idx), by = 2)]

  mixup_run$Condition[flip1] = "Condition2"
  mixup_run$Condition[flip2] = "Condition1"

  mixup_run$New_run = mixup_run$R.FileName
  r1 = mixup_run$R.FileName[flip1]
  r2 = mixup_run$R.FileName[flip2]
  mixup_run$New_run[flip2] = r1
  mixup_run$New_run[flip1] = r2

  insig = insig %>%
    select(-Condition, -Order) %>%
    merge(
      as.data.frame(mixup_run),
      by = "R.FileName",
      all.x = TRUE,
      all.y = FALSE
    )

  insig$R.FileName = insig$New_run
  insig$New_run = NULL

  feature_input = rbindlist(list(sig, insig), use.names = TRUE)
  feature_input$Feature = paste(
    feature_input$EG.PrecursorId,
    feature_input$F.FrgIon,
    feature_input$F.Charge,
    sep = "_"
  )

  feature_input %>%
    filter(F.ExcludedFromQuantification == FALSE) %>%
    select(PG.ProteinGroups, Feature, R.FileName, F.PeakArea) %>%
    group_by(PG.ProteinGroups, Feature, R.FileName) %>%
    dplyr::summarize(F.PeakArea = max(F.PeakArea), .groups = "drop") %>%
    filter(PG.ProteinGroups != "") %>%
    pivot_wider(
      id_cols = c("PG.ProteinGroups", "Feature"),
      names_from = "R.FileName",
      values_from = "F.PeakArea"
    ) %>%
    as.data.frame()
}

data_folder = "."
output_dir = file.path(data_folder, "nonmsstats_preprocessed")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

data_file = "20250130_163144_CSF dilutions Jan 2025 no normalization_Report.tsv"
raw_input = fread(file.path(data_folder, data_file), sep = "\t")
annotation = fread(file.path(data_folder, "CSF_annotation.csv"))

raw_input = raw_input[tolower(raw_input$R.Condition) != "blank",]
annotation = annotation[tolower(annotation$Condition) != "blank",]
annotation$Run = annotation$R.FileName

merged_input = merge(
  raw_input,
  annotation,
  by = "R.FileName",
  all.x = TRUE,
  all.y = FALSE
)

protein_swap_list = fread(file = "CSF_protein_swap_list.csv")
true_positives = unlist(
  protein_swap_list[protein_swap_list$Label == "Positive", "Protein"]
)
true_negatives = unlist(
  protein_swap_list[protein_swap_list$Label == "Negative", "Protein"]
)

truth_labels = rbindlist(list(
  data.table(Protein = true_positives, Label = "Positive"),
  data.table(Protein = true_negatives, Label = "Negative")
))

msqrob_input = prepare_data_for_msqrob(
  merged_input,
  true_positives,
  true_negatives
)
msqrob_wide = dcast(
  msqrob_input,
  ProteinName + Fragment ~ Run,
  value.var = "F.PeakArea",
  fun.aggregate = max,
  fill = NA
)

limma_input = prepare_data_for_limma(
  merged_input,
  true_positives,
  true_negatives
)

feature_peakarea = prepare_feature_peakarea_matrix(
  merged_input,
  true_positives,
  true_negatives
)

fwrite(annotation, file.path(output_dir, "annotation.csv"))
fwrite(truth_labels, file.path(output_dir, "truth_labels.csv"))
fwrite(msqrob_wide, file.path(output_dir, "msqrob_feature_peakarea_wide.csv"))
fwrite(limma_input, file.path(output_dir, "limma_feature_peakarea_long.csv"))
fwrite(feature_peakarea, file.path(output_dir, "feature_peakarea_wide.csv"))

writeLines(
  c(
    "# CSF Spectronaut non-MSstats preprocessing",
    "",
    "This folder contains raw, filtered, benchmark-swapped intermediate inputs for non-MSstats reruns.",
    "",
    "No models are fitted here.",
    "No MSstats or MSstats+ outputs are created or modified here.",
    "",
    "Files:",
    "- `annotation.csv`: non-blank sample annotations.",
    "- `truth_labels.csv`: positive/negative protein labels.",
    "- `msqrob_feature_peakarea_wide.csv`: raw feature peak-area matrix for msqrob2.",
    "- `limma_feature_peakarea_long.csv`: long raw feature table for MaxLFQ + limma.",
    "- `feature_peakarea_wide.csv`: raw feature peak-area matrix for limpa and DEqMS."
  ),
  file.path(output_dir, "README.md")
)

print(data.table(
  file = c(
    "annotation.csv",
    "truth_labels.csv",
    "msqrob_feature_peakarea_wide.csv",
    "limma_feature_peakarea_long.csv",
    "feature_peakarea_wide.csv"
  ),
  rows = c(
    nrow(annotation),
    nrow(truth_labels),
    nrow(msqrob_wide),
    nrow(limma_input),
    nrow(feature_peakarea)
  )
))
