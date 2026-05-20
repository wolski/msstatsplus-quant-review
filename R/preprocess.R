## preprocess.R - Spectronaut TSV reader and per-package input shaping.
##
## Lifted from benchmark_experiments_functions.R. The in-script run/condition
## swap that lived inside each prepare_data_for_* in the original code is
## REMOVED here: the shared pipeline reads pre-swapped TSVs produced by the
## Python helpers in quant/src/ (protein swap: swap_spectronaut_report.py;
## sample swap: swap_spectronaut_report_samples.py - §2). The author-faithful
## CSF_Spectronaut scripts continue to call the originals from
## benchmark_experiments_functions.R.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

## Read the raw Spectronaut TSV and merge with annotation. Drops Blanks and
## any conditions in exclude_dilutions. Adds Order from the "...SeqNN" suffix
## in R.FileName, matching the convention used across the original scripts.
read_spectronaut <- function(report_path, annotation_path,
                              exclude_dilutions = character()) {
  raw <- fread(report_path, sep = "\t")
  annotation <- fread(annotation_path)
  raw <- raw[tolower(R.Condition) != "blank"]
  annotation <- annotation[tolower(Condition) != "blank"]
  if (length(exclude_dilutions) > 0) {
    drop_runs <- raw[R.Condition %in% exclude_dilutions, unique(R.FileName)]
    raw <- raw[!R.Condition %in% exclude_dilutions]
    annotation <- annotation[!R.FileName %in% drop_runs]
  }
  merged <- merge(raw, annotation, by = "R.FileName", all.x = FALSE)
  merged[, Order := as.integer(str_split_i(R.FileName, "Seq", 2))]
  merged
}

.filter_spectronaut <- function(input) {
  input %>% filter(
    F.ExcludedFromQuantification == FALSE &
      !grepl("Blank", R.Condition) &
      F.PeakArea > 1 &
      PG.Qvalue < .01 &
      EG.Qvalue < .01
  )
}

## Fragment-level long table for msqrob2. Output columns: ProteinName,
## Fragment, Run, F.PeakArea.
prepare_data_for_msqrob <- function(input) {
  d <- .filter_spectronaut(input) %>%
    select(PG.ProteinGroups, EG.PrecursorId, F.FrgIon,
           FG.Charge, F.Charge, R.FileName, F.PeakArea)
  d$Fragment <- paste(d$EG.PrecursorId, d$F.Charge, d$F.FrgIon,
                      d$FG.Charge, sep = "_")
  d <- d %>%
    group_by(PG.ProteinGroups, R.FileName, Fragment) %>%
    dplyr::summarize(F.PeakArea = max(F.PeakArea), .groups = "drop") %>%
    filter(PG.ProteinGroups != "") %>%
    setnames(c("PG.ProteinGroups", "R.FileName"), c("ProteinName", "Run"))
  as.data.table(as.data.frame(d)[, c("ProteinName", "Fragment", "Run", "F.PeakArea")])
}

## Fragment-level long table for limma / MaxLFQ. Output keeps all the
## Spectronaut columns plus a Feature key.
prepare_data_for_limma <- function(input) {
  d <- .filter_spectronaut(input)
  d$Feature <- paste(d$EG.PrecursorId, d$F.FrgIon, d$F.Charge, sep = "_")
  as.data.table(d)
}

## Wide log2 matrix for limpa. Rows = (ProteinGroup, Feature); columns = runs.
prepare_data_for_limpa <- function(input) {
  d <- .filter_spectronaut(input)
  d$Feature <- paste(d$EG.PrecursorId, d$F.FrgIon, d$F.Charge, sep = "_")
  d$LogIntensities <- log2(d$F.PeakArea)
  d <- d %>%
    group_by(PG.ProteinGroups, Feature, R.FileName) %>%
    dplyr::summarize(LogIntensities = max(LogIntensities), .groups = "drop")
  d <- pivot_wider(d,
                   id_cols = c("PG.ProteinGroups", "Feature"),
                   names_from = "R.FileName",
                   values_from = "LogIntensities")
  as.data.frame(d)
}

## Wide log2 matrix for DEqMS. Same shape as limpa input.
prepare_data_for_deqms <- function(input) {
  d <- .filter_spectronaut(input)
  d[d$F.PeakArea < 1, "F.PeakArea"] <- NA
  d$Feature <- paste(d$EG.PrecursorId, d$F.FrgIon, d$F.Charge, sep = "_")
  d$LogIntensities <- log2(d$F.PeakArea)
  d <- d %>%
    select(PG.ProteinGroups, Feature, R.FileName, LogIntensities) %>%
    group_by(PG.ProteinGroups, Feature, R.FileName) %>%
    dplyr::summarize(LogIntensities = max(LogIntensities), .groups = "drop") %>%
    filter(PG.ProteinGroups != "")
  d <- pivot_wider(d,
                   id_cols = c("PG.ProteinGroups", "Feature"),
                   names_from = "R.FileName",
                   values_from = "LogIntensities")
  as.data.frame(d)
}

## MSstats long format via SpectronauttoMSstatsFormat. `plus = TRUE` keeps the
## anomaly-score features used by MSstats+; `plus = FALSE` matches the vanilla
## MSstats path.
prepare_data_for_msstats <- function(report_path, annotation_path,
                                      plus = FALSE,
                                      exclude_dilutions = character()) {
  raw <- fread(report_path, sep = "\t")
  annotation <- fread(annotation_path)
  raw <- raw[tolower(R.Condition) != "blank"]
  annotation <- annotation[tolower(Condition) != "blank"]
  if (length(exclude_dilutions) > 0) {
    drop_runs <- raw[R.Condition %in% exclude_dilutions, unique(R.FileName)]
    raw <- raw[!R.Condition %in% exclude_dilutions]
    annotation <- annotation[!R.FileName %in% drop_runs]
  }
  anomaly_features <- if (plus) {
    c("FGShapeQualityScoreMS2", "FGShapeQualityScoreMS1", "EGDeltaRT")
  } else {
    NULL
  }
  MSstats::SpectronauttoMSstatsFormat(
    input = raw,
    annotation = annotation,
    removeFewMeasurements = TRUE,
    use_log_file = FALSE,
    features = anomaly_features
  )
}
