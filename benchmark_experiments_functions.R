library(tidyverse)
library(data.table)

#' Swap condition labels for a subset of runs
#'
#' This function introduces controlled label swapping between two conditions
#' ("Condition1" and "Condition2") at the run level, while preserving the
#' original labels for a subset of proteins. Specifically, every second run
#' within each condition is reassigned to the opposite condition. The swapping
#' is applied only to proteins specified as non-significant, while significant
#' proteins retain their original condition labels.
#'
#' @param input A data.frame or data.table containing at least the columns:
#'   \code{originalRUN}, \code{GROUP}, \code{Order}, and \code{Protein}.
#' @param sample_proteins_sig A character vector of protein identifiers that
#'   should retain their original condition labels (treated as "signal").
#' @param sample_proteins_insig A character vector of protein identifiers that
#'   will have their condition labels swapped (treated as "noise").
#'
#' @return A data.table with the same structure as \code{input}, where condition
#'   labels (\code{GROUP}) have been selectively swapped for non-significant
#'   proteins according to the alternating run pattern.
swap_condition_labels = function(input, sample_proteins_sig, sample_proteins_insig){

  mixup_run = input %>% 
    distinct(originalRUN, GROUP, Order) %>% arrange(Order)
  
  # Get the indices of each condition
  cond1_idx = which(mixup_run$GROUP == "Condition1")
  cond2_idx = which(mixup_run$GROUP == "Condition2")
  
  # Flip every other instance
  flip1 = cond1_idx[seq(2, length(cond1_idx), by = 2)]
  flip2 = cond2_idx[seq(2, length(cond2_idx), by = 2)]
  
  # Swap values
  mixup_run$GROUP[flip1] = "Condition2"
  mixup_run$GROUP[flip2] = "Condition1"
  
  sig = input[(input$Protein %in% sample_proteins_sig),]
  insig = input[(input$Protein %in% sample_proteins_insig),]
  insig = insig %>% dplyr::select(-GROUP, -Order) %>% merge(
    as.data.frame(mixup_run), by="originalRUN", all.x=TRUE, all.y=FALSE
  )
  input = rbindlist(list(sig, insig), use.names=TRUE)

  return(input)
}

#' Swap condition labels on MSstats long-format input (pre-dataProcess)
#'
#' Equivalent to \code{swap_condition_labels} but operates on the
#' SpectronauttoMSstatsFormat / DIANNtoMSstatsFormat long table, before
#' summarization. Required columns: \code{Run}, \code{Condition},
#' \code{ProteinName}. \code{Order} is derived from \code{Run} as
#' \code{as.integer(str_split_i(Run, "Seq", 2))} to match the convention used
#' elsewhere in this benchmark.
#'
#' @param input data.frame/data.table in MSstats long format.
#' @param sample_proteins_sig Proteins that retain their original Condition.
#' @param sample_proteins_insig Proteins whose Condition is swapped on every
#'   second run within each condition.
swap_condition_labels_msstats = function(input,
                                          sample_proteins_sig,
                                          sample_proteins_insig){
  input = as.data.table(input)

  run_table = unique(input[, .(Run, Condition)])
  run_table[, Order := as.integer(str_split_i(Run, "Seq", 2))]
  setorder(run_table, Order)

  cond1_idx = which(run_table$Condition == "Condition1")
  cond2_idx = which(run_table$Condition == "Condition2")
  flip1 = cond1_idx[seq(2, length(cond1_idx), by = 2)]
  flip2 = cond2_idx[seq(2, length(cond2_idx), by = 2)]

  r1 = run_table$Run[flip1]
  r2 = run_table$Run[flip2]

  # Swap the Run identifier itself for paired runs (mirrors the msqrob2
  # variant in prepare_data_for_msqrob). Each Run name then still maps to
  # exactly one Condition, so dataProcess sees consistent (Run, Condition)
  # pairs. Condition is read from the original (Run -> Condition) table.
  partner = setNames(run_table$Run, run_table$Run)
  partner[r1] = r2
  partner[r2] = r1

  cond_for_run = setNames(run_table$Condition, run_table$Run)

  sig = input[ProteinName %in% sample_proteins_sig]
  insig = input[ProteinName %in% sample_proteins_insig]
  insig[, Run := partner[Run]]
  insig[, Condition := cond_for_run[Run]]

  rbindlist(list(sig, insig), use.names = TRUE)
}

#' Prepare benchmark data for MSqRob analysis with controlled label/run swapping
#'
#' Filters and reshapes benchmark  data into fragment-level input for MSqRob,
#' while introducing structured label and run swapping for a subset of proteins.
#'
#' @param input A data.frame or data.table containing DIA-NN/Spectronaut-style
#'   evidence data with required columns (e.g., \code{R.FileName},
#'   \code{Condition}, \code{PG.ProteinGroups}, \code{F.PeakArea}, etc.).
#' @param sample_proteins_sig A character vector of protein identifiers that
#'   should retain their original condition labels (treated as "signal").
#' @param sample_proteins_insig A character vector of protein identifiers that
#'   will have their condition labels swapped (treated as "noise").
#'
#' @return A data.table with columns \code{ProteinName}, \code{Fragment},
#'   \code{Run}, and \code{F.PeakArea}, suitable for MSqRob input.
prepare_data_for_msqrob = function(input, sample_proteins_sig, sample_proteins_insig){
  
  evidence = copy(input)
  evidence = evidence %>% filter(R.Condition != "Blank")
  evidence$Order = as.integer(
    str_split_i(evidence$R.FileName, "Seq", 2))
  
  # Mix up runs
  mixup_run = evidence %>%
    distinct(R.FileName, Condition, Order) %>%
    arrange(Order)

  # Get the indices of each condition
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

  sig = evidence[(evidence$PG.ProteinGroups %in% sample_proteins_sig), ]
  insig = evidence[(evidence$PG.ProteinGroups %in% sample_proteins_insig), ]

  insig = insig %>% select(-Condition, -Order) %>% merge(
    as.data.frame(mixup_run), by="R.FileName", all.x=TRUE, all.y=FALSE
  )
  insig$R.FileName = insig$New_run
  insig$New_run = NULL
  evidence = rbindlist(list(sig, insig), use.names=TRUE)
  
  # Apply similar filtering as in MSstats
  input_data = evidence %>% 
    filter(F.ExcludedFromQuantification == FALSE & 
             !grepl("Blank", R.Condition) &
             F.PeakArea > 1 & 
             PG.Qvalue < .01 &
             EG.Qvalue < .01) %>% 
    select(PG.ProteinGroups, EG.PrecursorId, F.FrgIon, 
           FG.Charge, F.Charge, R.FileName, F.PeakArea)
  input_data$Fragment = paste(input_data$EG.PrecursorId, 
                              input_data$F.Charge,
                              input_data$F.FrgIon, 
                              input_data$FG.Charge, sep="_")
  input_data = input_data %>% group_by(PG.ProteinGroups, R.FileName, Fragment) %>% 
    dplyr::summarize(F.PeakArea=max(F.PeakArea))
  input_data = input_data %>% filter(PG.ProteinGroups != "")
  input_data = setnames(input_data, 
                        c("PG.ProteinGroups", "R.FileName"),
                        c("ProteinName", "Run"))
  
  df_LFQ = as.data.frame(input_data) %>% 
    select(ProteinName, Fragment, Run, F.PeakArea)
  df_LFQ = as.data.table(df_LFQ)
  
  return(df_LFQ)
}

#' Prepare benchmark proteomics data for limma with structured label/run swapping
#'
#' Filters and formats fragment-level proteomics data for limma analysis,
#' while introducing controlled condition and run swapping for non-significant
#' proteins.
#'
#' @param input A data.frame or data.table containing proteomics evidence data.
#' @param sample_proteins_sig Character vector of proteins treated as signal
#'   (labels unchanged).
#' @param sample_proteins_insig Character vector of proteins treated as noise
#'   (labels and runs may be swapped).
#'
#' @return A data.table with fragment-level features and updated run/condition
#'   labels, suitable for downstream summarization and limma analysis.
prepare_data_for_limma = function(input, sample_proteins_sig, sample_proteins_insig){
  
  # Similar filtering to MSstats
  limma_input = input %>% filter(
    F.ExcludedFromQuantification == FALSE & 
      !grepl("Blank", R.Condition) &
      F.PeakArea > 1 & 
      PG.Qvalue < .01 &
      EG.Qvalue < .01)
  
  sig = limma_input[(limma_input$PG.ProteinGroups %in% sample_proteins_sig), ]
  insig = limma_input[(limma_input$PG.ProteinGroups %in% sample_proteins_insig), ]

  # Mix up runs
  mixup_run = limma_input %>%
    distinct(R.FileName, Condition, Order) %>%
    filter(Condition != "Blank") %>%
    arrange(Order)

  # Get the indices of each condition
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

  insig = insig %>% select(-Condition, -Order) %>% merge(
    as.data.frame(mixup_run), by="R.FileName",
    all.x=TRUE, all.y=FALSE)

  insig$R.FileName = insig$New_run
  insig$New_run = NULL

  limma_input = rbindlist(list(sig, insig), use.names=TRUE)
  
  # Prep data for summarization
  limma_input$Feature = paste(limma_input$EG.PrecursorId, 
                              limma_input$F.FrgIon,
                              limma_input$F.Charge, sep="_")
  
  return(limma_input)
}

#' Prepare benchmark proteomics data for limpa with randomized label/run swapping
#'
#' Filters and reshapes proteomics data into a wide, log-intensity matrix for
#' limpa, while introducing randomized condition and run swapping for a subset
#' of non-significant proteins.
#'
#' @param input A data.frame or data.table containing proteomics evidence data.
#' @param sample_proteins_sig Character vector of proteins treated as signal
#'   (labels unchanged).
#' @param sample_proteins_insig Character vector of proteins treated as noise
#'   (labels and runs may be swapped).
#'   
#' @return A data.frame in wide format with proteins and features as rows and
#'   runs as columns, containing log2-transformed intensities.
prepare_data_for_limpa = function(input, sample_proteins_sig, sample_proteins_insig){
  
  # Similar filtering to MSstats
  limpa_input = input %>% filter(
    F.ExcludedFromQuantification == FALSE & 
      !grepl("Blank", R.Condition) &
      F.PeakArea > 1 & 
      PG.Qvalue < .01 &
      EG.Qvalue < .01)
  
  sig = limpa_input[(limpa_input$PG.ProteinGroups %in% sample_proteins_sig), ]
  insig = limpa_input[(limpa_input$PG.ProteinGroups %in% sample_proteins_insig), ]

  # Mix up runs
  mixup_run = limpa_input %>%
    distinct(R.FileName, Condition, Order) %>%
    filter(Condition != "Blank") %>%
    arrange(Order)

  # Get the indices of each condition
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

  insig = insig %>% select(-Condition, -Order) %>% merge(
    as.data.frame(mixup_run), by="R.FileName",
    all.x=TRUE, all.y=FALSE)

  insig$R.FileName = insig$New_run
  insig$New_run = NULL

  limpa_input = rbindlist(list(sig, insig), use.names=TRUE)
  
  # Prep data for summarization
  limpa_input$Feature = paste(limpa_input$EG.PrecursorId, 
                              limpa_input$F.FrgIon,
                              limpa_input$F.Charge, sep="_")
  
  limpa_input = limpa_input %>% 
    filter(F.ExcludedFromQuantification == FALSE)
  
  limpa_input$LogIntensities = log2(limpa_input$F.PeakArea)
  
  limpa_input = limpa_input %>% 
    group_by(PG.ProteinGroups, Feature, R.FileName) %>% 
    summarize(LogIntensities=max(LogIntensities))
  
  limpa_input = pivot_wider(limpa_input,
                            id_cols=c("PG.ProteinGroups", "Feature"),
                            names_from="R.FileName",
                            values_from="LogIntensities")
  
  limpa_input = as.data.frame(limpa_input)
  
  return(limpa_input)
}

#' Prepare benchmark proteomics data for DEqMS with structured label/run swapping
#'
#' Filters and reshapes proteomics data into a wide log-intensity matrix for
#' DEqMS, while introducing alternating condition and run swaps for
#' non-significant proteins.
#'
#' @param input A data.frame or data.table containing proteomics evidence data.
#' @param sample_proteins_sig Character vector of proteins treated as signal
#'   (labels unchanged).
#' @param sample_proteins_insig Character vector of proteins treated as noise
#'   (labels and runs may be swapped).
#'
#' @return A data.frame in wide format with proteins and features as rows and
#'   runs as columns, containing log2-transformed intensities.
prepare_data_for_deqms = function(input, sample_proteins_sig, sample_proteins_insig){
  
  deqms_input = input %>% filter(
    F.ExcludedFromQuantification == FALSE & 
      !grepl("Blank", R.Condition) &
      F.PeakArea > 1 & 
      PG.Qvalue < .01 &
      EG.Qvalue < .01)
  deqms_input[deqms_input$F.PeakArea < 1, "F.PeakArea"] = NA
  
  sig = deqms_input[(deqms_input$PG.ProteinGroups %in% sample_proteins_sig), ]
  insig = deqms_input[(deqms_input$PG.ProteinGroups %in% sample_proteins_insig), ]

  # Mix up runs
  mixup_run = deqms_input %>%
    distinct(R.FileName, Condition, Order) %>%
    filter(Condition != "Blank") %>%
    arrange(Order)

  # Get the indices of each condition
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

  insig = insig %>% select(-Condition, -Order) %>% merge(
    as.data.frame(mixup_run), by="R.FileName",
    all.x=TRUE, all.y=FALSE)

  insig$R.FileName = insig$New_run
  insig$New_run = NULL

  deqms_input = rbindlist(list(sig, insig), use.names=TRUE)
  
  # Prep data for summarization
  deqms_input$Feature = paste(deqms_input$EG.PrecursorId, 
                              deqms_input$F.FrgIon,
                              deqms_input$F.Charge, sep="_")
  
  deqms_input = deqms_input %>% 
    filter(F.ExcludedFromQuantification == FALSE) %>% 
    select(PG.ProteinGroups,Feature, R.FileName, F.PeakArea)
  deqms_input$LogIntensities = log2(deqms_input$F.PeakArea)
  
  deqms_input = deqms_input %>% group_by(PG.ProteinGroups, Feature, R.FileName) %>% 
    dplyr::summarize(LogIntensities=max(LogIntensities))
  deqms_input = deqms_input %>% filter(PG.ProteinGroups != "")
  
  deqms_input = pivot_wider(deqms_input,
                            id_cols=c("PG.ProteinGroups", "Feature"),
                            names_from="R.FileName",
                            values_from="LogIntensities")
  return(deqms_input)
}

#' Prepare DIA-NN protein-level data for limma
#'
#' Reshapes DIA-NN protein-level output into a log2-transformed wide matrix for
#' limma analysis, with optional condition and run swapping for a subset of
#' non-significant proteins.
#'
#' @param input A data.frame or data.table containing DIA-NN protein-level output.
#' @param annotation A data.frame with run-level annotation, including run,
#'   condition, and order information.
#' @param swap_conditions Logical; whether to swap condition labels and run
#'   identities for non-significant proteins.
#' @param sample_proteins_sig Optional character vector of proteins treated as
#'   signal (unchanged).
#' @param sample_proteins_insig Optional character vector of proteins treated as
#'   noise (may be swapped).
#'
#' @return A wide data.frame with proteins as rows, runs as columns, and
#'   log2-transformed intensities.
prepare_diann_data_for_limma = function(input,
                                        annotation,
                                        swap_conditions=FALSE, 
                                        sample_proteins_sig=NULL, 
                                        sample_proteins_insig=NULL){
  
  colnames(input) = basename(colnames(input))
  
  # Identify run columns (those ending in .mzML)
  run_cols = grep("\\.mzML$", names(input), value = TRUE)
  
  input = melt(
    input,
    id.vars = setdiff(names(input), run_cols),
    measure.vars = run_cols,
    variable.name = "Run",
    value.name = "Intensity"
  )
  input[, Run := sub("\\.mzML$", "", Run)]
  
  input = input %>% merge(annotation, by="Run", all.x=TRUE, all.y=FALSE)
  
  if (swap_conditions){
    sig = input[(input$Protein.Group %in% sample_proteins_sig), ]
    insig = input[(input$Protein.Group %in% sample_proteins_insig), ]
    
    # Mix up runs
    mixup_run = input %>% 
      distinct(R.FileName, Condition, Order) %>% 
      filter(Condition != "Blank") %>% 
      arrange(Order)
    
    # Get the indices of each condition
    cond1_idx = which(mixup_run$Condition == "Condition1")
    cond2_idx = which(mixup_run$Condition == "Condition2")
    
    flip1 = sample(cond1_idx[1:15], 10)
    flip2 = sample(cond2_idx[1:15], 10)
    
    mixup_run$Condition[flip1] = "Condition2"
    mixup_run$Condition[flip2] = "Condition1"
    
    mixup_run$New_run = mixup_run$R.FileName
    r1 = mixup_run$R.FileName[flip1]
    r2 = mixup_run$R.FileName[flip2]
    mixup_run$New_run[flip2] = r1
    mixup_run$New_run[flip1] = r2
    
    insig = insig %>% select(-Condition, -Order) %>% merge(
      as.data.frame(mixup_run), by="R.FileName", 
      all.x=TRUE, all.y=FALSE)
    
    insig$Run = insig$New_run
    insig$New_run = NULL
    
    input = rbindlist(list(sig, insig), use.names=TRUE)
  }
  
  input$Intensity = log2(input$Intensity)
  
  input = input %>% pivot_wider(id_cols="Protein.Group",
                                names_from="Run",
                                values_from="Intensity",
                                values_fn=max) %>% as.data.frame()
  
  input= as.data.frame(input)
  rownames(input) = input$Protein.Group
  input$Protein.Group = NULL
  input$Protein.Names = NULL
  input$Genes = NULL
  input$First.Protein.Description  = NULL
  input$N.Sequences = NULL
  input$N.Proteotypic.Sequences = NULL
  
  return(input)
}
