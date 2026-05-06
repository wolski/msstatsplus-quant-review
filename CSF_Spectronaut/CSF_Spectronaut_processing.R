## Packages --------------------------------------------------------------------
# Data analysis
library(data.table)
library(tidyverse)

# Comparison packages
library(MSstats)
library(MSstatsConvert)
library(QFeatures)
library(msqrob2)
library(iq)
library(limma)
library(limpa)
library(DEqMS)

# Functions for data conversion/method comparison
source("benchmark_experiments_functions.R")

## Load + prepare inputs -------------------------------------------------------
# Load data
data_folder = ""
data_file = "20250130_163144_CSF dilutions Jan 2025 no normalization_Report.tsv"
raw_input = fread(paste(data_folder, data_file, sep="/"), sep="\t")

annotation = fread(paste(data_folder, "CSF_annotation.csv", sep="/"))

# Prepare annotation info
raw_input = raw_input[raw_input$R.Condition != "Blank",]
annotation = annotation[annotation$Condition != "Blank",]
annotation$Run = annotation$R.FileName
run_order = unique(annotation[, .(Run, Order)])

merged_input = merge(raw_input, annotation, by="R.FileName",
                     all.x=TRUE, all.y=FALSE)

protein_swap_list = fread(file="CSF_protein_swap_list.csv")

true_positives = unlist(
  protein_swap_list[protein_swap_list$Label == "Positive", "Protein"])
true_negatives = unlist(
  protein_swap_list[protein_swap_list$Label == "Negative", "Protein"])

# MSstats+ ---------------------------------------------------------------------
msstats_input = MSstatsConvert::SpectronauttoMSstatsFormat(
  raw_input,
  annotation,
  intensity = 'PeakArea', 
  excludedFromQuantificationFilter = TRUE,
  filter_with_Qvalue = TRUE,
  calculateAnomalyScores=TRUE, 
  anomalyModelFeatures=c("FGShapeQualityScore(MS2)",
                         "FGShapeQualityScore(MS1)",
                         "EGDeltaRT"),
  anomalyModelFeatureTemporal=c("mean_decrease",
                                "mean_decrease",
                                "dispersion_increase"),
  removeMissingFeatures=.75,
  runOrder=run_order,
  max_depth="auto",
  numberOfCores=12)
fwrite(msstats_input, file=paste(data_folder, "MSstats+",
                                 "MSstats+_input.csv", sep="/"))

summarized = dataProcess(msstats_input,
                         normalization=FALSE,
                         featureSubset = "topN",
                         n_top_feature = 100,
                         MBimpute=TRUE,
                         summaryMethod="linear",
                         numberOfCores = 12)

save(summarized, file=paste(data_folder, "MSstats+",
                            "MSstats+_summarized.rda", sep="/"))

# Swap condition labels for TP/TN analysis
weighted_input = summarized$ProteinLevelData
weighted_input$Order = as.integer(str_split_i(
  weighted_input$originalRUN, "Seq", 2))

weighted_input = swap_condition_labels(weighted_input, 
                                       true_positives, 
                                       true_negatives)

summarized$ProteinLevelData = weighted_input
summarized$ProteinLevelData$SUBJECT = as.numeric(as.factor(
  paste0(summarized$ProteinLevelData$originalRUN,
         summarized$ProteinLevelData$GROUP)))

comparison = matrix(c(-1,1),nrow=1)
row.names(comparison) = "Condition2-Condition1"
colnames(comparison) = c("Condition1", "Condition2")

msstatsplus_model = groupComparison(comparison, summarized, 
                                 numberOfCores=12)

msstatsplus_model$ComparisonResult$Label = ifelse(
  msstatsplus_model$ComparisonResult$Protein %in% true_positives,
  "Positive", "Negative")

msstatsplus_model = msstatsplus_model$ComparisonResult
fwrite(msstatsplus_model, file=paste(data_folder, "MSstats+", 
                                     "MSstats+_model_swap.csv", sep="/"))

print("MSstats+ finished")

# MSstats ----------------------------------------------------------------------
base_msstats_input = MSstatsConvert::SpectronauttoMSstatsFormat(
  raw_input, annotation, intensity = 'PeakArea', 
  excludedFromQuantificationFilter = TRUE,
  filter_with_Qvalue = TRUE)
fwrite(base_msstats_input, 
       file=paste(data_folder, "MSstats", "MSstats_input.csv", sep="/"))

base_msstats_input = as.data.frame(base_msstats_input) %>% 
  filter(Condition != "Blank")

base_msstats_summarized = dataProcess(base_msstats_input, 
                                 normalization=FALSE,
                                 featureSubset = "topN",
                                 n_top_feature = 100,
                                 MBimpute=TRUE,
                                 summaryMethod="TMP",
                                 numberOfCores = 12)
save(base_msstats_summarized, 
     file=paste(data_folder, "MSstats", "MSstats_summarized.rda", sep="/"))

# Swap condition labels for TP/TN analysis
weighted_input = base_msstats_summarized$ProteinLevelData
weighted_input$Order = as.integer(
  str_split_i(weighted_input$originalRUN, "Seq", 2))

weighted_input = swap_condition_labels(weighted_input, protein_swap_list)

weighted_input$Order = NULL
base_msstats_summarized$ProteinLevelData = weighted_input

base_msstats_summarized$ProteinLevelData$Variance = NA
base_msstats_summarized$ProteinLevelData$SUBJECT = as.numeric(as.factor(
  paste0(base_msstats_summarized$ProteinLevelData$originalRUN,
         base_msstats_summarized$ProteinLevelData$GROUP)))

levels(base_msstats_summarized$ProteinLevelData$GROUP)
comparison = matrix(c(-1,1),nrow=1)
row.names(comparison) = "Condition2-Condition1"
colnames(comparison) = c("Condition1", "Condition2")

msstats_model = groupComparison(comparison, base_msstats_summarized, 
                                numberOfCores=12)

msstats_model$ComparisonResult$Label = ifelse(
  msstats_model$ComparisonResult$Protein %in% sample_proteins_sig,
  "Positive", "Negative")

msstats_model = msstats_model$ComparisonResult

fwrite(msstats_model, 
       file=paste(data_folder, "MSstats", "MSstats_model.csv", sep="/"))

print("MSstats finished")
# MSqRob ----------------------------------------------------------------------
df_LFQ = prepare_data_for_msqrob(merged_input, true_positives, true_negatives)

# Convert into MSqRob2 format
df_LFQ = dcast(df_LFQ, ProteinName + Fragment ~ Run, 
               value.var = "F.PeakArea", fun.aggregate = max, fill=NA)

# Prepare data for modeling
colnames(df_LFQ)[3:length(colnames(df_LFQ))] = paste(
  "F.PeakArea", colnames(df_LFQ)[3:length(colnames(df_LFQ))], sep="_")

ecols = grep("F.PeakArea", colnames(df_LFQ))

pe = readQFeatures(df_LFQ, fnames = 2, quantCols = ecols,  name = "peptideRaw")

colData(pe)$condition = annotation$Condition[
  match(str_remove(colnames(df_LFQ)[grep("F.PeakArea", colnames(df_LFQ))], 
                   "F.PeakArea_"), annotation$R.FileName)] %>% 
  as.factor()

rowData(pe[["peptideRaw"]])$nNonZero = rowSums(assay(pe[["peptideRaw"]]) > 0, 
                                               na.rm = TRUE)
pe = zeroIsNA(pe, "peptideRaw") # convert 0 to NA

pe = logTransform(pe, base = 2, i = "peptideRaw", name = "peptideLog")

Protein_filter = rowData(
  pe[["peptideLog"]])$ProteinName %in% smallestUniqueGroups(
    rowData(pe[["peptideLog"]])$ProteinName)
pe = pe[Protein_filter,]

pe = filterFeatures(pe, ~ nNonZero >= 2)

pe = QFeatures::impute(pe, i = "peptideLog",
                       name = "peptideImp",
                       method = "QRILC")

# Summarization
pe = aggregateFeatures(pe,
                       i = "peptideImp", fcol = "ProteinName",
                       name = "protein")

# Model
pe = msqrob(object = pe, i = "protein", formula = ~condition)

L = makeContrast("conditionCondition2=0", 
                 parameterNames = c("conditionCondition2",
                                    "conditionCondition1"))
pe = hypothesisTest(object = pe, i = "protein", contrast = L)
save(pe, file=paste(data_folder, "msqrob2", "msqrob_obj.rda", sep="/"))

msqrob2_model = rowData(pe[["protein"]])$`conditionCondition2`

msqrob2_model$Label = ifelse(
  rownames(msqrob2_model) %in% sample_proteins_sig, "Positive", "Negative")

msqrob2_model$Protein = rownames(msqrob2_model)

fwrite(msqrob2_model, file=paste(data_folder, "msqrob2", 
                                 "msqrob2_model.csv", sep="/"))

print("msqrob finished")

# Limma ------------------------------------------------------------------------
limma_input = prepare_data_for_limma(merged_input, true_positives, true_negatives)

maxlfq_input = preprocess(
  limma_input,
  primary_id = "PG.ProteinGroups",
  secondary_id = c("Feature"),
  sample_id = "R.FileName",
  intensity_col = "F.PeakArea",
  median_normalization = FALSE,
  log2_intensity_cutoff = 0,
  pdf_out = "qc-plots.pdf",
  pdf_width = 12,
  pdf_height = 8,
  intensity_col_sep = NULL,
  intensity_col_id = NULL,
  na_string = "0")

maxlfq_summarized = fast_MaxLFQ(maxlfq_input)
maxlfq_summarized = maxlfq_summarized$estimate

class = annotation$Condition[
  match(colnames(maxlfq_summarized), annotation$R.FileName)] %>% 
  as.factor()
design = model.matrix(~0+class) # fitting without intercept

fit1 = lmFit(maxlfq_summarized, design = design)
cont = makeContrasts(classCondition2-classCondition1, levels = design)
fit2 = contrasts.fit(fit1, contrasts = cont)
fit3 = eBayes(fit2)

limma_model = data.frame("Protein" = rownames(fit3$coefficients),
                         "logFC" = as.numeric(fit3$coefficients),
                         "SE" = as.numeric(sqrt(fit3$s2.post) * fit3$stdev.unscaled),
                         "DF" = fit3$df.total,
                         "pvalue" = as.numeric(fit3$p.value)
)

rownames(limma_model) = NULL
limma_model$adj.pvalue = p.adjust(limma_model$pvalue, method = 'BH')

limma_model$Label = ifelse(
  limma_model$Protein %in% true_positives, "Positive", "Negative")

fwrite(limma_model, file=paste(data_folder, "limma", 
                               "limma_model.csv", sep="/"))

print("limma finished")

# Limpa ------------------------------------------------------------------------
limpa_input = prepare_data_for_limpa(merged_input, true_positives, true_negatives)

mapper = limpa_input[c("PG.ProteinGroups", "Feature")]
row.names(mapper) = mapper$Feature
limpa_dt = copy(limpa_input)
row.names(limpa_dt) = limpa_dt$Feature
limpa_dt = limpa_dt[, !colnames(limpa_dt) %in% c("PG.ProteinGroups", "Feature")]

annotation = annotation[annotation$Condition != "Blank"]
targets = as.data.frame(annotation)[c("R.FileName", "Condition")]
row.names(targets) = targets$R.FileName
limpa_elist = new("EList", list(E = limpa_dt,
                                genes = mapper,
                                targets = targets))

dpcfit = dpc(limpa_elist)
dpcfit$dpc
plotDPC(dpcfit)

## dpcQuant
y.protein = dpcQuant(limpa_elist, "PG.ProteinGroups", 
                     dpc = dpcfit)

class = annotation$Condition[
  match(colnames(y.protein$E), annotation$R.FileName)] %>% 
  as.factor()
design = model.matrix(~0+class) # fitting without intercept

## dpcDE
fit = dpcDE(y.protein, design, plot=TRUE)
cont = makeContrasts(classCondition2-classCondition1, levels = design)
fit = contrasts.fit(fit, contrasts = cont)
fit = eBayes(fit)

limpa_model = data.frame("Protein" = rownames(fit$coefficients),
                         "logFC" = as.numeric(fit$coefficients),
                         "SE" = as.numeric(sqrt(fit$s2.post) * fit$stdev.unscaled),
                         "DF" = fit$df.total,
                         "pvalue" = as.numeric(fit$p.value)
)

rownames(limpa_model) = NULL
limpa_model$adj.pvalue = p.adjust(limpa_model$pvalue, method = 'BH')

limpa_model$Label = ifelse(
  limpa_model$Protein %in% sample_proteins_sig, "Positive", "Negative")

fwrite(limpa_model, file=paste(data_folder, "limpa",
                               "limpa_model.csv", sep="/"))

print("limpa finished")
# DEqMS ------------------------------------------------------------------------
deqms_input = prepare_data_for_deqms(merged_input, true_positives, true_negatives)

# Adjusted DEqMS::medianSummary function without ref_col
summarize = function (dat, group_col = 2) 
{
  dat.ratio = dat
  dat.ratio[, 3:ncol(dat)] = dat.ratio[, 3:ncol(dat)] - 
    matrixStats::rowMedians(as.matrix(dat.ratio[, 3:ncol(dat)]), na.rm = TRUE)
  dat.summary = plyr::ddply(dat.ratio, colnames(dat)[group_col], 
                function(x) matrixStats::colMedians(as.matrix(x[, 3:ncol(dat)]), 
                                                    na.rm = TRUE))
  colnames(dat.summary)[2:ncol(dat.summary)] = colnames(dat)[3:ncol(dat)]
  dat.new = dat.summary[, -1]
  rownames(dat.new) = dat.summary[, 1]
  return(dat.new)
}

deqms_summarized = summarize(deqms_input, group_col=1)

pep_count = deqms_input %>% group_by(PG.ProteinGroups)%>%
  summarise(count=n_distinct(Feature),.groups = 'drop') %>%
  as.data.frame()

# Minimum peptide count of some proteins can be 0
# add pseudocount 1 to all proteins
pep_count$count = pep_count$count+1
pep_count = as.data.frame(pep_count)
rownames(pep_count) = pep_count$PG.ProteinGroups
pep_count$PG.ProteinGroups = NULL

annotation = annotation %>% filter(Condition != "Blank")

class = annotation$Condition[
  match(colnames(deqms_summarized), annotation$R.FileName)] %>% 
  as.factor()
design = model.matrix(~0+class) # fitting without intercept

fit1 = lmFit(deqms_summarized, design = design)
cont = makeContrasts(classCondition2-classCondition1, levels = design)
fit2 = contrasts.fit(fit1, contrasts = cont)
fit3 = eBayes(fit2)

fit3$count = pep_count[rownames(fit3$coefficients), "count"]

fit4 = spectraCounteBayes(fit3)

## Analyze results
deqms_model = data.frame("Protein" = rownames(fit4$coefficients),
                         "logFC" = as.numeric(fit4$coefficients),
                         "SE" = as.numeric(sqrt(fit4$s2.post) * fit4$stdev.unscaled),
                         "DF" = fit4$df.total,
                         "pvalue" = as.numeric(fit4$p.value))
rownames(deqms_model) = NULL
deqms_model$adj.pvalue = p.adjust(deqms_model$pvalue, method = 'BH')

deqms_model$Label = ifelse(
  deqms_model$Protein %in% sample_proteins_sig, "Positive", "Negative")

fwrite(deqms_model, file=paste(data_folder, "DEqMS", 
                               "deqms_model.csv", sep="/"))

print("DeqMS finished")
# mapDIA -----------------------------------------------------------------------
# Run externally (outside of R). Here we just prepare the input.
# Use same data processesing as for limma
mapdia_input = prepare_data_for_limma(merged_input, true_positives, true_negatives)
mapdia_input$Fragment = paste(mapdia_input$FG.Charge, 
                              mapdia_input$F.FrgIon,
                              mapdia_input$F.Charge, sep="_")
mapdia_input$Sample_info = paste(mapdia_input$Condition, 
                                 mapdia_input$BioReplicate, sep="_")

rt_data = mapdia_input %>% group_by(EG.ModifiedSequence, Fragment) %>% 
  summarize(mean_rt = mean(EG.ApexRT))

mapdia_input = mapdia_input %>% select(PG.ProteinGroups, EG.ModifiedSequence, 
                                       Fragment, Sample_info, F.PeakArea) %>%
  pivot_wider(
    names_from = Sample_info,
    values_from = F.PeakArea,
    values_fn = max,
    values_fill = NA_real_
  )
id_cols = c("PG.ProteinGroups", "EG.ModifiedSequence", "Fragment")
mapdia_input = mapdia_input %>%
  select(id_cols, sort(setdiff(names(.), id_cols)))

fwrite(mapdia_input, 
       file=paste(data_folder, "mapDIA", "CSF_input_file.txt", sep="/"), 
       sep="\t")

# Extract model results
mapdia_model = fread(
  file=paste(data_folder, "mapDIA", "analysis_output.txt", sep="/"),sep="\t")

mapdia_model$Label = ifelse(
  mapdia_model$Protein %in% true_positives, "Positive", "Negative")

fwrite(mapdia_model, file=paste(data_folder, "mapDIA", 
                                "mapdia_model.csv", sep="/"))

print("mapDIA finished")