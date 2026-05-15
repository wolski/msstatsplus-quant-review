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

## Load + prepare inputs -------------------------------------------------------
# Load data
data_folder = ""
data_file = "20250422_140629_OKE_April22_2025_Olsen_astral_benchmark_mix_200ng_NE_report_nodecoy.tsv"
raw_input = fread(paste(data_folder, data_file, sep="/"), sep="\t")
annotation = fread(paste(data_folder, "Mix_of_Proteome_annotation.csv", sep="/"))

merged_input = merge(raw_input, annotation, by="R.FileName", 
                     all.x=TRUE, all.y=FALSE)

# MSstats+ ---------------------------------------------------------------------
# Extract order of run acquisition
annotation$Run = annotation$R.FileName
run_order = raw_input %>% 
  distinct(R.FileName, `R.Run Date (Formatted)`) %>% 
  mutate(`R.Run Date (Formatted)` = mdy_hms(`R.Run Date (Formatted)`)) %>% 
  arrange(`R.Run Date (Formatted)`)
run_order$Order = 1:nrow(run_order)
run_order$Run = run_order$R.FileName
run_order = run_order[, c("Run", "Order")]

msstats_input = MSstatsConvert::SpectronauttoMSstatsFormat(
  raw_input, annotation, intensity = 'PeakArea', 
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
fwrite(msstats_input, 
       file=paste(data_folder, "MSstats+", "MSstats+_input.csv", sep="/"))

summarized = dataProcess(msstats_input,
                         normalization=FALSE,
                         featureSubset = "topN",
                         n_top_feature = 100,
                         MBimpute=TRUE,
                         summaryMethod="linear",
                         numberOfCores = 12)
save(summarized, 
     file=paste(data_folder, "MSstats+", "MSstats+_summarized.rda", sep="/"))

msstatsplus_model = groupComparison("pairwise", summarized, 
                                 numberOfCores=12)
msstatsplus_model = msstatsplus_model$ComparisonResult
fwrite(msstatsplus_model, 
       file=paste(data_folder, "MSstats+", "MSstats+_model.csv", sep="/"))

print("MSstats+ finished")

# MSstats ----------------------------------------------------------------------
annotation$Run = annotation$R.FileName

msstats_input = MSstatsConvert::SpectronauttoMSstatsFormat(
  raw_input, annotation, intensity = 'PeakArea', 
  excludedFromQuantificationFilter = TRUE,
  filter_with_Qvalue = TRUE)
fwrite(msstats_input, 
       file=paste(data_folder, "MSstats", "msstats_input.csv", sep="/"))

msstats_input$AnomalyScores = NA
msstats_summarized = dataProcess(msstats_input, 
                                 normalization=FALSE,
                                 featureSubset = "topN",
                                 n_top_feature = 100,
                                 MBimpute=TRUE,
                                 summaryMethod="TMP",
                                 numberOfCores = 12)
save(msstats_summarized, 
     file=paste(data_folder, "MSstats", "msstats_summarized.rda", sep="/"))

msstats_model = groupComparison("pairwise", msstats_summarized, numberOfCores=12)
msstats_model = msstats_model$ComparisonResult

fwrite(msstats_model, 
       file=paste(data_folder, "MSstats", "msstats_model.csv", sep="/"))

print("MSstats finished")

# MSqRob ----------------------------------------------------------------------
input_data = merged_input %>% 
  filter(F.ExcludedFromQuantification == FALSE & 
           R.Condition != "blank" & 
           F.PeakArea > 1 & 
           PG.Qvalue < .01 &
           EG.Qvalue < .01) %>% 
  select(PG.ProteinGroups, EG.PrecursorId, F.FrgIon, 
         FG.Charge, F.Charge, R.FileName, F.PeakArea)

input_data$Fragment = paste(input_data$EG.PrecursorId, 
                            input_data$F.Charge,
                            input_data$F.FrgIon, 
                            input_data$FG.Charge, sep="_")

# Remove duplicate fragments
input_data = input_data %>% group_by(PG.ProteinGroups, R.FileName, Fragment) %>% 
  dplyr::summarize(F.PeakArea=max(F.PeakArea))
input_data = input_data %>% filter(PG.ProteinGroups != "")
input_data = setnames(input_data, 
                      c("PG.ProteinGroups", "R.FileName"),
                      c("ProteinName", "Run"))

## Take logsum of proteins
df.LFQ = as.data.frame(input_data) %>% 
  select(ProteinName, Fragment, Run, F.PeakArea)
df.LFQ = as.data.table(df.LFQ)

## Convert into MSqRob2 format
df.LFQ = dcast(df.LFQ, ProteinName + Fragment ~ Run, 
               value.var = "F.PeakArea", fun.aggregate = max, fill=NA)
colnames(df.LFQ)[3:length(colnames(df.LFQ))] = paste(
  "F.PeakArea", colnames(df.LFQ)[3:length(colnames(df.LFQ))], sep="_")

ecols = grep("F.PeakArea", colnames(df.LFQ))

pe = readQFeatures(df.LFQ, fnames = 2, quantCols = ecols,  name = "peptideRaw")

colData(pe)$condition = annotation$Condition[
  match(str_remove(colnames(df.LFQ)[grep("F.PeakArea", colnames(df.LFQ))], 
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

pe = aggregateFeatures(pe,
                       i = "peptideImp", fcol = "ProteinName",
                       name = "protein")

pe = msqrob(object = pe, i = "protein", formula = ~condition)

# Specify pairwise contrast
L = makeContrast(c("conditionE20H50Y30=0", 
                   "conditionE30H50Y20=0",
                   "conditionE40H50Y10=0",
                   "conditionE45H50Y5=0",
                   "conditionE5H50Y45=0",
                   "conditionE20H50Y30-conditionE30H50Y20=0",
                   "conditionE20H50Y30-conditionE40H50Y10=0",
                   "conditionE20H50Y30-conditionE45H50Y5=0",
                   "conditionE20H50Y30-conditionE5H50Y45=0",
                   "conditionE30H50Y20-conditionE40H50Y10=0",
                   "conditionE30H50Y20-conditionE45H50Y5=0",
                   "conditionE30H50Y20-conditionE5H50Y45=0",
                   "conditionE40H50Y10-conditionE45H50Y5=0",
                   "conditionE40H50Y10-conditionE5H50Y45=0",
                   "conditionE45H50Y5-conditionE5H50Y45=0"),
                 parameterNames = c(
                   "(Intercept)", "conditionE20H50Y30",
                   "conditionE30H50Y20", "conditionE40H50Y10", 
                   "conditionE45H50Y5", "conditionE5H50Y45"))
pe = hypothesisTest(object = pe, i = "protein", contrast = L)
save(pe, file=paste(data_folder, "msqrob2", "msqrob_obj.rda", sep="/"))

msqrob_comps = list()
comps = c("conditionE20H50Y30", "conditionE30H50Y20", "conditionE40H50Y10",
          "conditionE45H50Y5", "conditionE5H50Y45", 
          "conditionE20H50Y30 - conditionE30H50Y20",
          "conditionE20H50Y30 - conditionE40H50Y10", 
          "conditionE20H50Y30 - conditionE45H50Y5", 
          "conditionE20H50Y30 - conditionE5H50Y45", 
          "conditionE30H50Y20 - conditionE40H50Y10", 
          "conditionE30H50Y20 - conditionE45H50Y5",
          "conditionE30H50Y20 - conditionE5H50Y45", 
          "conditionE40H50Y10 - conditionE45H50Y5", 
          "conditionE40H50Y10 - conditionE5H50Y45", 
          "conditionE45H50Y5 - conditionE5H50Y45" )
for (c in comps){
  temp_comp = rowData(pe[["protein"]])[[c]]
  temp_comp$Protein = rownames(temp_comp)
  temp_comp$Label = c
  msqrob_comps[[c]] = temp_comp
}

msqrob2_model = rbindlist(msqrob_comps)
fwrite(msqrob2_model, file=paste(data_folder, "msqrob2", 
                                 "msqrob2_model.csv", sep="/"))
print("msqrob finished")

# Limma ------------------------------------------------------------------------
limma_input = merged_input %>% filter(
  F.ExcludedFromQuantification == FALSE & 
    R.Condition != "blank" & 
    F.PeakArea > 1 & 
    PG.Qvalue < .01 &
    EG.Qvalue < .01)

limma_input$Feature = paste(limma_input$EG.PrecursorId, 
                            limma_input$F.FrgIon,
                            limma_input$F.Charge, sep="_")


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
cont = makeContrasts(classE10H50Y40-classE20H50Y30, 
                     classE10H50Y40-classE30H50Y20,
                     classE10H50Y40-classE40H50Y10,
                     classE10H50Y40-classE45H50Y5,
                     classE10H50Y40-classE5H50Y45,
                     classE20H50Y30-classE30H50Y20,
                     classE20H50Y30-classE40H50Y10,
                     classE20H50Y30-classE45H50Y5,
                     classE20H50Y30-classE5H50Y45,
                     classE30H50Y20-classE40H50Y10,
                     classE30H50Y20-classE45H50Y5,
                     classE30H50Y20-classE5H50Y45,
                     classE40H50Y10-classE45H50Y5,
                     classE40H50Y10-classE5H50Y45,
                     classE45H50Y5-classE5H50Y45,
                     levels = design)

fit2 = contrasts.fit(fit1, contrasts = cont)
fit3 = eBayes(fit2)

# Loop over each contrast
all_results = lapply(1:ncol(fit3$coefficients), function(i) {
  res = topTable(fit3, coef=i, number=Inf, sort.by="none")
  
  # Calculate standard error: SE = logFC / t
  res$SE = res$logFC / res$t
  
  # Move rownames (features/proteins) into a column
  res = res %>%
    tibble::rownames_to_column(var = "Protein") %>%
    mutate(Contrast = colnames(fit3$coefficients)[i])
  
  return(res)
})

# Bind into one giant table
limma_model = bind_rows(all_results)

fwrite(limma_model, file=paste(data_folder, "limma", 
                               "limma_model.csv", sep="/"))

print("limma finished")

# Limpa ------------------------------------------------------------------------
limpa_input = merged_input %>% filter(
  F.ExcludedFromQuantification == FALSE & 
    R.Condition != "blank" & 
    F.PeakArea > 1 & 
    PG.Qvalue < .01 &
    EG.Qvalue < .01)

limpa_input$Feature = paste(limpa_input$EG.PrecursorId, 
                            limpa_input$F.FrgIon,
                            limpa_input$F.Charge, sep="_")

limpa_input$LogIntensities = log2(limpa_input$F.PeakArea)

limpa_input = limpa_input %>% 
  group_by(PG.ProteinGroups, Feature, R.FileName) %>% 
  summarize(LogIntensities=max(LogIntensities))

limpa_input = pivot_wider(limpa_input,
                          id_cols=c("PG.ProteinGroups", "Feature"),
                          names_from="R.FileName",
                          values_from="LogIntensities")

limpa_input = as.data.frame(limpa_input)
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

## dpcQuant
y.protein = dpcQuant(limpa_elist, "PG.ProteinGroups", 
                     dpc = dpcfit)

class = annotation$Condition[
  match(colnames(y.protein$E), annotation$R.FileName)] %>% 
  as.factor()
design = model.matrix(~0+class) # fitting without intercept

## dpcDE
fit = dpcDE(y.protein, design, plot=TRUE)
cont = makeContrasts(classE10H50Y40-classE20H50Y30, 
                     classE10H50Y40-classE30H50Y20,
                     classE10H50Y40-classE40H50Y10,
                     classE10H50Y40-classE45H50Y5,
                     classE10H50Y40-classE5H50Y45,
                     classE20H50Y30-classE30H50Y20,
                     classE20H50Y30-classE40H50Y10,
                     classE20H50Y30-classE45H50Y5,
                     classE20H50Y30-classE5H50Y45,
                     classE30H50Y20-classE40H50Y10,
                     classE30H50Y20-classE45H50Y5,
                     classE30H50Y20-classE5H50Y45,
                     classE40H50Y10-classE45H50Y5,
                     classE40H50Y10-classE5H50Y45,
                     classE45H50Y5-classE5H50Y45,
                     levels = design)
fit = contrasts.fit(fit, contrasts = cont)
fit = eBayes(fit)

# Loop over each contrast
all_results = lapply(1:ncol(fit$coefficients), function(i) {
  res = topTable(fit, coef=i, number=Inf, sort.by="none")
  
  # Calculate standard error: SE = logFC / t
  res$SE = res$logFC / res$t
  
  # Move rownames (features/proteins) into a column
  res = res %>%
    tibble::rownames_to_column(var = "Protein") %>%
    mutate(Contrast = colnames(fit$coefficients)[i])
  
  return(res)
})

# Bind into one giant table
limpa_model = bind_rows(all_results)

fwrite(limpa_model, file=paste(data_folder, "limpa",
                               "limpa_model.csv", sep="/"))

print("limpa finished")

# DEqMS ------------------------------------------------------------------------
deqms_input = merged_input %>% filter(
  F.ExcludedFromQuantification == FALSE & 
    R.Condition != "blank" & 
    F.PeakArea > 1 & 
    PG.Qvalue < .01 &
    EG.Qvalue < .01)

# Prep data for summarization
deqms_input$Feature = paste(deqms_input$EG.PrecursorId, 
                            deqms_input$F.FrgIon,
                            deqms_input$F.Charge, sep="_")

deqms_input = deqms_input %>% 
  select(PG.ProteinGroups,Feature, R.FileName, F.PeakArea)
deqms_input$LogIntensities = log2(deqms_input$F.PeakArea)

deqms_input = deqms_input %>% group_by(PG.ProteinGroups, Feature, R.FileName) %>% 
  summarize(LogIntensities=max(LogIntensities))
deqms_input = deqms_input %>% filter(PG.ProteinGroups != "")

deqms_input = pivot_wider(deqms_input,
                          id_cols=c("PG.ProteinGroups", "Feature"),
                          names_from="R.FileName",
                          values_from="LogIntensities")

# Adjusted DEqMS::medianSummary function without ref_col
summarize = function (dat, group_col = 2) 
{
  dat.ratio = dat
  dat.ratio[, 3:ncol(dat)] = dat.ratio[, 3:ncol(dat)] - 
    matrixStats::rowMedians(as.matrix(dat.ratio[, 3:ncol(dat)]), na.rm = TRUE)
  dat.summary = plyr::ddply(dat.ratio, colnames(dat)[group_col], 
                            function(x) matrixStats::colMedians(
                              as.matrix(x[, 3:ncol(dat)]), na.rm = TRUE))
  colnames(dat.summary)[2:ncol(dat.summary)] = colnames(dat)[3:ncol(dat)]
  dat.new = dat.summary[, -1]
  rownames(dat.new) = dat.summary[, 1]
  return(dat.new)
}

deqms_summarized = summarize(deqms_input, group_col=1)
pep_count = deqms_input %>% group_by(PG.ProteinGroups) %>%
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
cont = makeContrasts(classE10H50Y40-classE20H50Y30, 
                     classE10H50Y40-classE30H50Y20,
                     classE10H50Y40-classE40H50Y10,
                     classE10H50Y40-classE45H50Y5,
                     classE10H50Y40-classE5H50Y45,
                     classE20H50Y30-classE30H50Y20,
                     classE20H50Y30-classE40H50Y10,
                     classE20H50Y30-classE45H50Y5,
                     classE20H50Y30-classE5H50Y45,
                     classE30H50Y20-classE40H50Y10,
                     classE30H50Y20-classE45H50Y5,
                     classE30H50Y20-classE5H50Y45,
                     classE40H50Y10-classE45H50Y5,
                     classE40H50Y10-classE5H50Y45,
                     classE45H50Y5-classE5H50Y45,
                     levels = design)
fit2 = contrasts.fit(fit1, contrasts = cont)
fit3 = eBayes(fit2)

fit3$count = pep_count[rownames(fit3$coefficients), "count"]

fit4 = spectraCounteBayes(fit3)

## Analyze results
# Loop over each contrast
all_results = lapply(1:ncol(fit4$coefficients), function(i) {
  res = topTable(fit4, coef=i, number=Inf, sort.by="none")
  
  # Calculate standard error: SE = logFC / t
  res$SE = res$logFC / res$t
  
  # Move rownames (features/proteins) into a column
  res = res %>%
    tibble::rownames_to_column(var = "Protein") %>%
    mutate(Contrast = colnames(fit4$coefficients)[i])
  
  return(res)
})

# Bind into one giant table
deqms_model = bind_rows(all_results)

fwrite(deqms_model, file=paste(data_folder, "DEqMS", 
                               "deqms_model.csv", sep="/"))

print("DeqMS finished")

# mapDIA -----------------------------------------------------------------------
# Run externally (outside of R). Here we just prepare the input.
# Use same data processesing as for limma
mapdia_input = merged_input %>% filter(
  F.ExcludedFromQuantification == FALSE & 
    R.Condition != "blank" & 
    F.PeakArea > 1 & 
    PG.Qvalue < .01 &
    EG.Qvalue < .01)

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

mapdia_input = mapdia_input %>% merge(rt_data, by=c("EG.ModifiedSequence", "Fragment"))
mapdia_input = mapdia_input %>%
  select(id_cols, sort(setdiff(names(.), id_cols)))

fwrite(mapdia_input, 
       file=paste(data_folder, "mapDIA", "Mix_input_file.txt", sep="/"), 
       sep="\t")

# Extract model results
mapdia_model = fread(
  file=paste(data_folder, "mapDIA", "analysis_output.txt", sep="/"),sep="\t")

fwrite(mapdia_model, file=paste(data_folder, "mapDIA", 
                                "mapdia_model.csv", sep="/"))

print("mapDIA finished")

# Analysis ---------------------------------------------------------------------
# Load data
msstats_input = fread(file=paste(data_folder, "MSstats+", "MSstats+_input.csv", sep="/"))

load(file=paste(data_folder, "MSstats+", "MSstats+_summarized.rda", sep="/"))

protein_mapping = fread(file=paste(data_folder, "idmapping.tsv", sep="/"),
                        sep="\t")

weighted_model = fread(file=paste(data_folder, 
                                  "MSstats+_model.csv", sep="/"))
weighted_model = merge(weighted_model,
                       protein_mapping %>% select(Entry, Organism),
                       all.x=TRUE, all.y=FALSE,
                       by.x="Protein", by.y="Entry")

msstats_model = fread(file=paste(data_folder, "MSstats", 
                                 "msstats_model.csv", sep="/"))
msstats_model = merge(msstats_model,
                      protein_mapping %>% select(Entry, Organism),
                      all.x=TRUE, all.y=FALSE, 
                      by.x="Protein", by.y="Entry")

msqrob2_model = fread(file=paste(data_folder, "msqrob2", 
                          "msqrob2_model.csv", sep="/"))
msqrob2_model = merge(msqrob2_model,
                      protein_mapping %>% select(Entry, Organism),
                      all.x=TRUE, all.y=FALSE, 
                      by.x="Protein", by.y="Entry")

limma_model = fread(file=paste(data_folder, "limma", 
                        "limma_model.csv", sep="/"))
limma_model = merge(limma_model,
                    protein_mapping %>% select(Entry, Organism),
                    all.x=TRUE, all.y=FALSE, 
                    by.x="Protein", by.y="Entry")

limpa_model = fread(file=paste(data_folder, "limpa",
                        "limpa_model.csv", sep="/"))
limpa_model = merge(limpa_model,
                    protein_mapping %>% select(Entry, Organism),
                    all.x=TRUE, all.y=FALSE, 
                    by.x="Protein", by.y="Entry")

deqms_model = fread(file=paste(data_folder, "DEqMS", 
                               "deqms_model.csv", sep="/"))
deqms_model = merge(deqms_model,
                    protein_mapping %>% select(Entry, Organism),
                    all.x=TRUE, all.y=FALSE, 
                    by.x="Protein", by.y="Entry")

mapdia_model = fread(file=paste(data_folder, "mapDIA", 
                               "mapdia_model.csv", sep="/"))
mapdia_model = merge(mapdia_model,
                    protein_mapping %>% select(Entry, Organism),
                    all.x=TRUE, all.y=FALSE, 
                    by.x="Protein", by.y="Entry")

weighted_model %>% filter(!is.na(Organism)) %>% group_by(Organism) %>% 
  summarize(n_tested = n_distinct(Protein))

msstats_model %>% filter(!is.na(Organism)) %>% group_by(Organism) %>% 
  summarize(n_tested = n_distinct(Protein))

msqrob2_model %>% filter(!is.na(Organism)) %>% group_by(Organism) %>% 
  summarize(n_tested = n_distinct(Protein))

limma_model %>% filter(!is.na(Organism)) %>% group_by(Organism) %>% 
  summarize(n_tested = n_distinct(Protein))

limpa_model %>% filter(!is.na(Organism)) %>% group_by(Organism) %>% 
  summarize(n_tested = n_distinct(Protein))

deqms_model %>% filter(!is.na(Organism)) %>% group_by(Organism) %>% 
  summarize(n_tested = n_distinct(Protein))

mapdia_model %>% filter(!is.na(Organism)) %>% group_by(Organism) %>% 
  summarize(n_tested = n_distinct(Protein))


# Summary plots

health_info = MSstatsConvert::CheckDataHealth(msstats_input)
skew_score = health_info[[2]]
ggplot(skew_score, aes(x = skew)) + 
  geom_histogram(fill = "#009E73", color = "black", binwidth = 0.2) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1.5) + 
  theme_minimal(base_size = 16) + 
  labs(
    x = "Pearson's moment coefficient of skewness",
    y = "Count"
  ) +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18)
  ) +
  xlim(-1, 7) +
  scale_y_continuous(expand = c(0, 0))

msstats_input %>%
  ggplot(aes(x = log2(Intensity), y = AnomalyScores)) +
  geom_hex(bins = 50) +
  scale_fill_gradientn(
    colours = c("grey90", "#fcae91", "#fb6a4a", "#cb181d"), 
    values = scales::rescale(c(0, 0.3, 0.6, 1)),  # More contrast in low-mid range
    trans = "sqrt",
    guide = "none"
  ) +
  theme_minimal(base_size = 16) +
  labs(
    x = expression(log[2] * " Intensity"),
    y = "Anomaly Score",
    title = ""
  ) +
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 26),
    panel.grid.major = element_line(size = 0.3),
    panel.grid.minor = element_blank()
  )


# Calculate summary statistics
calculate_fdr = function(dt, pval_col, label_col){
  dt = dt[!is.na(dt$Organism),]
  dt[, is_human := grepl("Homo sapiens", Organism)]
  dt[, is_significant := get(pval_col) < 0.05]  
  
  fdr_by_label = dt[, {
    n_total = .N
    n_discovered = sum(is_significant, na.rm=TRUE)
    n_false_positive = sum(is_human & is_significant, na.rm=TRUE)
    fdr = n_false_positive / max(n_discovered, 1)
    .(n_total, n_discovered, n_false_positive, fdr)
  }, by = get(label_col)]
  
  return(fdr_by_label)
}

calculate_tpr = function(dt, pval_col, label_col){
  dt = dt[!is.na(dt$Organism),]
  dt[, p := !grepl("Homo sapiens", Organism)]
  dt[, tp := (get(pval_col) < 0.05) & (!grepl("Homo sapiens", Organism))]
  
  tpr_by_label = dt[, {
    n_positive = sum(p, na.rm=TRUE)
    n_tp = sum(tp, na.rm=TRUE)
    tpr = n_tp/n_positive
    .(n_positive, n_tp, tpr)
  }, by = get(label_col)]
  
  return(tpr_by_label)
}

prop_fdr = calculate_fdr(as.data.table(weighted_model), "adj.pvalue", "Label")
msstats_fdr = calculate_fdr(msstats_model, "adj.pvalue", "Label")
msqrob2_fdr = calculate_fdr(msqrob2_model, "adjPval", "Label")
limma_fdr = calculate_fdr(limma_model, "adj.P.Val", "Contrast")
limpa_fdr = calculate_fdr(limpa_model, "adj.P.Val", "Contrast")
deqms_fdr = calculate_fdr(as.data.table(deqms_model), "adj.P.Val", "Contrast")
mapdia_fdr = calculate_fdr(as.data.table(mapdia_model), "FDR", "Label2")

1 - prop_fdr %>% dplyr::summarize(sum(n_false_positive) / sum(n_discovered))
1 - msstats_fdr %>% dplyr::summarize(sum(n_false_positive) / sum(n_discovered))
1 - msqrob2_fdr %>% dplyr::summarize(sum(n_false_positive) / sum(n_discovered))
1 - limma_fdr %>% dplyr::summarize(sum(n_false_positive) / sum(n_discovered))
1 - limpa_fdr %>% dplyr::summarize(sum(n_false_positive) / sum(n_discovered))
1 - deqms_fdr %>% dplyr::summarize(sum(n_false_positive) / sum(n_discovered))
1 - mapdia_fdr %>% dplyr::summarize(sum(n_false_positive) / sum(n_discovered))


prop_tpr = calculate_tpr(as.data.table(weighted_model), "adj.pvalue", "Label")
msstats_tpr = calculate_tpr(msstats_model, "adj.pvalue", "Label")
msqrob2_tpr = calculate_tpr(msqrob2_model, "adjPval", "Label")
limma_tpr = calculate_tpr(limma_model, "adj.P.Val", "Contrast")
limpa_tpr = calculate_tpr(limpa_model, "adj.P.Val", "Contrast")
deqms_tpr = calculate_tpr(as.data.table(deqms_model), "adj.P.Val", "Contrast")
mapdia_tpr = calculate_tpr(as.data.table(mapdia_model), "FDR", "Label2")

prop_tpr %>% dplyr::summarize(sum(n_tp) / sum(n_positive))
msstats_tpr %>% dplyr::summarize(sum(n_tp) / sum(n_positive))
msqrob2_tpr %>% dplyr::summarize(sum(n_tp) / sum(n_positive))
limma_tpr %>% dplyr::summarize(sum(n_tp) / sum(n_positive))
limpa_tpr %>% dplyr::summarize(sum(n_tp) / sum(n_positive))
deqms_tpr %>% dplyr::summarize(sum(n_tp) / sum(n_positive))
mapdia_tpr %>% dplyr::summarize(sum(n_tp) / sum(n_positive))

# Plot example
protein = "P56134"#P00338"#P62861,P30405, O95674, O95084, O14618

weighted_model %>% 
  filter(Protein == protein)
msstats_model %>% 
  filter(Protein == protein)
limma_model %>% 
  filter(Protein == protein)
limpa_model %>% 
  filter(Protein == protein)
msqrob2_model %>% 
  filter(Protein == protein)
deqms_model %>% 
  filter(Protein == protein)
mapdia_model %>% 
  filter(Protein == protein)

ggplot() + 
  geom_boxplot(data = weighted_model %>% filter(Protein == protein),
               aes(x="MSstats+", y=SE)) + 
  geom_boxplot(data = msstats_model %>% filter(Protein == protein),
               aes(x="MSstats", y=SE)) + 
  geom_boxplot(data = limma_model %>% filter(Protein == protein),
               aes(x="limma", y=SE)) + 
  geom_boxplot(data = limpa_model %>% filter(Protein == protein),
               aes(x="limpa", y=SE)) + 
  geom_boxplot(data = msqrob2_model %>% filter(Protein == protein),
               aes(x="msqrob", y=se)) + 
  geom_boxplot(data = deqms_model %>% filter(Protein == protein),
               aes(x="DEqMS", y=SE)) + theme_bw()

unique_correct

# Create ordering based on Condition
peps = summarized$FeatureLevelData %>%
  filter(PROTEIN == protein) %>% distinct(PEPTIDE) %>% unlist()
charge = summarized$FeatureLevelData %>%
  filter(PROTEIN == protein) %>% distinct(PEPTIDE) %>% unlist() %>% 
  str_split_i("_",i=2)

plot_data = msstats_input %>% 
  mutate(PEPTIDE =paste(PeptideSequence, PrecursorCharge, sep="_"), 
         FEATURE = paste(PeptideSequence, PrecursorCharge, FragmentIon, 
                         ProductCharge,
                         sep="_")) %>% 
  filter(PEPTIDE %in% peps) %>%
  arrange(Condition, Run) %>%
  mutate(Run = factor(Run, levels = unique(Run)))

plot_data = as.data.table(plot_data)
plot_data[, Condition := gsub(
  "([A-Z])(\\d+)", "\\1\\2_", Condition)]
plot_data[, Condition := sub("_$", "", Condition)]


# Get midpoint for condition labels
# Total number of runs
n_runs = length(levels(plot_data$Run))

# Compute x-intercepts at 3.5, 6.5, 9.5, ...
vline_positions = seq(3.5, n_runs, by = 3)

# Midpoints for condition labels
condition_labels = plot_data %>%
  group_by(Condition) %>%
  summarize(mid = mean(as.numeric(Run)), .groups = "drop")

  # Plot
ggplot(plot_data%>% filter(PEPTIDE == "ASVGEC[Carbamidomethyl (C)]PAPVPVK_2")) +# 
    geom_line(aes(x = Run, y = AnomalyScores, group = PEPTIDE),
              color="darkorchid4", linewidth = 1.25) +
    geom_point(aes(x = Run, y = AnomalyScores),
               color="darkorchid4", size = 3.5) +
    geom_vline(xintercept = vline_positions, 
               linetype = "dotted", color = "black") +
    geom_text(data = condition_labels, aes(x = mid, y = Inf, label = Condition),
              vjust = 1.5, size = 6, fontface = "bold") +
    theme_bw(base_size=22) +
    ylim(.2, .65) +
    theme(
      legend.position = "none",
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16), 
      axis.text.x = element_blank(),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 24),
      plot.title = element_text(size = 32)
    ) + 
    labs(title = "", x = "", y = "") #+ ylim(.15, .7)

ggplot(plot_data %>% filter(PEPTIDE == "ASVGEC[Carbamidomethyl (C)]PAPVPVK_2")) +
  geom_line(aes(x = Run, `FGShapeQualityScore(MS2)`, group = PrecursorCharge),
            color = "darkred", linewidth = 1.25) +
  geom_point(aes(x = Run, y = `FGShapeQualityScore(MS2)`), 
             color = "darkred", size = 3.5) +
  geom_vline(xintercept = vline_positions, 
             linetype = "dotted", color = "black") +
  geom_text(data = condition_labels, aes(x = mid, y = Inf, label = Condition),
            vjust = 1.5, size = 6, fontface = "bold", angle=0) +
  theme_bw(base_size=22) +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16), 
    axis.text.x = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 24),
    plot.title = element_text(size = 32)
  ) + 
  labs(title = "", x = "", y = "") + ylim(.2, .75)

ggplot(plot_data %>% filter(PEPTIDE == "ASVGEC[Carbamidomethyl (C)]PAPVPVK_2")) +
  geom_line(aes(x = Run, `EGDeltaRT`, group = PrecursorCharge),
            color = "darkred", linewidth = 1.25) +
  geom_point(aes(x = Run, y = `EGDeltaRT`), 
             color = "darkred", size = 3.5) +
  geom_vline(xintercept = vline_positions, 
             linetype = "dotted", color = "black") +
  geom_text(data = condition_labels, aes(x = mid, y = Inf, label = Condition),
            vjust = 1.5, size = 6, fontface = "bold", angle=0) +
  theme_bw(base_size=22) +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16), 
    axis.text.x = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 24),
    plot.title = element_text(size = 32)
  ) + 
  labs(title = "", x = "", y = "") + ylim(-.05, .35)

plot = plot_data %>% filter(PEPTIDE == "ASVGEC[Carbamidomethyl (C)]PAPVPVK_2") %>% 
  merge(summarized$FeatureLevelData, all.x = TRUE, all.y = FALSE, 
        by.x = c("FEATURE", "Run"), by.y = c("FEATURE", "originalRUN")) %>% 
  ggplot() + 
  geom_line(aes(x = Run, y = newABUNDANCE, group = FragmentIon, 
                color = FragmentIon), linewidth = 1.25) +
  geom_point(aes(x = Run, y = newABUNDANCE, group = FragmentIon, 
                 color = FragmentIon), size = 3.5) + 
  geom_point(aes(x = Run, y = predicted, group = FragmentIon), 
             color = "white", size = 2) +
  geom_vline(xintercept = vline_positions, 
             linetype = "dotted", color = "black") +
  geom_text(data = condition_labels, aes(x = mid, y = Inf, label = Condition),
            vjust = 1.5, size = 6, fontface = "bold") +
  scale_color_viridis_d(option = "D") + 
  theme_bw(base_size=22) + 
  theme(
    # legend.position = "none",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16), 
    axis.text.x = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 24),
    plot.title = element_text(size = 32)
  ) + 
  labs(title = "", x = "", y = "") #+ ylim(4.75, 9.25)
print(plot)

plot_data = summarized$ProteinLevelData %>%
  filter(Protein == protein) %>%
  arrange(GROUP, RUN) %>%
  mutate(RUN = factor(RUN, levels = unique(RUN)))

# Vertical lines at 3.5, 6.5, 9.5, ...
vlines = seq(3.5, length(unique(plot_data$RUN)), by = 3)

# Midpoints for condition labels
group_labels = plot_data %>%
  group_by(GROUP) %>%
  summarize(mid = mean(as.numeric(RUN)), .groups = "drop")


plot_data %>% ggplot() + geom_boxplot(aes(x=GROUP, y=Variance))

plot =
  ggplot() + 
  geom_line(data = summarized$FeatureLevelData %>% filter(PROTEIN == protein),
            aes(x = RUN, y = newABUNDANCE, group = FEATURE),
            color = "grey", linewidth = 1.25) +
  geom_point(data = summarized$FeatureLevelData %>% filter(PROTEIN == protein), 
             aes(x = RUN, y = newABUNDANCE, group = FEATURE), 
             size = 3.5, color = "grey") + 
  geom_line(data = plot_data,
            aes(x = RUN, y = LogIntensities, group = Protein),
            color = "darkred", linewidth = 1.25) +
  geom_point(data = plot_data, 
             aes(x = RUN, y = LogIntensities), 
             size = 3.5, color = "darkred") + 
  geom_errorbar(data=plot_data,
                aes(x=RUN, ymin = LogIntensities-sqrt(Variance),
                    ymax = LogIntensities+sqrt(Variance)),
                width = .75, linewidth=1.1, color="darkred") + 
  geom_vline(xintercept = vline_positions, 
             linetype = "dotted", color = "black") +
  geom_text(data = group_labels, aes(x = mid, y = Inf, label = GROUP),
            vjust = 1.5, size = 6, fontface = "bold") +
  scale_color_viridis_d(option = "D") + 
  theme_bw(base_size=22) +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16), 
    axis.text.x = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 24),
    plot.title = element_text(size = 32)
  ) + 
  labs(title = "", x = "", y ="") + ylim(0, 16)
print(plot)

# distribution of fold change vizualization
human_ratios = rep(0,15)
ecoli_ratios = log2(c(10/20,10/30,10/40,10/45,10/5,20/30,20/40,
                      20/45,20/5,30/40,30/45,30/5,40/45,40/5,45/5))
yeast_ratios = log2(c(40/30,40/20,40/10,40/5,40/45,30/20,30/10,
                      30/5,30/45,20/10,20/5,20/45,10/5,10/45,5/45))


comparisons = c(
  "E10H50Y40 vs E20H50Y30",
  "E10H50Y40 vs E30H50Y20",
  "E10H50Y40 vs E40H50Y10",
  "E10H50Y40 vs E45H50Y5",
  "E10H50Y40 vs E5H50Y45",
  "E20H50Y30 vs E30H50Y20",
  "E20H50Y30 vs E40H50Y10",
  "E20H50Y30 vs E45H50Y5",
  "E20H50Y30 vs E5H50Y45",
  "E30H50Y20 vs E40H50Y10",
  "E30H50Y20 vs E45H50Y5",
  "E30H50Y20 vs E5H50Y45",
  "E40H50Y10 vs E45H50Y5",
  "E40H50Y10 vs E5H50Y45",
  "E45H50Y5 vs E5H50Y45"
)

limma_comparisons = c(
  "classE10H50Y40 - classE20H50Y30",
  "classE10H50Y40 - classE30H50Y20",
  "classE10H50Y40 - classE40H50Y10",
  "classE10H50Y40 - classE45H50Y5",
  "classE10H50Y40 - classE5H50Y45",
  "classE20H50Y30 - classE30H50Y20",
  "classE20H50Y30 - classE40H50Y10",
  "classE20H50Y30 - classE45H50Y5",
  "classE20H50Y30 - classE5H50Y45",
  "classE30H50Y20 - classE40H50Y10",
  "classE30H50Y20 - classE45H50Y5",
  "classE30H50Y20 - classE5H50Y45",
  "classE40H50Y10 - classE45H50Y5",
  "classE40H50Y10 - classE5H50Y45",
  "classE45H50Y5 - classE5H50Y45"
)

msqrob_contrasts = c(
  "conditionE20H50Y30",
  "conditionE30H50Y20",
  "conditionE40H50Y10",
  "conditionE45H50Y5",
  "conditionE5H50Y45",
  "conditionE20H50Y30 - conditionE30H50Y20",
  "conditionE20H50Y30 - conditionE40H50Y10",
  "conditionE20H50Y30 - conditionE45H50Y5",
  "conditionE20H50Y30 - conditionE5H50Y45",
  "conditionE30H50Y20 - conditionE40H50Y10",
  "conditionE30H50Y20 - conditionE45H50Y5",
  "conditionE30H50Y20 - conditionE5H50Y45",
  "conditionE40H50Y10 - conditionE45H50Y5",
  "conditionE40H50Y10 - conditionE5H50Y45",
  "conditionE45H50Y5 - conditionE5H50Y45"
)

# Combine into a single data.frame
comparison_map = data.frame(
  Label = comparisons,
  limma = limma_comparisons,
  msqrob = msqrob_contrasts,
  human_ratios = human_ratios,
  ecoli_ratios = ecoli_ratios,
  yeast_ratios = yeast_ratios,
  stringsAsFactors = FALSE
)


# Merge step-by-step on matching keys
weighted_model = merge(weighted_model, comparison_map, by = "Label", all.x = TRUE)
weighted_model = weighted_model %>% select(Label, Organism, Protein, log2FC, SE, 
                                           adj.pvalue, human_ratios, 
                                           ecoli_ratios, yeast_ratios)
weighted_model$Model = "MSstats+"

msstats_model = merge(msstats_model, comparison_map, by = "Label", all.x = TRUE)
msstats_model = msstats_model %>% select(Label, Organism, Protein, log2FC, SE,
                                         adj.pvalue, human_ratios, 
                                         ecoli_ratios, yeast_ratios)
msstats_model$Model = "MSstats"

msqrob2_model = merge(msqrob2_model, comparison_map, by.x = "Label", by.y="msqrob", all.x = TRUE)
msqrob2_model = msqrob2_model %>% select(Label.y, Organism, Protein, logFC, se,
                                         adjPval, human_ratios, 
                                         ecoli_ratios, yeast_ratios)
setnames(msqrob2_model, 
         c("Label.y", "logFC", "se", "adjPval"), 
         c("Label", "log2FC", "SE", "adj.pvalue"))
msqrob2_model$Model = "msqrob2"
msqrob2_model$log2FC = ifelse(msqrob2_model$Label %in% c(  "E10H50Y40 vs E20H50Y30",
                                                           "E10H50Y40 vs E30H50Y20",
                                                           "E10H50Y40 vs E40H50Y10",
                                                           "E10H50Y40 vs E45H50Y5",
                                                           "E10H50Y40 vs E5H50Y45"),
                              -msqrob2_model$log2FC,
                              msqrob2_model$log2FC)
                              

limma_model = merge(limma_model, comparison_map, by.x="Contrast", by.y = "limma", all.x = TRUE)
limma_model = limma_model %>% select(Label, Organism, Protein, logFC, SE, adj.P.Val,
                                     human_ratios, ecoli_ratios, yeast_ratios)
setnames(limma_model, c("logFC", "adj.P.Val"), c("log2FC", "adj.pvalue"))
limma_model$Model = "limma"

limpa_model = merge(limpa_model, comparison_map, by.x="Contrast", by.y = "limma", all.x = TRUE)
limpa_model = limpa_model %>% select(Label, Organism, Protein, logFC, SE,adj.P.Val,
                                     human_ratios, ecoli_ratios, yeast_ratios)
setnames(limpa_model, c("logFC", "adj.P.Val"), c("log2FC", "adj.pvalue"))
limpa_model$Model = "limpa"

deqms_model = merge(as.data.table(deqms_model), comparison_map, by.x="Contrast", by.y = "limma", all.x = TRUE)
deqms_model = deqms_model %>% select(Label, Organism, Protein, logFC, SE, adj.P.Val,
                                     human_ratios, ecoli_ratios, yeast_ratios)
setnames(deqms_model, c("logFC", "adj.P.Val"), c("log2FC", "adj.pvalue"))
deqms_model$Model = "DEqMS"

comparison_df = rbindlist(list(
  weighted_model, msstats_model, msqrob2_model, limma_model, limpa_model, deqms_model
))

comparison_df$error = ifelse(grepl("Homo", comparison_df$Organism), 
                             comparison_df$log2FC - comparison_df$human_ratios,
                             ifelse(grepl("yeast", comparison_df$Organism),
                                    comparison_df$log2FC - comparison_df$yeast_ratios,
                                    comparison_df$log2FC - comparison_df$ecoli_ratios))

# Your palette
color_p <- c("#E69F00", "#56B4E9", "#009E73",
             "#F0E442", "#0072B2", "#D55E00",
             "#CC79A7")

pal_named <- c(
  "MSstats+"         = color_p[1],
  "MSstats"          = color_p[2],
  "msqrob2"          = color_p[3],
  "limma"            = color_p[4],
  "DEqMS"            = color_p[5],
  "limpa"            = color_p[6]
  # color_p[7] reserved if you add another method
)

comparison_df %>%
  filter(!is.na(Organism)) %>%
  ggplot(aes(x = Organism, y = abs(error), fill = Model)) +
  geom_boxplot(width = 0.75,   outliers = FALSE) +
  scale_fill_manual(values = pal_named, name = "Model") +
  ylim(0, 2) +
  theme_minimal(base_size = 30) +
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

comparison_df %>%
  filter(!is.na(Organism)) %>%
  ggplot(aes(x = Organism, y = SE, fill = Model)) +
  geom_boxplot(width = 0.75,   outliers = FALSE) +
  scale_fill_manual(values = pal_named, name = "Model") +
  ylim(0, 1) +
  theme_minimal(base_size = 30) +
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

# AUROC curves -----------------------------------------------------------------
calc_tpr_fpr_curve <- function(dt,
                               alpha_seq = seq(0, 1.0, by = 0.0001)) {
  stopifnot(is.data.table(dt))
  
  dt = copy(dt)
  
  pos_org = c(
    "Escherichia coli (strain K12)",
    "Saccharomyces cerevisiae (strain ATCC 204508 / S288c) (Baker's yeast)"
  )
  neg_org = "Homo sapiens (Human)"
  
  # truth: TRUE = should be DA (Ecoli/Yeast), FALSE = null (Human)
  dt[, truth := fifelse(
    Organism %chin% pos_org, TRUE,
    fifelse(Organism == neg_org, FALSE, NA)
  )]
  
  # keep only defined truth + non-missing adj.pvalue
  dt = dt[!is.na(truth) & !is.na(adj.pvalue)]
  
  # For strict "< alpha" with findInterval (which effectively counts <= x),
  # subtract a tiny epsilon from alpha.
  eps = .Machine$double.eps
  
  out = dt[, {
    p_pos = sort(adj.pvalue[truth])
    p_neg = sort(adj.pvalue[!truth])
    
    P = length(p_pos)
    N = length(p_neg)
    
    TP = findInterval(alpha_seq - eps, p_pos)
    FP = findInterval(alpha_seq - eps, p_neg)
    
    TPR = if (P > 0) TP / P else rep(NA_real_, length(alpha_seq))
    FPR = if (N > 0) FP / N else rep(NA_real_, length(alpha_seq))
    FDP = FP / pmax(TP + FP, 1)
    
    data.table(
      alpha = alpha_seq,
      TP = TP,
      FP = FP,
      FN = P - TP,
      TN = N - FP,
      TPR = TPR,
      FPR = FPR,
      Pos = P,
      Neg = N,
      FDP=FDP
    )
  }, by = .(Model, Label)]
  
  setorder(out, Model, Label, alpha)
  out
}


curve_tbl = calc_tpr_fpr_curve(comparison_df %>% filter(is.finite(log2FC)))

curve_tbl %>%
  ggplot(aes(FDP, TPR, color = Model)) +
  
  ## ROC curves
  geom_line(linewidth = 1.2) +
  
  ## Random classifier reference
  geom_abline(
    slope = 1, intercept = 0,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.8
  ) +
  
  ## Equal axes (CRITICAL for ROC)
  coord_equal(xlim = c(0,1), ylim = c(0,1), expand = TRUE) +
  
  facet_wrap(~Label, ncol = 5) +
  
  scale_color_manual(values = pal_named) +
  
  labs(
    x = "False Discovery Proportion (FDP)",
    y = "True Positive Rate (TPR)",
    color = NULL
  ) +
  
  theme_bw(base_size = 20) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold"),
    
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    aspect.ratio = .8,
    legend.text = element_text(size = 22)
  )
calc_auroc <- function(curve_tbl) {
  
  stopifnot(is.data.table(curve_tbl))
  
  curve_tbl[
    ,
    {
      dt = copy(.SD)
      
      dt = dt[!is.na(FPR) & !is.na(TPR)]
      dt = dt[order(FPR, TPR)]
      
      # keep best TPR for duplicated FPR values
      dt = dt[, .(TPR = max(TPR)), by = FPR]
      
      # add endpoints if missing
      if (nrow(dt) == 0 || dt$FPR[1] > 0) {
        dt = rbind(
          data.table(FPR = 0, TPR = 0),
          dt
        )
      }
      if (dt$FPR[nrow(dt)] < 1) {
        dt = rbind(
          dt,
          data.table(FPR = 1, TPR = 1)
        )
      }
      
      dt = dt[order(FPR)]
      
      auc = sum(
        diff(dt$FPR) *
          (head(dt$TPR, -1) + tail(dt$TPR, -1)) / 2
      )
      
      .(AUROC = auc)
    },
    by = .(Model, Label)
  ]
}

auroc_tbl = calc_auroc(curve_tbl)

auroc_summary = auroc_tbl[
  ,
  .(
    mean_AUROC = mean(AUROC, na.rm = TRUE),
    sd_AUROC   = sd(AUROC, na.rm = TRUE)
  ),
  by = Model
][order(-mean_AUROC)]
