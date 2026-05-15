# Mix of Proteome v1_log2 rerun
#
# Variant of Mixture_of_proteomes_processing.R that re-runs all 6 methods
# (MSstats+, MSstats, msqrob2, MaxLFQ+limma, limpa, DEqMS) with log2 as the
# per-feature intensity transform (the existing default) and writes outputs
# to V1_log2/. Reproduces the Mixture-of-species column of Table 1.
#
# Code, filters, contrasts and thresholds are taken verbatim from
# Mixture_of_proteomes_processing.R; only the output paths are redirected
# and a clean comparison-table builder replaces the inline analysis block
# at the bottom of the original.

## Packages --------------------------------------------------------------------
library(data.table)
library(tidyverse)
library(MSstats)
library(MSstatsConvert)
library(QFeatures)
library(msqrob2)
library(iq)
library(limma)
library(limpa)
library(DEqMS)

## Setup -----------------------------------------------------------------------
data_folder = "."
output_dir  = file.path(data_folder, "V1_log2")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

data_file = "20250422_140629_OKE_April22_2025_Olsen_astral_benchmark_mix_200ng_NE_report_nodecoy.tsv"
raw_input = fread(paste(data_folder, data_file, sep="/"), sep="\t")
annotation = fread(paste(data_folder, "Mix_of_Proteome_annotation.csv", sep="/"))

merged_input = merge(raw_input, annotation, by="R.FileName",
                     all.x=TRUE, all.y=FALSE)

# MSstats+ ---------------------------------------------------------------------
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
  numberOfCores=4)
fwrite(msstats_input,
       file=file.path(output_dir, "MSstatsplus_input.csv"))

summarized = dataProcess(msstats_input,
                         normalization=FALSE,
                         featureSubset = "topN",
                         n_top_feature = 100,
                         MBimpute=TRUE,
                         summaryMethod="linear",
                         numberOfCores=4)
save(summarized,
     file=file.path(output_dir, "MSstatsplus_summarized.rda"))

msstatsplus_model = groupComparison("pairwise", summarized,
                                    numberOfCores=4)
msstatsplus_model = msstatsplus_model$ComparisonResult
fwrite(msstatsplus_model,
       file=file.path(output_dir, "MSstatsplus_model.csv"))

print("MSstats+ finished")

# MSstats ----------------------------------------------------------------------
annotation$Run = annotation$R.FileName

msstats_input2 = MSstatsConvert::SpectronauttoMSstatsFormat(
  raw_input, annotation, intensity = 'PeakArea',
  excludedFromQuantificationFilter = TRUE,
  filter_with_Qvalue = TRUE)
fwrite(msstats_input2,
       file=file.path(output_dir, "msstats_input.csv"))

msstats_input2$AnomalyScores = NA
msstats_summarized = dataProcess(msstats_input2,
                                 normalization=FALSE,
                                 featureSubset = "topN",
                                 n_top_feature = 100,
                                 MBimpute=TRUE,
                                 summaryMethod="TMP",
                                 numberOfCores=4)
save(msstats_summarized,
     file=file.path(output_dir, "msstats_summarized.rda"))

msstats_model = groupComparison("pairwise", msstats_summarized, numberOfCores=4)
msstats_model = msstats_model$ComparisonResult

fwrite(msstats_model,
       file=file.path(output_dir, "msstats_model.csv"))

print("MSstats finished")

# MSqRob2 ---------------------------------------------------------------------
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

input_data = input_data %>% group_by(PG.ProteinGroups, R.FileName, Fragment) %>%
  dplyr::summarize(F.PeakArea=max(F.PeakArea))
input_data = input_data %>% filter(PG.ProteinGroups != "")
input_data = setnames(input_data,
                      c("PG.ProteinGroups", "R.FileName"),
                      c("ProteinName", "Run"))

df.LFQ = as.data.frame(input_data) %>%
  select(ProteinName, Fragment, Run, F.PeakArea)
df.LFQ = as.data.table(df.LFQ)

df.LFQ = dcast(df.LFQ, ProteinName + Fragment ~ Run,
               value.var = "F.PeakArea", fun.aggregate = max, fill=NA)
colnames(df.LFQ)[3:length(colnames(df.LFQ))] = paste(
  "F.PeakArea", colnames(df.LFQ)[3:length(colnames(df.LFQ))], sep="_")

ecols = grep("F.PeakArea", colnames(df.LFQ))

pe = readQFeatures(df.LFQ, fnames = 2, quantCols = ecols, name = "peptideRaw")

colData(pe)$condition = annotation$Condition[
  match(str_remove(colnames(df.LFQ)[grep("F.PeakArea", colnames(df.LFQ))],
                   "F.PeakArea_"), annotation$R.FileName)] %>%
  as.factor()

rowData(pe[["peptideRaw"]])$nNonZero = rowSums(assay(pe[["peptideRaw"]]) > 0,
                                               na.rm = TRUE)
pe = zeroIsNA(pe, "peptideRaw")

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
save(pe, file=file.path(output_dir, "msqrob_obj.rda"))

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
          "conditionE45H50Y5 - conditionE5H50Y45")
for (c in comps){
  temp_comp = rowData(pe[["protein"]])[[c]]
  temp_comp$Protein = rownames(temp_comp)
  temp_comp$Label = c
  msqrob_comps[[c]] = temp_comp
}

msqrob2_model = rbindlist(msqrob_comps)
fwrite(msqrob2_model, file=file.path(output_dir, "msqrob2_model.csv"))
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
  pdf_out = file.path(output_dir, "qc-plots.pdf"),
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
design = model.matrix(~0+class)

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

all_results = lapply(1:ncol(fit3$coefficients), function(i) {
  res = topTable(fit3, coef=i, number=Inf, sort.by="none")
  res$SE = res$logFC / res$t
  res = res %>%
    tibble::rownames_to_column(var = "Protein") %>%
    mutate(Contrast = colnames(fit3$coefficients)[i])
  return(res)
})

limma_model = bind_rows(all_results)

fwrite(limma_model, file=file.path(output_dir, "limma_model.csv"))

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

annotation_nb = annotation[annotation$Condition != "Blank"]
targets = as.data.frame(annotation_nb)[c("R.FileName", "Condition")]
row.names(targets) = targets$R.FileName
limpa_elist = new("EList", list(E = limpa_dt,
                                genes = mapper,
                                targets = targets))

pdf(file.path(output_dir, "limpa_dpc_plots.pdf"), width = 12, height = 8)
dpcfit = dpc(limpa_elist)

y.protein = dpcQuant(limpa_elist, "PG.ProteinGroups",
                     dpc = dpcfit)

class = annotation_nb$Condition[
  match(colnames(y.protein$E), annotation_nb$R.FileName)] %>%
  as.factor()
design = model.matrix(~0+class)

fit = dpcDE(y.protein, design, plot=TRUE)
dev.off()

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

all_results = lapply(1:ncol(fit$coefficients), function(i) {
  res = topTable(fit, coef=i, number=Inf, sort.by="none")
  res$SE = res$logFC / res$t
  res = res %>%
    tibble::rownames_to_column(var = "Protein") %>%
    mutate(Contrast = colnames(fit$coefficients)[i])
  return(res)
})

limpa_model = bind_rows(all_results)

fwrite(limpa_model, file=file.path(output_dir, "limpa_model.csv"))

print("limpa finished")

# DEqMS ------------------------------------------------------------------------
deqms_input = merged_input %>% filter(
  F.ExcludedFromQuantification == FALSE &
    R.Condition != "blank" &
    F.PeakArea > 1 &
    PG.Qvalue < .01 &
    EG.Qvalue < .01)

deqms_input$Feature = paste(deqms_input$EG.PrecursorId,
                            deqms_input$F.FrgIon,
                            deqms_input$F.Charge, sep="_")

deqms_input = deqms_input %>%
  select(PG.ProteinGroups, Feature, R.FileName, F.PeakArea)
deqms_input$LogIntensities = log2(deqms_input$F.PeakArea)

deqms_input = deqms_input %>% group_by(PG.ProteinGroups, Feature, R.FileName) %>%
  summarize(LogIntensities=max(LogIntensities))
deqms_input = deqms_input %>% filter(PG.ProteinGroups != "")

deqms_input = pivot_wider(deqms_input,
                          id_cols=c("PG.ProteinGroups", "Feature"),
                          names_from="R.FileName",
                          values_from="LogIntensities")

# Adjusted DEqMS::medianSummary function without ref_col
deqms_summarize = function (dat, group_col = 2)
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

deqms_summarized = deqms_summarize(deqms_input, group_col=1)
pep_count = deqms_input %>% group_by(PG.ProteinGroups) %>%
  summarise(count=n_distinct(Feature),.groups = 'drop') %>%
  as.data.frame()

pep_count$count = pep_count$count+1
pep_count = as.data.frame(pep_count)
rownames(pep_count) = pep_count$PG.ProteinGroups
pep_count$PG.ProteinGroups = NULL

annotation_nb = annotation %>% filter(Condition != "Blank")

class = annotation_nb$Condition[
  match(colnames(deqms_summarized), annotation_nb$R.FileName)] %>%
  as.factor()
design = model.matrix(~0+class)

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

all_results = lapply(1:ncol(fit4$coefficients), function(i) {
  res = topTable(fit4, coef=i, number=Inf, sort.by="none")
  res$SE = res$logFC / res$t
  res = res %>%
    tibble::rownames_to_column(var = "Protein") %>%
    mutate(Contrast = colnames(fit4$coefficients)[i])
  return(res)
})

deqms_model = bind_rows(all_results)

fwrite(deqms_model, file=file.path(output_dir, "deqms_model.csv"))

print("DEqMS finished")

## Comparison table ------------------------------------------------------------
# Verbatim from Mixture_of_proteomes_processing.R lines 668-731.

protein_mapping = fread(file=file.path(data_folder, "idmapping.tsv"), sep="\t")

attach_organism = function(model_dt, by_x = "Protein"){
  merge(as.data.table(model_dt),
        protein_mapping %>% select(Entry, Organism),
        all.x = TRUE, all.y = FALSE,
        by.x = by_x, by.y = "Entry")
}

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

pooled_metrics = function(model_dt, pval_col, label_col, method, source_label){
  dt = attach_organism(model_dt)
  fdr = calculate_fdr(copy(dt), pval_col, label_col)
  tpr = calculate_tpr(copy(dt), pval_col, label_col)
  PPV = 1 - sum(fdr$n_false_positive) / max(sum(fdr$n_discovered), 1)
  TPR = sum(tpr$n_tp)            / max(sum(tpr$n_positive), 1)
  data.table(Method = method, Source = source_label,
             Mix_TPR = TPR, Mix_PPV = PPV)
}

method_order = c("MSstats+", "MSstats", "limpa",
                 "MaxLFQ + limma", "msqrob2", "DEqMS")

comparison_table = rbindlist(list(
  pooled_metrics(msstatsplus_model, "adj.pvalue", "Label",
                 "MSstats+", "log2 rerun"),
  pooled_metrics(msstats_model,     "adj.pvalue", "Label",
                 "MSstats",  "log2 rerun"),
  pooled_metrics(limpa_model,       "adj.P.Val",  "Contrast",
                 "limpa",    "log2 rerun"),
  pooled_metrics(limma_model,       "adj.P.Val",  "Contrast",
                 "MaxLFQ + limma", "log2 rerun"),
  pooled_metrics(msqrob2_model,     "adjPval",    "Label",
                 "msqrob2",  "log2 rerun"),
  pooled_metrics(deqms_model,       "adj.P.Val",  "Contrast",
                 "DEqMS",    "log2 rerun")
))

comparison_table[, Method := factor(Method, levels = method_order)]
setorder(comparison_table, Method)

rounded = copy(comparison_table)
rounded[, c("Mix_TPR", "Mix_PPV") := lapply(.SD, round, 3),
        .SDcols = c("Mix_TPR", "Mix_PPV")]

fwrite(rounded, file.path(output_dir, "Mix_of_Proteome_comparison_table.csv"))
writeLines(capture.output(print(rounded)),
           file.path(output_dir, "Mix_of_Proteome_comparison_table.txt"))

writeLines(c(
  "# Mix of Proteome v1_log2 rerun",
  "",
  "Variant label: log2 rerun",
  "",
  "All six methods (MSstats+, MSstats, msqrob2, MaxLFQ+limma, limpa, DEqMS)",
  "are re-fit from raw input here, using log2 as the per-feature intensity",
  "transform. Reproduces the Mixture-of-species column of Table 1.",
  "",
  "Output files:",
  "- MSstatsplus_model.csv, msstats_model.csv",
  "- msqrob2_model.csv, msqrob_obj.rda",
  "- limma_model.csv, qc-plots.pdf",
  "- limpa_model.csv, limpa_dpc_plots.pdf",
  "- deqms_model.csv",
  "- Mix_of_Proteome_comparison_table.csv / .txt",
  "",
  paste0("- R: ", getRversion()),
  paste0("- MSstats: ", as.character(utils::packageVersion("MSstats"))),
  paste0("- msqrob2: ", as.character(utils::packageVersion("msqrob2"))),
  paste0("- limma: ",   as.character(utils::packageVersion("limma"))),
  paste0("- limpa: ",   as.character(utils::packageVersion("limpa"))),
  paste0("- DEqMS: ",   as.character(utils::packageVersion("DEqMS")))
), file.path(output_dir, "README.md"))

print(rounded)
