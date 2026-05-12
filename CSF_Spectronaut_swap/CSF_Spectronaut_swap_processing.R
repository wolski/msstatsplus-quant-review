## CSF Spectronaut precursor-swap benchmark — processing
##
## Adapted from ../CSF_Spectronaut/CSF_Spectronaut_processing.R. The swap is
## already baked into the input TSV (intensities swapped between paired
## precursors in G2 runs), so we DO NOT call swap_condition_labels() or pass
## proteins to the run-relabel logic inside the non-MSstats prepare_*()
## helpers. We simply compare Condition1 (G1) vs Condition2 (G2).
##
## Environment variables that select the run:
##   REPORT_PATH        : path to TSV (defaults to swapped report in this dir)
##   VARIANT            : "V1_log2" (default; log2 normalization) or "v2_vsn"
##                        (vsn::justvsn for non-MSstats methods; MSstats
##                        variants are skipped because they normalize
##                        internally).
##   OUT_SUFFIX         : "" (post-swap) or "_preswap" (against original TSV)
##   EXCLUDE_DILUTIONS  : comma-separated R.Condition values to drop from
##                        analysis (e.g. "1to32,1to64"). Default: none.
##   OUT_TAG            : suffix appended to the variant directory name so
##                        results from filtered analyses don't overwrite
##                        the canonical V1_log2/ / v2_vsn/ outputs
##                        (e.g. "_no_high_dilutions").
##
## Outputs land in <VARIANT><OUT_TAG>/<METHOD><OUT_SUFFIX>/. Run from this dir.

## Packages --------------------------------------------------------------------
suppressPackageStartupMessages({
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
  library(prolfqua)
  library(vsn)
})

source("../benchmark_experiments_functions.R")
source("run_prolfqua_step.R")

## Configuration ---------------------------------------------------------------
report_path = Sys.getenv(
  "REPORT_PATH",
  unset = "20250130_163144_CSF dilutions Jan 2025 no normalization_Report.tsv"
)
variant    = Sys.getenv("VARIANT",    unset = "V1_log2")
out_suffix = Sys.getenv("OUT_SUFFIX", unset = "")
out_tag    = Sys.getenv("OUT_TAG",    unset = "")
exclude_dilutions = Sys.getenv("EXCLUDE_DILUTIONS", unset = "")
exclude_dilutions = if (nchar(exclude_dilutions) > 0) {
  trimws(strsplit(exclude_dilutions, ",", fixed = TRUE)[[1]])
} else {
  character(0)
}

stopifnot(variant %in% c("V1_log2", "v2_vsn"))
apply_vsn = (variant == "v2_vsn")

cat(sprintf("[config] report      = %s\n", report_path))
cat(sprintf("[config] variant     = %s (vsn=%s)\n", variant, apply_vsn))
cat(sprintf("[config] out_suffix  = '%s'\n", out_suffix))
cat(sprintf("[config] out_tag     = '%s'\n", out_tag))
cat(sprintf("[config] exclude     = %s\n",
            if (length(exclude_dilutions) > 0)
              paste(exclude_dilutions, collapse = ", ") else "(none)"))

data_folder = "."
variant_dir = paste0(variant, out_tag)
dir.create(variant_dir, recursive = TRUE, showWarnings = FALSE)
out_dir = function(method) {
  d = file.path(variant_dir, paste0(method, out_suffix))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

## vsn helper ------------------------------------------------------------------
vsn_normalize_matrix = function(m) {
  m = as.matrix(m)
  storage.mode(m) = "double"
  dn = dimnames(m)
  m[!is.finite(m)] = NA_real_
  out = vsn::justvsn(m)
  dimnames(out) = dn
  out
}

## Load + prepare inputs -------------------------------------------------------
raw_input  = fread(report_path, sep = "\t")
annotation = fread(file.path(data_folder, "CSF_annotation.csv"))

raw_input  = raw_input[tolower(raw_input$R.Condition) != "blank", ]
annotation = annotation[tolower(annotation$Condition) != "blank", ]
if (length(exclude_dilutions) > 0) {
  drop_runs = raw_input[R.Condition %in% exclude_dilutions, unique(R.FileName)]
  raw_input  = raw_input[!R.Condition %in% exclude_dilutions, ]
  annotation = annotation[!R.FileName %in% drop_runs, ]
  cat(sprintf("[config] dropped %d runs in dilutions {%s}\n",
              length(drop_runs), paste(exclude_dilutions, collapse = ", ")))
}
annotation$Run = annotation$R.FileName
run_order = unique(annotation[, .(Run, Order)])

merged_input = merge(raw_input, annotation, by = "R.FileName",
                     all.x = TRUE, all.y = FALSE)

protein_swap_list = fread(file = "CSF_protein_swap_list.csv")
true_positives = protein_swap_list[Label == "Positive", Protein]
true_negatives = protein_swap_list[Label == "Negative", Protein]

# `insig = character(0)` disables the run-label flip in the helpers; passing
# all proteins as `sig` keeps them all in the output.
all_proteins = protein_swap_list$Protein
no_swap = character(0)

label_proteins = function(df, protein_col = "Protein") {
  df$Label = ifelse(df[[protein_col]] %in% true_positives,
                     "Positive", "Negative")
  df
}

## MSstats+ (V1_log2 only) -----------------------------------------------------
if (!apply_vsn) {
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
    normalization = FALSE,
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

  comparison = matrix(c(-1, 1), nrow = 1)
  rownames(comparison) = "Condition2-Condition1"
  colnames(comparison) = c("Condition1", "Condition2")

  msstatsplus_model = groupComparison(comparison, summarized,
                                       numberOfCores = 12)$ComparisonResult
  msstatsplus_model = label_proteins(msstatsplus_model)
  fwrite(msstatsplus_model, file = file.path(out_dir("MSstats+"),
                                                "MSstats+_model.csv"))
  message("MSstats+ finished")

  ## MSstats (V1_log2 only) ----------------------------------------------------
  base_msstats_input = MSstatsConvert::SpectronauttoMSstatsFormat(
    raw_input, annotation,
    intensity = "PeakArea",
    excludedFromQuantificationFilter = TRUE,
    filter_with_Qvalue = TRUE
  )
  fwrite(base_msstats_input, file = file.path(out_dir("MSstats"),
                                                  "MSstats_input.csv"))

  base_msstats_input = as.data.frame(base_msstats_input) %>%
    filter(Condition != "Blank")

  base_msstats_summarized = dataProcess(
    base_msstats_input,
    normalization = FALSE,
    featureSubset = "topN",
    n_top_feature = 100,
    MBimpute = TRUE,
    summaryMethod = "TMP",
    numberOfCores = 12
  )
  save(base_msstats_summarized, file = file.path(out_dir("MSstats"),
                                                    "MSstats_summarized.rda"))

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
}

## msqrob2 --------------------------------------------------------------------
df_LFQ = prepare_data_for_msqrob(merged_input, all_proteins, no_swap)
df_LFQ = dcast(df_LFQ, ProteinName + Fragment ~ Run,
                value.var = "F.PeakArea", fun.aggregate = max, fill = NA)

colnames(df_LFQ)[3:ncol(df_LFQ)] = paste("F.PeakArea",
                                          colnames(df_LFQ)[3:ncol(df_LFQ)],
                                          sep = "_")
ecols = grep("F.PeakArea", colnames(df_LFQ))

pe = readQFeatures(df_LFQ, fnames = 2, quantCols = ecols, name = "peptideRaw")
colData(pe)$condition = annotation$Condition[
  match(str_remove(colnames(df_LFQ)[ecols], "F.PeakArea_"),
        annotation$R.FileName)
] %>% as.factor()
rowData(pe[["peptideRaw"]])$nNonZero = rowSums(assay(pe[["peptideRaw"]]) > 0,
                                                  na.rm = TRUE)
pe = zeroIsNA(pe, "peptideRaw")
pe = logTransform(pe, base = 2, i = "peptideRaw", name = "peptideLog")
if (apply_vsn) {
  assay(pe[["peptideLog"]]) = vsn_normalize_matrix(assay(pe[["peptideRaw"]]))
}

Protein_filter = rowData(pe[["peptideLog"]])$ProteinName %in%
  smallestUniqueGroups(rowData(pe[["peptideLog"]])$ProteinName)
pe = pe[Protein_filter, ]
pe = filterFeatures(pe, ~ nNonZero >= 2)
pe = QFeatures::impute(pe, i = "peptideLog", name = "peptideImp",
                        method = "QRILC")
pe = aggregateFeatures(pe, i = "peptideImp", fcol = "ProteinName",
                        name = "protein")

pe = msqrob(object = pe, i = "protein", formula = ~condition)
L = makeContrast("conditionCondition2=0",
                  parameterNames = c("conditionCondition2",
                                      "conditionCondition1"))
pe = hypothesisTest(object = pe, i = "protein", contrast = L)
save(pe, file = file.path(out_dir("msqrob2"), "msqrob_obj.rda"))

msqrob2_model = rowData(pe[["protein"]])$`conditionCondition2`
msqrob2_model$Protein = rownames(msqrob2_model)
msqrob2_model = label_proteins(msqrob2_model)
fwrite(msqrob2_model, file = file.path(out_dir("msqrob2"),
                                          "msqrob2_model.csv"))
message("msqrob2 finished")

## limma -----------------------------------------------------------------------
limma_input = prepare_data_for_limma(merged_input, all_proteins, no_swap)
maxlfq_input = preprocess(
  limma_input,
  primary_id = "PG.ProteinGroups",
  secondary_id = c("Feature"),
  sample_id = "R.FileName",
  intensity_col = "F.PeakArea",
  median_normalization = FALSE,
  log2_intensity_cutoff = 0,
  pdf_out = file.path(out_dir("limma"), "qc-plots.pdf"),
  pdf_width = 12, pdf_height = 8,
  intensity_col_sep = NULL,
  intensity_col_id = NULL,
  na_string = "0"
)
maxlfq_summarized = fast_MaxLFQ(maxlfq_input)$estimate
if (apply_vsn) {
  # MaxLFQ output is log2; v2_vsn expects raw-scale input to vsn -> reverse
  # the log2 to feed vsn the linear intensities (matches Mix_of_Proteome
  # v2_vsn convention).
  maxlfq_summarized = vsn_normalize_matrix(2 ^ maxlfq_summarized)
}

class = annotation$Condition[
  match(colnames(maxlfq_summarized), annotation$R.FileName)
] %>% as.factor()
design = model.matrix(~ 0 + class)

fit1 = lmFit(maxlfq_summarized, design = design)
cont = makeContrasts(classCondition2 - classCondition1, levels = design)
fit2 = contrasts.fit(fit1, contrasts = cont)
fit3 = eBayes(fit2)

limma_model = data.frame(
  Protein = rownames(fit3$coefficients),
  logFC   = as.numeric(fit3$coefficients),
  SE      = as.numeric(sqrt(fit3$s2.post) * fit3$stdev.unscaled),
  DF      = fit3$df.total,
  pvalue  = as.numeric(fit3$p.value)
)
limma_model$adj.pvalue = p.adjust(limma_model$pvalue, method = "BH")
limma_model = label_proteins(limma_model)
fwrite(limma_model, file = file.path(out_dir("limma"), "limma_model.csv"))
message("limma finished")

## limpa -----------------------------------------------------------------------
limpa_input = prepare_data_for_limpa(merged_input, all_proteins, no_swap)
mapper = limpa_input[c("PG.ProteinGroups", "Feature")]
row.names(mapper) = mapper$Feature

limpa_dt = copy(limpa_input)
row.names(limpa_dt) = limpa_dt$Feature
limpa_dt = limpa_dt[, !colnames(limpa_dt) %in% c("PG.ProteinGroups", "Feature")]
if (apply_vsn) {
  limpa_dt = as.data.frame(vsn_normalize_matrix(limpa_dt))
}

targets = as.data.frame(annotation)[c("R.FileName", "Condition")]
row.names(targets) = targets$R.FileName

limpa_elist = new("EList",
                   list(E = limpa_dt, genes = mapper, targets = targets))
dpcfit = dpc(limpa_elist)
y.protein = dpcQuant(limpa_elist, "PG.ProteinGroups", dpc = dpcfit)

class = annotation$Condition[
  match(colnames(y.protein$E), annotation$R.FileName)
] %>% as.factor()
design = model.matrix(~ 0 + class)

fit = dpcDE(y.protein, design, plot = FALSE)
cont = makeContrasts(classCondition2 - classCondition1, levels = design)
fit = contrasts.fit(fit, contrasts = cont)
fit = eBayes(fit)

limpa_model = data.frame(
  Protein = rownames(fit$coefficients),
  logFC   = as.numeric(fit$coefficients),
  SE      = as.numeric(sqrt(fit$s2.post) * fit$stdev.unscaled),
  DF      = fit$df.total,
  pvalue  = as.numeric(fit$p.value)
)
limpa_model$adj.pvalue = p.adjust(limpa_model$pvalue, method = "BH")
limpa_model = label_proteins(limpa_model)
fwrite(limpa_model, file = file.path(out_dir("limpa"), "limpa_model.csv"))
message("limpa finished")

## DEqMS -----------------------------------------------------------------------
deqms_input = prepare_data_for_deqms(merged_input, all_proteins, no_swap)

summarize_deqms_no_ref_col = function(dat, group_col = 2) {
  dat.ratio = dat
  dat.ratio[, 3:ncol(dat)] = dat.ratio[, 3:ncol(dat)] -
    matrixStats::rowMedians(as.matrix(dat.ratio[, 3:ncol(dat)]), na.rm = TRUE)
  dat.summary = plyr::ddply(
    dat.ratio, colnames(dat)[group_col],
    function(x) matrixStats::colMedians(as.matrix(x[, 3:ncol(dat)]),
                                          na.rm = TRUE)
  )
  colnames(dat.summary)[2:ncol(dat.summary)] = colnames(dat)[3:ncol(dat)]
  dat.new = dat.summary[, -1]
  rownames(dat.new) = dat.summary[, 1]
  dat.new
}

deqms_summarized = summarize_deqms_no_ref_col(deqms_input, group_col = 1)
if (apply_vsn) {
  # prepare_data_for_deqms applies log2 upstream; vsn expects raw scale.
  raw_protein = 2 ^ as.matrix(deqms_summarized)
  deqms_summarized = as.data.frame(vsn_normalize_matrix(raw_protein))
}
pep_count = deqms_input %>%
  group_by(PG.ProteinGroups) %>%
  summarise(count = n_distinct(Feature) + 1, .groups = "drop") %>%
  as.data.frame()
rownames(pep_count) = pep_count$PG.ProteinGroups
pep_count$PG.ProteinGroups = NULL

class = annotation$Condition[
  match(colnames(deqms_summarized), annotation$R.FileName)
] %>% as.factor()
design = model.matrix(~ 0 + class)

fit1 = lmFit(deqms_summarized, design = design)
cont = makeContrasts(classCondition2 - classCondition1, levels = design)
fit2 = contrasts.fit(fit1, contrasts = cont)
fit3 = eBayes(fit2)
fit3$count = pep_count[rownames(fit3$coefficients), "count"]
fit4 = spectraCounteBayes(fit3)

deqms_model = data.frame(
  Protein = rownames(fit4$coefficients),
  logFC   = as.numeric(fit4$coefficients),
  SE      = as.numeric(sqrt(fit4$s2.post) * fit4$stdev.unscaled),
  DF      = fit4$df.total,
  pvalue  = as.numeric(fit4$p.value)
)
deqms_model$adj.pvalue = p.adjust(deqms_model$pvalue, method = "BH")
deqms_model = label_proteins(deqms_model)
fwrite(deqms_model, file = file.path(out_dir("DEqMS"), "deqms_model.csv"))
message("DEqMS finished")

## prolfqua --------------------------------------------------------------------
# Default lm model with imputation via ContrastsLMImputeFacade: lm fit per
# protein; proteins with NA coefficients (all missing in one group) are
# re-fit after LOD imputation with a borrowed covariance matrix, then
# moderated. See prolfqua/R/ContrastsFacades.R.
prolfqua_model = run_prolfqua_step(merged_input, annotation, all_proteins,
                                    no_swap, apply_vsn,
                                    vsn_func = vsn_normalize_matrix)
prolfqua_model = label_proteins(prolfqua_model)
fwrite(prolfqua_model, file = file.path(out_dir("prolfqua"),
                                          "prolfqua_model.csv"))
message("prolfqua finished")
