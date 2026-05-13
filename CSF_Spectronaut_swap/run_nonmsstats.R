## Non-MSstats runner: msqrob2, limma, limpa, DEqMS, prolfqua.
##
## NORMALIZATION values that this script understands:
##   "none"     -> log2 transform only (no inter-sample normalization)
##   "vsn"      -> vsn::justvsn on the raw intensity scale
##   "quantile" -> log2 first, then limma::normalizeBetweenArrays(method="quantile")
## "equalizeMedians" is MSstats-only.
suppressPackageStartupMessages({
  library(QFeatures)
  library(msqrob2)
  library(iq)
  library(limma)
  library(limpa)
  library(DEqMS)
  library(prolfqua)
  library(vsn)
})
source("run_step_common.R")
source("run_prolfqua_step.R")

if (normalization == "equalizeMedians") {
  stop("run_nonmsstats.R does not implement NORMALIZATION=equalizeMedians ",
       "(MSstats-only). Use 'none', 'vsn', or 'quantile'.")
}
apply_quantile = (normalization == "quantile")

## msqrob2 (hurdle) ------------------------------------------------------------
# Port of the prolfquabenchmark workflow (vignettes/Benchmark_msqrob2.Rmd):
# precursor LFQData -> normalization at precursor scale ->
# LFQDataToSummarizedExperiment -> QFeatures -> aggregateFeatures with
# robustSummary (QFeatures default) -> msqrobHurdle (rlm + glm).
# preprocess timing covers prep + normalization + aggregation;
# model timing covers msqrobHurdle + hypothesisTestHurdle.
t_pre = tic()
ms2_input = prepare_data_for_limma(merged_input, all_proteins, no_swap)
ms2_input = as.data.frame(ms2_input)
ms2_input = ms2_input[is.finite(ms2_input$F.PeakArea) &
                        ms2_input$F.PeakArea > 0, , drop = FALSE]

config_ms2 = prolfqua::AnalysisConfiguration$new()
config_ms2$file_name = "R.FileName"
config_ms2$factors["group_"] = "Condition"
config_ms2$hierarchy[["protein_Id"]]   = "PG.ProteinGroups"
config_ms2$hierarchy[["precursor_Id"]] = "EG.PrecursorId"
config_ms2$hierarchy[["fragment_Id"]]  = c("Feature", "F.FrgLossType")
config_ms2$hierarchy_depth = 1
config_ms2$set_response("F.PeakArea")
adata_ms2 = prolfqua::setup_analysis(ms2_input, config_ms2)
lfq_pep = prolfqua::LFQData$new(adata_ms2, config_ms2)

tr_pep = lfq_pep$get_Transformer()
if (apply_vsn) {
  tr_pep$intensity_matrix(.func = vsn_normalize_matrix)
} else if (apply_quantile) {
  tr_pep$log2()
  # force=TRUE because prolfqua's transformer refuses to run on already-
  # transformed data; without it the quantile step silently no-ops.
  tr_pep$intensity_matrix(.func = quantile_normalize_log2_matrix,
                          force = TRUE)
} else {
  tr_pep$log2()
  tr_pep$robscale()
}
lfq_pep_norm = tr_pep$lfq

# precursor LFQData -> SummarizedExperiment -> QFeatures
se = prolfqua::LFQDataToSummarizedExperiment(lfqdata = lfq_pep_norm)
pe = QFeatures::QFeatures(list(peptide = se), colData = colData(se))

# Use the QFeatures default (MsCoreUtils::robustSummary) — vanilla msqrob2
# vignette behaviour. Robust-regression peptide -> protein rollup, no
# explicit medianPolish override.
pe = QFeatures::aggregateFeatures(
  pe, i = "peptide", fcol = "protein_Id", name = "protein"
)

pre_s = toc(t_pre)
t_mod = tic()
# Drop intercept so both factor levels get their own coefficient — this
# matches the prolfquabenchmark contrast spec (parameters named after
# the factor levels directly, not "(Intercept)" + offset).
prlm = msqrob2::msqrobHurdle(pe, i = "protein", formula = ~ 0 + group_,
                              overwrite = TRUE)
L_ms2 = msqrob2::makeContrast(
  "group_Condition2 - group_Condition1=0",
  parameterNames = c("group_Condition1", "group_Condition2")
)
prlm = msqrob2::hypothesisTestHurdle(prlm, i = "protein", L_ms2,
                                      overwrite = TRUE)
mod_s = toc(t_mod)
save(prlm, file = file.path(out_dir("msqrob2"), "msqrob_obj.rda"))

# Extract the (single) hurdle contrast DataFrame.
# Columns from hypothesisTestHurdle (msqrob2):
#   intensity (rlm) : logFC, logFCse, logFCdf, logFCt, logFCpval
#   count    (glm)  : logOR, logORse, logORdf, logORt, logORpval
#   combined        : fisher, fisherDf, fisherPval, fisherAdjPval
# Prefer the intensity model where it could fit; fall back to the count
# model for proteins observed in only one group.
xx = SummarizedExperiment::rowData(prlm[["protein"]])
hurdle_cols = grep("^hurdle_", names(xx), value = TRUE)
stopifnot(length(hurdle_cols) == 1)
hdf = as.data.frame(xx[[hurdle_cols[1]]])
hdf$Protein = rownames(xx)

use_intensity = !is.na(hdf$logFC)
msqrob2_model = data.frame(
  Protein = c(hdf$Protein[use_intensity], hdf$Protein[!use_intensity]),
  logFC   = c(hdf$logFC[use_intensity],   hdf$logOR[!use_intensity]),
  SE      = c(hdf$logFCse[use_intensity], hdf$logORse[!use_intensity]),
  DF      = c(hdf$logFCdf[use_intensity], hdf$logORdf[!use_intensity]),
  pvalue  = c(hdf$logFCpval[use_intensity], hdf$logORpval[!use_intensity]),
  source  = c(rep("intensity", sum(use_intensity)),
              rep("count",     sum(!use_intensity))),
  stringsAsFactors = FALSE
)
msqrob2_model$adj.pvalue = p.adjust(msqrob2_model$pvalue, method = "BH")
msqrob2_model = label_proteins(msqrob2_model)
fwrite(msqrob2_model, file = file.path(out_dir("msqrob2"),
                                          "msqrob2_model.csv"))
write_timing("msqrob2", out_dir("msqrob2"), pre_s, mod_s)
message("msqrob2 finished")

## limma -----------------------------------------------------------------------
t_pre = tic()
limma_long = prepare_data_for_limma(merged_input, all_proteins, no_swap)
maxlfq_input = preprocess(
  limma_long,
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
  # MaxLFQ output is log2; un-log2 to feed vsn linear intensities.
  maxlfq_summarized = vsn_normalize_matrix(2 ^ maxlfq_summarized)
} else if (apply_quantile) {
  # MaxLFQ output is already log2 -> quantile-normalize directly.
  maxlfq_summarized = quantile_normalize_log2_matrix(maxlfq_summarized)
}

class = annotation$Condition[
  match(colnames(maxlfq_summarized), annotation$R.FileName)
] |> as.factor()
design = model.matrix(~ 0 + class)

pre_s = toc(t_pre)
t_mod = tic()
fit1 = lmFit(maxlfq_summarized, design = design)
cont = makeContrasts(classCondition2 - classCondition1, levels = design)
fit2 = contrasts.fit(fit1, contrasts = cont)
fit3 = eBayes(fit2)
mod_s = toc(t_mod)

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
write_timing("limma", out_dir("limma"), pre_s, mod_s)
message("limma finished")

## limpa -----------------------------------------------------------------------
t_pre = tic()
limpa_input_wide = prepare_data_for_limpa(merged_input, all_proteins, no_swap)
mapper = limpa_input_wide[c("PG.ProteinGroups", "Feature")]
row.names(mapper) = mapper$Feature

limpa_dt = copy(limpa_input_wide)
row.names(limpa_dt) = limpa_dt$Feature
limpa_dt = limpa_dt[, !colnames(limpa_dt) %in% c("PG.ProteinGroups", "Feature")]
# shared_input_limpa is on the LINEAR scale (raw F.PeakArea pivoted wide).
if (apply_vsn) {
  limpa_dt = as.data.frame(vsn_normalize_matrix(limpa_dt))
} else if (apply_quantile) {
  # log2 the raw matrix, then quantile-normalize on log2.
  limpa_dt = as.data.frame(
    quantile_normalize_log2_matrix(log2(as.matrix(limpa_dt)))
  )
}

targets = as.data.frame(annotation)[c("R.FileName", "Condition")]
row.names(targets) = targets$R.FileName

limpa_elist = new("EList",
                   list(E = limpa_dt, genes = mapper, targets = targets))
pre_s = toc(t_pre)
t_mod = tic()
dpcfit = tryCatch(dpc(limpa_elist),
                   error = function(e) {
                     warning("limpa dpc() failed under NORMALIZATION=",
                             normalization, ": ", conditionMessage(e),
                             "; writing empty model output.")
                     NULL
                   })
if (!is.null(dpcfit)) {
  y.protein = dpcQuant(limpa_elist, "PG.ProteinGroups", dpc = dpcfit)
  class = annotation$Condition[
    match(colnames(y.protein$E), annotation$R.FileName)
  ] |> as.factor()
  design = model.matrix(~ 0 + class)
  fit = dpcDE(y.protein, design, plot = FALSE)
  cont = makeContrasts(classCondition2 - classCondition1, levels = design)
  fit = contrasts.fit(fit, contrasts = cont)
  fit = eBayes(fit)
  mod_s = toc(t_mod)
  limpa_model = data.frame(
    Protein = rownames(fit$coefficients),
    logFC   = as.numeric(fit$coefficients),
    SE      = as.numeric(sqrt(fit$s2.post) * fit$stdev.unscaled),
    DF      = fit$df.total,
    pvalue  = as.numeric(fit$p.value)
  )
  limpa_model$adj.pvalue = p.adjust(limpa_model$pvalue, method = "BH")
  limpa_model = label_proteins(limpa_model)
  fwrite(limpa_model, file = file.path(out_dir("limpa"),
                                          "limpa_model.csv"))
  write_timing("limpa", out_dir("limpa"), pre_s, mod_s)
  message("limpa finished")
} else {
  mod_s = toc(t_mod)
  write_timing("limpa", out_dir("limpa"), pre_s, mod_s)
  message("limpa skipped (dpc failed)")
}

## DEqMS -----------------------------------------------------------------------
t_pre = tic()
shared_input_deqms = prepare_data_for_deqms(merged_input, all_proteins, no_swap)
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

deqms_summarized = summarize_deqms_no_ref_col(shared_input_deqms, group_col = 1)
# summarize_deqms_no_ref_col emits a log2-scale, median-summarized protein
# matrix (prepare_data_for_deqms applies log2 upstream).
if (apply_vsn) {
  # vsn wants raw -> undo the log2 first then feed vsn.
  raw_protein = 2 ^ as.matrix(deqms_summarized)
  deqms_summarized = as.data.frame(vsn_normalize_matrix(raw_protein))
} else if (apply_quantile) {
  # Already log2 -> quantile-normalize directly.
  deqms_summarized = as.data.frame(
    quantile_normalize_log2_matrix(deqms_summarized)
  )
}
pep_count = shared_input_deqms |>
  dplyr::group_by(PG.ProteinGroups) |>
  dplyr::summarise(count = dplyr::n_distinct(Feature) + 1, .groups = "drop") |>
  as.data.frame()
rownames(pep_count) = pep_count$PG.ProteinGroups
pep_count$PG.ProteinGroups = NULL

class = annotation$Condition[
  match(colnames(deqms_summarized), annotation$R.FileName)
] |> as.factor()
design = model.matrix(~ 0 + class)

pre_s = toc(t_pre)
t_mod = tic()
fit1 = lmFit(deqms_summarized, design = design)
cont = makeContrasts(classCondition2 - classCondition1, levels = design)
fit2 = contrasts.fit(fit1, contrasts = cont)
fit3 = eBayes(fit2)
fit3$count = pep_count[rownames(fit3$coefficients), "count"]
fit4 = spectraCounteBayes(fit3)
mod_s = toc(t_mod)

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
write_timing("DEqMS", out_dir("DEqMS"), pre_s, mod_s)
message("DEqMS finished")

## prolfqua --------------------------------------------------------------------
# run_prolfqua_step.R is the single source of truth for the prolfqua DE block
# (config -> log+medpolish+exp -> protein-level normalization ->
# ContrastsLMImputeFacade). Returns model + per-phase timings.
prolfqua_res = run_prolfqua_step(
  merged_input, annotation, all_proteins, no_swap,
  normalization = normalization,
  vsn_func      = vsn_normalize_matrix,
  quantile_func = quantile_normalize_log2_matrix
)
prolfqua_model = label_proteins(prolfqua_res$model)
fwrite(prolfqua_model,
        file = file.path(out_dir("prolfqua"), "prolfqua_model.csv"))
write_timing("prolfqua", out_dir("prolfqua"),
             prolfqua_res$preprocess_seconds, prolfqua_res$model_seconds)
message("prolfqua finished")
