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

## Pre-compute the shared inputs ----------------------------------------------
shared_input_msqrob = prepare_data_for_msqrob(merged_input, all_proteins, no_swap)
shared_input_limma  = prepare_data_for_limma(merged_input, all_proteins, no_swap)
shared_input_limpa  = prepare_data_for_limpa(merged_input, all_proteins, no_swap)
shared_input_deqms  = prepare_data_for_deqms(merged_input, all_proteins, no_swap)

## msqrob2 --------------------------------------------------------------------
df_LFQ = dcast(shared_input_msqrob, ProteinName + Fragment ~ Run,
                value.var = "F.PeakArea", fun.aggregate = max, fill = NA)

colnames(df_LFQ)[3:ncol(df_LFQ)] = paste("F.PeakArea",
                                          colnames(df_LFQ)[3:ncol(df_LFQ)],
                                          sep = "_")
ecols = grep("F.PeakArea", colnames(df_LFQ))

pe = readQFeatures(df_LFQ, fnames = 2, quantCols = ecols, name = "peptideRaw")
colData(pe)$condition = annotation$Condition[
  match(sub("^F\\.PeakArea_", "", colnames(df_LFQ)[ecols]),
        annotation$R.FileName)
] |> as.factor()
rowData(pe[["peptideRaw"]])$nNonZero = rowSums(assay(pe[["peptideRaw"]]) > 0,
                                                  na.rm = TRUE)
pe = zeroIsNA(pe, "peptideRaw")
pe = logTransform(pe, base = 2, i = "peptideRaw", name = "peptideLog")
if (apply_vsn) {
  # vsn::justvsn applied to raw (linear) precursors.
  assay(pe[["peptideLog"]]) = vsn_normalize_matrix(assay(pe[["peptideRaw"]]))
} else if (apply_quantile) {
  # quantile on the log2-scale precursor matrix.
  assay(pe[["peptideLog"]]) = quantile_normalize_log2_matrix(
    assay(pe[["peptideLog"]])
  )
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
maxlfq_input = preprocess(
  shared_input_limma,
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
mapper = shared_input_limpa[c("PG.ProteinGroups", "Feature")]
row.names(mapper) = mapper$Feature

limpa_dt = copy(shared_input_limpa)
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
  message("limpa finished")
} else {
  message("limpa skipped (dpc failed)")
}

## DEqMS -----------------------------------------------------------------------
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
prolfqua_model = run_prolfqua_step(
  merged_input, annotation, all_proteins, no_swap,
  normalization = normalization,
  vsn_func      = vsn_normalize_matrix,
  quantile_func = quantile_normalize_log2_matrix
)
prolfqua_model = label_proteins(prolfqua_model)
fwrite(prolfqua_model,
        file = file.path(out_dir("prolfqua"), "prolfqua_model.csv"))
message("prolfqua finished")
