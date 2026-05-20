## Non-MSstats runner for the original CSF Spectronaut publication benchmark.
## Methods: MaxLFQ + limma, limpa, DEqMS, prolfqua. Normalization is log2 only.
## msqrob2 is skipped by default because it is slow; set RUN_MSQROB2=1 to run it.

suppressPackageStartupMessages({
  library(iq)
  library(limma)
  library(limpa)
  library(DEqMS)
  library(prolfqua)
})
source("run_step_common.R")
source("run_prolfqua_step.R")

if (variant != "V1_log2" || normalization != "none") {
  stop("CSF_Spectronaut/run_nonmsstats.R currently supports only ",
       "VARIANT=V1_log2 and NORMALIZATION=none.")
}

true_negatives = protein_swap_list[Label == "Negative", Protein]

if (identical(Sys.getenv("RUN_MSQROB2"), "1")) {
  suppressPackageStartupMessages({
    library(QFeatures)
    library(msqrob2)
  })
  source("run_msqrob2_step.R")
} else {
  message("msqrob2 skipped; set RUN_MSQROB2=1 to run it")
}

## limma -----------------------------------------------------------------------
t_pre = tic()
limma_input = prepare_data_for_limma(merged_input, true_positives,
                                     true_negatives)
maxlfq_input = preprocess(
  limma_input,
  primary_id = "PG.ProteinGroups",
  secondary_id = c("Feature"),
  sample_id = "R.FileName",
  intensity_col = "F.PeakArea",
  median_normalization = FALSE,
  log2_intensity_cutoff = 0,
  pdf_out = file.path(out_dir("MaxLFQ_limma"), "qc-plots.pdf"),
  pdf_width = 12,
  pdf_height = 8,
  intensity_col_sep = NULL,
  intensity_col_id = NULL,
  na_string = "0"
)
maxlfq_summarized = fast_MaxLFQ(maxlfq_input)$estimate

class = annotation$Condition[
  match(colnames(maxlfq_summarized), annotation$R.FileName)
] |> as.factor()
design = model.matrix(~ 0 + class)

pre_s = toc(t_pre)
t_mod = tic()
fit1 = limma::lmFit(maxlfq_summarized, design = design)
cont = limma::makeContrasts(classCondition2 - classCondition1,
                            levels = design)
fit2 = limma::contrasts.fit(fit1, contrasts = cont)
fit3 = limma::eBayes(fit2)
mod_s = toc(t_mod)

limma_model = data.frame(
  Protein = rownames(fit3$coefficients),
  logFC = as.numeric(fit3$coefficients),
  SE = as.numeric(sqrt(fit3$s2.post) * fit3$stdev.unscaled),
  DF = fit3$df.total,
  pvalue = as.numeric(fit3$p.value)
)
limma_model$adj.pvalue = p.adjust(limma_model$pvalue, method = "BH")
limma_model = label_proteins(limma_model)
fwrite(limma_model, file = file.path(out_dir("MaxLFQ_limma"), "limma_model.csv"))
write_timing("limma", out_dir("MaxLFQ_limma"), pre_s, mod_s)
message("limma finished")

## limpa -----------------------------------------------------------------------
t_pre = tic()
limpa_input = prepare_data_for_limpa(merged_input, true_positives,
                                     true_negatives)
mapper = limpa_input[c("PG.ProteinGroups", "Feature")]
row.names(mapper) = mapper$Feature

limpa_dt = data.table::copy(limpa_input)
row.names(limpa_dt) = limpa_dt$Feature
limpa_dt = limpa_dt[, !colnames(limpa_dt) %in% c("PG.ProteinGroups",
                                                 "Feature")]

targets = as.data.frame(annotation)[c("R.FileName", "Condition")]
row.names(targets) = targets$R.FileName
limpa_elist = new("EList", list(E = limpa_dt, genes = mapper,
                                targets = targets))

pre_s = toc(t_pre)
t_mod = tic()
dpcfit = limpa::dpc(limpa_elist)
y_protein = limpa::dpcQuant(limpa_elist, "PG.ProteinGroups", dpc = dpcfit)
class = annotation$Condition[
  match(colnames(y_protein$E), annotation$R.FileName)
] |> as.factor()
design = model.matrix(~ 0 + class)
fit = limpa::dpcDE(y_protein, design, plot = TRUE)
cont = limma::makeContrasts(classCondition2 - classCondition1,
                            levels = design)
fit = limma::contrasts.fit(fit, contrasts = cont)
fit = limma::eBayes(fit)
mod_s = toc(t_mod)

limpa_model = data.frame(
  Protein = rownames(fit$coefficients),
  logFC = as.numeric(fit$coefficients),
  SE = as.numeric(sqrt(fit$s2.post) * fit$stdev.unscaled),
  DF = fit$df.total,
  pvalue = as.numeric(fit$p.value)
)
limpa_model$adj.pvalue = p.adjust(limpa_model$pvalue, method = "BH")
limpa_model = label_proteins(limpa_model)
fwrite(limpa_model, file = file.path(out_dir("limpa"), "limpa_model.csv"))
write_timing("limpa", out_dir("limpa"), pre_s, mod_s)
message("limpa finished")

## DEqMS -----------------------------------------------------------------------
t_pre = tic()
deqms_input = prepare_data_for_deqms(merged_input, true_positives,
                                     true_negatives)

summarize_deqms_no_ref_col = function(dat, group_col = 2) {
  dat_ratio = dat
  dat_ratio[, 3:ncol(dat)] = dat_ratio[, 3:ncol(dat)] -
    matrixStats::rowMedians(as.matrix(dat_ratio[, 3:ncol(dat)]),
                            na.rm = TRUE)
  dat_summary = plyr::ddply(
    dat_ratio,
    colnames(dat)[group_col],
    function(x) matrixStats::colMedians(as.matrix(x[, 3:ncol(dat)]),
                                        na.rm = TRUE)
  )
  colnames(dat_summary)[2:ncol(dat_summary)] = colnames(dat)[3:ncol(dat)]
  dat_new = dat_summary[, -1]
  rownames(dat_new) = dat_summary[, 1]
  dat_new
}

deqms_summarized = summarize_deqms_no_ref_col(deqms_input, group_col = 1)
pep_count = deqms_input |>
  dplyr::group_by(PG.ProteinGroups) |>
  dplyr::summarise(count = dplyr::n_distinct(Feature) + 1,
                   .groups = "drop") |>
  as.data.frame()
rownames(pep_count) = pep_count$PG.ProteinGroups
pep_count$PG.ProteinGroups = NULL

class = annotation$Condition[
  match(colnames(deqms_summarized), annotation$R.FileName)
] |> as.factor()
design = model.matrix(~ 0 + class)

pre_s = toc(t_pre)
t_mod = tic()
fit1 = limma::lmFit(deqms_summarized, design = design)
cont = limma::makeContrasts(classCondition2 - classCondition1,
                            levels = design)
fit2 = limma::contrasts.fit(fit1, contrasts = cont)
fit3 = limma::eBayes(fit2)
fit3$count = pep_count[rownames(fit3$coefficients), "count"]
fit4 = DEqMS::spectraCounteBayes(fit3)
mod_s = toc(t_mod)

deqms_model = data.frame(
  Protein = rownames(fit4$coefficients),
  logFC = as.numeric(fit4$coefficients),
  SE = as.numeric(sqrt(fit4$s2.post) * fit4$stdev.unscaled),
  DF = fit4$df.total,
  pvalue = as.numeric(fit4$p.value)
)
deqms_model$adj.pvalue = p.adjust(deqms_model$pvalue, method = "BH")
deqms_model = label_proteins(deqms_model)
fwrite(deqms_model, file = file.path(out_dir("DEqMS"), "deqms_model.csv"))
write_timing("DEqMS", out_dir("DEqMS"), pre_s, mod_s)
message("DEqMS finished")

## prolfqua --------------------------------------------------------------------
prolfqua_res = run_prolfqua_step(
  merged_input, annotation, true_positives, true_negatives,
  normalization = normalization
)
prolfqua_model = label_proteins(prolfqua_res$model)
fwrite(prolfqua_model,
       file = file.path(out_dir("prolfqua"), "prolfqua_model.csv"))
write_timing("prolfqua", out_dir("prolfqua"),
             prolfqua_res$preprocess_seconds, prolfqua_res$model_seconds)
message("prolfqua finished")
