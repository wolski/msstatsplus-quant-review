## Standalone msqrob2 runner — re-fits only the msqrob2 (hurdle) block
## without touching limma / limpa / DEqMS / prolfqua / MSstats. Same
## env-var contract as run_nonmsstats.R.
suppressPackageStartupMessages({
  library(QFeatures)
  library(msqrob2)
  library(SummarizedExperiment)
  library(prolfqua)
  library(vsn)
  library(limma)
})
source("run_step_common.R")

## msqrob2 (hurdle) ------------------------------------------------------------
# Mirror of the msqrob2 block in run_nonmsstats.R:
#   precursors -> log2
#             -> SE -> QFeatures -> aggregateFeatures(robustSummary)
#             -> normalize protein matrix (none | vsn | quantile)
#             -> msqrobHurdle + hypothesisTestHurdle
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
tr_pep$log2()
lfq_pep_log = tr_pep$lfq

se = prolfqua::LFQDataToSummarizedExperiment(lfqdata = lfq_pep_log)
pe = QFeatures::QFeatures(list(peptide = se), colData = colData(se))
pe = QFeatures::normalize(pe, i = "peptide", method = "center.median",
                          name = "peptide_norm")
pe = QFeatures::aggregateFeatures(
  pe, i = "peptide_norm", fcol = "protein_Id", name = "protein"
)

apply_vsn      = (normalization == "vsn")
apply_quantile = (normalization == "quantile")
prot_mat = SummarizedExperiment::assay(pe[["protein"]])
if (apply_vsn) {
  prot_mat = vsn_normalize_matrix(2 ^ prot_mat)
} else if (apply_quantile) {
  prot_mat = quantile_normalize_log2_matrix(prot_mat)
}
SummarizedExperiment::assay(pe[["protein"]]) = prot_mat

pre_s = toc(t_pre)
t_mod = tic()
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
fwrite(msqrob2_model, file = file.path(out_dir("msqrob2"), "msqrob2_model.csv"))
write_timing("msqrob2", out_dir("msqrob2"), pre_s, mod_s)
message("msqrob2 finished")
