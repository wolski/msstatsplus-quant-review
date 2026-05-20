## Optional msqrob2 step for the standalone protein-swap non-MSstats runner.
## Sourced from CSF_Spectronaut_protein_swap/run_nonmsstats.R only when
## RUN_MSQROB2=1.

## msqrob2 (hurdle) ------------------------------------------------------------
# Pipeline:
#   precursor LFQData -> log2 -> SummarizedExperiment -> QFeatures
#     -> aggregateFeatures(robustSummary)  (vanilla msqrob2 vignette default)
#     -> normalize the PROTEIN-LEVEL matrix (vsn / quantile / robscale)
#     -> msqrobHurdle + hypothesisTestHurdle
# Normalization is applied AT THE PROTEIN LEVEL, matching prolfqua and the
# established proteomics convention. Quantile / vsn on the precursor
# matrix breaks robustSummary's per-protein rlm rollup, so the only
# place where it's safe (and meaningful) to switch the normalization in
# this pipeline is post-aggregation.
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

# log2 at precursor level only (variant-specific normalization happens
# at the protein level below).
tr_pep = lfq_pep$get_Transformer()
tr_pep$log2()
lfq_pep_log = tr_pep$lfq

# precursor LFQData -> SummarizedExperiment -> QFeatures
se = prolfqua::LFQDataToSummarizedExperiment(lfqdata = lfq_pep_log)
pe = QFeatures::QFeatures(list(peptide = se), colData = colData(se))

# Center each sample column on its median (vanilla msqrob2 vignette step).
# Without this, aggregateFeatures(robustSummary) bails on the raw log2
# precursor matrix during per-protein rlm fits.
pe = QFeatures::normalize(pe, i = "peptide", method = "center.median",
                          name = "peptide_norm")

# Aggregate centered log2 precursors to protein (QFeatures default
# = MsCoreUtils::robustSummary, per-protein rlm rollup).
pe = QFeatures::aggregateFeatures(
  pe, i = "peptide_norm", fcol = "protein_Id", name = "protein"
)

# Normalize at PROTEIN level.
prot_mat = SummarizedExperiment::assay(pe[["protein"]])
if (apply_vsn) {
  # vsn expects raw scale; protein matrix is log2 after aggregation.
  prot_mat = vsn_normalize_matrix(2 ^ prot_mat)
} else if (apply_quantile) {
  prot_mat = quantile_normalize_log2_matrix(prot_mat)
} else if (apply_median) {
  prot_mat = median_normalize_log2_matrix(prot_mat)
}
# else: keep the log2 protein matrix as-is (V1_log2 baseline).
SummarizedExperiment::assay(pe[["protein"]]) = prot_mat

pre_s = toc(t_pre)
t_mod = tic()
# Drop intercept so both factor levels get their own coefficient -- this
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
