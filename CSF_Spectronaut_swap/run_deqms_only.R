## Standalone DEqMS runner — produces <variant>/DEqMS<suffix>/deqms_model.csv.
suppressPackageStartupMessages({
  library(limma)
  library(DEqMS)
  library(vsn)
})
source("run_step_common.R")

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
  deqms_summarized = as.data.frame(vsn_normalize_matrix(deqms_summarized))
}
pep_count = deqms_input |>
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
message("deqms_only finished")
