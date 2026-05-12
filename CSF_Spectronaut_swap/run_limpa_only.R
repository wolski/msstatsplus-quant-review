## Standalone limpa runner — produces <variant>/limpa<suffix>/limpa_model.csv.
suppressPackageStartupMessages({
  library(limma)
  library(limpa)
  library(vsn)
})
source("run_step_common.R")

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
fwrite(limpa_model, file = file.path(out_dir("limpa"), "limpa_model.csv"))
message("limpa_only finished")
