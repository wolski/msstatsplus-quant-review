## Standalone prolfqua runner — produces <variant>/prolfqua<suffix>/prolfqua_model.csv.
suppressPackageStartupMessages({
  library(prolfqua)
  library(vsn)
})
source("run_step_common.R")
source("run_prolfqua_step.R")

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
message("prolfqua_only finished")
