## Standalone prolfqua runner — produces <variant>/prolfqua<suffix>/prolfqua_model.csv.
suppressPackageStartupMessages({
  library(prolfqua)
  library(vsn)
})
source("run_step_common.R")
source("run_prolfqua_step.R")

prolfqua_model = run_prolfqua_step(merged_input, annotation, all_proteins,
                                    no_swap, apply_vsn,
                                    vsn_func = vsn_normalize_matrix)
prolfqua_model = label_proteins(prolfqua_model)
fwrite(prolfqua_model,
        file = file.path(out_dir("prolfqua"), "prolfqua_model.csv"))
message("prolfqua_only finished")
