## Build the TPR/PPV comparison table for the precursor-swap benchmark.
## Reads model outputs from V1_log2/<METHOD>{,_preswap}/ and
## v2_vsn/<METHOD>{,_preswap}/ subdirs and emits a long-format table.
library(data.table)

spectronaut_metric = function(file, method, variant, swap_state, p_col, fc_col) {
  if (!file.exists(file)) {
    return(data.table(Method = method, Variant = variant, SwapState = swap_state,
                      TPR = NA_real_, PPV = NA_real_))
  }
  dt = fread(file)
  p     = dt[[p_col]]
  fc    = dt[[fc_col]]
  label = dt[["Label"]]
  finite_fc     = is.finite(fc)
  significant   = p < 0.05 & finite_fc
  n_discoveries = sum(significant, na.rm = TRUE)
  false_disc    = sum(significant & label == "Negative", na.rm = TRUE)

  data.table(
    Method    = method,
    Variant   = variant,
    SwapState = swap_state,
    TPR = sum(significant & label == "Positive", na.rm = TRUE) /
          sum(finite_fc & label == "Positive", na.rm = TRUE),
    PPV = if (n_discoveries > 0) 1 - false_disc / n_discoveries
          else NA_real_
  )
}

method_specs = list(
  # method label,           file name within method dir,    p_col,    fc_col
  c("MSstats+",            "MSstats+_model.csv",            "pvalue", "log2FC"),
  c("MSstats",             "MSstats_model.csv",             "pvalue", "log2FC"),
  c("limpa",               "limpa_model.csv",               "pvalue", "logFC"),
  c("MaxLFQ + limma",      "limma_model.csv",               "pvalue", "logFC"),
  c("msqrob2",             "msqrob2_model.csv",             "pval",   "logFC"),
  c("DEqMS",               "deqms_model.csv",               "pvalue", "logFC"),
  c("prolfqua",            "prolfqua_model.csv",            "pvalue", "logFC")
)

method_dirs = c(
  "MSstats+"       = "MSstats+",
  "MSstats"        = "MSstats",
  "limpa"          = "limpa",
  "MaxLFQ + limma" = "limma",
  "msqrob2"        = "msqrob2",
  "DEqMS"          = "DEqMS",
  "prolfqua"       = "prolfqua"
)

rows = list()
for (variant in c("V1_log2", "v2_vsn")) {
  for (swap_state in c("post", "pre")) {
    suffix = if (swap_state == "pre") "_preswap" else ""
    for (s in method_specs) {
      method   = s[[1]]
      file_nm  = s[[2]]
      p_col    = s[[3]]
      fc_col   = s[[4]]
      dir = method_dirs[[method]]
      file = file.path(variant, paste0(dir, suffix), file_nm)
      # MSstats+ / MSstats only computed for V1_log2
      if (variant == "v2_vsn" && method %in% c("MSstats+", "MSstats")) next
      rows[[length(rows) + 1]] = spectronaut_metric(file, method, variant,
                                                     swap_state, p_col, fc_col)
    }
  }
}
comparison = rbindlist(rows)

method_order = c("MSstats+", "MSstats", "limpa", "MaxLFQ + limma",
                 "msqrob2", "DEqMS", "prolfqua")
comparison[, Method  := factor(Method,  levels = method_order)]
comparison[, Variant := factor(Variant, levels = c("V1_log2", "v2_vsn"))]
comparison[, SwapState := factor(SwapState, levels = c("post", "pre"))]
setorder(comparison, SwapState, Variant, Method)

rounded = copy(comparison)
rounded[, c("TPR", "PPV") := lapply(.SD, round, 3),
        .SDcols = c("TPR", "PPV")]

fwrite(rounded, "CSF_Spectronaut_swap_comparison_table.csv")
writeLines(capture.output(print(rounded)),
            "CSF_Spectronaut_swap_comparison_table.txt")
print(rounded)
