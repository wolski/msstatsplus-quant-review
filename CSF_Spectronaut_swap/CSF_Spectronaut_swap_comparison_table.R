## Build the TPR/PPV comparison table for the precursor-swap benchmark.
## Reads model outputs from <variant><OUT_TAG>/<METHOD>{,_preswap}/ subdirs.
## Env vars:
##   OUT_TAG : suffix appended to each variant directory (default "" — the
##             canonical V1_log2/ and v2_vsn/ outputs). Use the same OUT_TAG
##             that was set when running the processing scripts.
## Output filenames also carry OUT_TAG to avoid overwriting the canonical
## table.
library(data.table)

out_tag = Sys.getenv("OUT_TAG", unset = "")

read_timing = function(timing_file) {
  if (!file.exists(timing_file)) return(NA_real_)
  t = fread(timing_file)
  if (!"seconds" %in% colnames(t) || nrow(t) == 0) return(NA_real_)
  as.numeric(t$seconds[1])
}

spectronaut_metric = function(file, method, variant, swap_state, p_col, fc_col,
                              timing_file = NULL) {
  if (!file.exists(file)) {
    return(data.table(Method = method, Variant = variant, SwapState = swap_state,
                      TPR = NA_real_, PPV = NA_real_,
                      seconds = read_timing(timing_file %||% "")))
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
          else NA_real_,
    seconds = read_timing(timing_file %||% "")
  )
}

`%||%` = function(a, b) if (is.null(a) || (is.character(a) && nchar(a) == 0)) b else a

method_specs = list(
  # method label,           file name within method dir,    p_col,    fc_col
  c("MSstats+",            "MSstats+_model.csv",            "pvalue", "log2FC"),
  c("MSstats",             "MSstats_model.csv",             "pvalue", "log2FC"),
  c("limpa",               "limpa_model.csv",               "pvalue", "logFC"),
  c("MaxLFQ + limma",      "limma_model.csv",               "pvalue", "logFC"),
  c("msqrob2",             "msqrob2_model.csv",             "pvalue", "logFC"),
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
for (variant in c("V1_log2", "v2_vsn", "v3_quantile")) {
  for (swap_state in c("post", "pre")) {
    suffix = if (swap_state == "pre") "_preswap" else ""
    for (s in method_specs) {
      method   = s[[1]]
      file_nm  = s[[2]]
      p_col    = s[[3]]
      fc_col   = s[[4]]
      dir = method_dirs[[method]]
      method_dir = file.path(paste0(variant, out_tag), paste0(dir, suffix))
      file = file.path(method_dir, file_nm)
      timing_file = file.path(method_dir, paste0(dir, "_timing.csv"))
      rows[[length(rows) + 1]] = spectronaut_metric(file, method, variant,
                                                     swap_state, p_col, fc_col,
                                                     timing_file = timing_file)
    }
  }
}
comparison = rbindlist(rows)

method_order = c("MSstats+", "MSstats", "limpa", "MaxLFQ + limma",
                 "msqrob2", "DEqMS", "prolfqua")
comparison[, Method  := factor(Method,  levels = method_order)]
comparison[, Variant := factor(Variant, levels = c("V1_log2", "v2_vsn", "v3_quantile"))]
comparison[, SwapState := factor(SwapState, levels = c("post", "pre"))]
setorder(comparison, SwapState, Variant, Method)

rounded = copy(comparison)
rounded[, c("TPR", "PPV") := lapply(.SD, round, 3),
        .SDcols = c("TPR", "PPV")]
rounded[, seconds := round(seconds, 1)]

out_csv = paste0("CSF_Spectronaut_swap_comparison_table", out_tag, ".csv")
out_txt = paste0("CSF_Spectronaut_swap_comparison_table", out_tag, ".txt")
fwrite(rounded, out_csv)
writeLines(capture.output(print(rounded)), out_txt)
cat(sprintf("[table] wrote %s and %s\n", out_csv, out_txt))
print(rounded)
