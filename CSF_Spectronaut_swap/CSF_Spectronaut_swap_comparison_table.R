## Build the TPR/PPV comparison table for the precursor-swap benchmark.
## Reads model outputs from <variant><OUT_TAG>/<METHOD>{,_preswap}/ subdirs.
## Env vars:
##   OUT_TAG : suffix appended to each variant directory (default "" — the
##             canonical V1_log2/ and v2_vsn/ outputs). Use the same OUT_TAG
##             that was set when running the processing scripts.
## Output filenames also carry OUT_TAG to avoid overwriting the canonical
## table.
library(data.table)

out_tag = Sys.getenv("OUT_TAG", unset = "all_dilutions")

read_timing = function(timing_file) {
  if (is.null(timing_file) || !file.exists(timing_file)) {
    return(list(preprocess_seconds = NA_real_, model_seconds = NA_real_))
  }
  t = fread(timing_file)
  if (nrow(t) == 0) {
    return(list(preprocess_seconds = NA_real_, model_seconds = NA_real_))
  }
  # New schema: separate preprocess and model columns.
  if (all(c("preprocess_seconds", "model_seconds") %in% colnames(t))) {
    return(list(
      preprocess_seconds = as.numeric(t$preprocess_seconds[1]),
      model_seconds      = as.numeric(t$model_seconds[1])
    ))
  }
  # Legacy schema: single 'seconds' column (treated as model time).
  if ("seconds" %in% colnames(t)) {
    return(list(preprocess_seconds = NA_real_,
                model_seconds = as.numeric(t$seconds[1])))
  }
  list(preprocess_seconds = NA_real_, model_seconds = NA_real_)
}

tpr_ppv = function(p, fc, label, threshold = 0.05) {
  finite_fc     = is.finite(fc)
  significant   = !is.na(p) & p < threshold & finite_fc
  n_disc        = sum(significant, na.rm = TRUE)
  false_disc    = sum(significant & label == "Negative", na.rm = TRUE)
  list(
    TPR = sum(significant & label == "Positive", na.rm = TRUE) /
          sum(finite_fc & label == "Positive", na.rm = TRUE),
    PPV = if (n_disc > 0) 1 - false_disc / n_disc else NA_real_
  )
}

spectronaut_metric = function(file, method, variant, swap_state, p_col, fc_col,
                              timing_file = NULL) {
  timing = read_timing(timing_file)
  if (!file.exists(file)) {
    return(data.table(Method = method, Variant = variant, SwapState = swap_state,
                      TPR_pval = NA_real_, PPV_pval = NA_real_,
                      TPR_FDR  = NA_real_, PPV_FDR  = NA_real_,
                      preprocess_seconds = timing$preprocess_seconds,
                      model_seconds      = timing$model_seconds))
  }
  dt = fread(file)
  fc    = dt[[fc_col]]
  label = dt[["Label"]]
  # raw p threshold at 0.05
  m_pval = tpr_ppv(dt[[p_col]], fc, label, threshold = 0.05)
  # FDR-adjusted threshold at 0.05; method already emits adj.pvalue.
  adj_col = if ("adj.pvalue" %in% colnames(dt)) "adj.pvalue"
            else if ("adj.P.Val" %in% colnames(dt)) "adj.P.Val"
            else if ("FDR" %in% colnames(dt)) "FDR"
            else NA_character_
  m_fdr  = if (!is.na(adj_col)) tpr_ppv(dt[[adj_col]], fc, label, threshold = 0.05)
           else list(TPR = NA_real_, PPV = NA_real_)

  data.table(
    Method    = method,
    Variant   = variant,
    SwapState = swap_state,
    TPR_pval = m_pval$TPR, PPV_pval = m_pval$PPV,
    TPR_FDR  = m_fdr$TPR,  PPV_FDR  = m_fdr$PPV,
    preprocess_seconds = timing$preprocess_seconds,
    model_seconds      = timing$model_seconds
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
for (variant in c("V1_log2", "v2_median", "v2_vsn", "v3_quantile")) {
  for (swap_state in c("post", "pre")) {
    suffix = if (swap_state == "pre") "_preswap" else ""
    for (s in method_specs) {
      method   = s[[1]]
      file_nm  = s[[2]]
      p_col    = s[[3]]
      fc_col   = s[[4]]
      dir = method_dirs[[method]]
      method_dir = file.path(out_tag, variant, paste0(dir, suffix))
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
comparison[, Variant := factor(Variant, levels = c("V1_log2", "v2_median", "v2_vsn", "v3_quantile"))]
comparison[, SwapState := factor(SwapState, levels = c("post", "pre"))]
setorder(comparison, SwapState, Variant, Method)

rounded = copy(comparison)
num_cols = c("TPR_pval", "PPV_pval", "TPR_FDR", "PPV_FDR")
rounded[, (num_cols) := lapply(.SD, round, 3), .SDcols = num_cols]
rounded[, preprocess_seconds := round(preprocess_seconds, 1)]
rounded[, model_seconds      := round(model_seconds, 1)]

# Cells where MSstats methods are *not meaningfully different from v2_median*
# or where they fail outright. We mark them so the table doesn't look like
# MSstats was run under vsn or quantile -- it wasn't.
rounded[, Note := NA_character_]
ms_methods = c("MSstats+", "MSstats")
rounded[Variant == "v2_vsn" & Method %in% ms_methods,
        `:=`(TPR_pval = NA_real_, PPV_pval = NA_real_,
             TPR_FDR = NA_real_,  PPV_FDR  = NA_real_,
             preprocess_seconds = NA_real_, model_seconds = NA_real_,
             Note = "see v2_median")]
rounded[Variant == "v3_quantile" & Method %in% ms_methods,
        `:=`(TPR_pval = NA_real_, PPV_pval = NA_real_,
             TPR_FDR = NA_real_,  PPV_FDR  = NA_real_,
             preprocess_seconds = NA_real_, model_seconds = NA_real_,
             Note = "fail")]

dir.create(out_tag, recursive = TRUE, showWarnings = FALSE)
out_csv = file.path(out_tag, "comparison_table.csv")
out_txt = file.path(out_tag, "comparison_table.txt")
fwrite(rounded, out_csv)
writeLines(capture.output(print(rounded)), out_txt)
cat(sprintf("[table] wrote %s and %s\n", out_csv, out_txt))

## Wide-format markdown summary -----------------------------------------------
# Same column structure as manuscript Table 1: one TPR and one PPV column per
# variant. Two tables per swap state — one thresholded on raw p, one on FDR
# (adj.pvalue). Both at threshold 0.05.
fmt_num = function(x) {
  if (is.na(x)) NA_character_ else sprintf("%.3f", x)
}
write_md_table = function(rows, swap_state, fh, threshold = c("pval", "FDR")) {
  threshold = match.arg(threshold)
  tpr_col = paste0("TPR_", threshold)
  ppv_col = paste0("PPV_", threshold)
  label   = if (threshold == "pval") "raw p < 0.05" else "FDR (adj.pvalue) < 0.05"
  cat(sprintf("\n### %s-swap, %s\n\n", swap_state, label), file = fh)
  var_levels = levels(rows$Variant)
  wide = dcast(rows[SwapState == swap_state], Method ~ Variant,
               value.var = c(tpr_col, ppv_col, "Note"))
  setorder(wide, Method)
  header = c("Method", as.vector(rbind(paste0("TPR_", var_levels),
                                        paste0("PPV_", var_levels))))
  cat("| ", paste(header, collapse = " | "), " |\n", sep = "", file = fh)
  cat("|", paste(rep("---", length(header)), collapse = "|"), "|\n",
      sep = "", file = fh)
  for (i in seq_len(nrow(wide))) {
    row = wide[i]
    cells = character(0)
    for (v in var_levels) {
      note = row[[paste0("Note_", v)]]
      tpr  = row[[paste0(tpr_col, "_", v)]]
      ppv  = row[[paste0(ppv_col, "_", v)]]
      if (!is.na(note)) {
        cells = c(cells, note, note)
      } else if (is.na(tpr) && is.na(ppv)) {
        cells = c(cells, "—", "—")
      } else {
        cells = c(cells, fmt_num(tpr), fmt_num(ppv))
      }
    }
    cat("| ", as.character(row$Method), " | ",
        paste(cells, collapse = " | "), " |\n",
        sep = "", file = fh)
  }
}
out_md = file.path(out_tag, "comparison_table.md")
fh = file(out_md, open = "w")
cat(sprintf("# %s — TPR / PPV per method × variant\n\n",
            out_tag), file = fh)
cat("Cells are TPR or PPV at p < 0.05. Column suffix is the data ",
    "transformation/normalization variant.\n",
    "`fail` = method errored out;  ",
    "`see v2_median` = same MSstats internal normalization as v2_median;  ",
    "`—` = method not run for that variant.\n",
    file = fh, sep = "")
for (state in levels(rounded$SwapState)) {
  write_md_table(rounded, state, fh, threshold = "pval")
  write_md_table(rounded, state, fh, threshold = "FDR")
}
close(fh)
cat(sprintf("[table] wrote %s\n", out_md))
print(rounded)
