library(data.table)

spectronaut_metric = function(file, method, p_col, fc_col) {
  dt = fread(file)
  p = dt[[p_col]]
  fc = dt[[fc_col]]
  label = dt[["Label"]]
  finite_fc = is.finite(fc)
  significant = p < 0.05 & finite_fc
  n_discoveries = sum(significant, na.rm = TRUE)
  false_discoveries = sum(significant & label == "Negative", na.rm = TRUE)

  data.table(
    Method = method,
    CSF_TPR = sum(significant & label == "Positive", na.rm = TRUE) /
      sum(finite_fc & label == "Positive", na.rm = TRUE),
    CSF_PPV = if (n_discoveries > 0) {
      1 - false_discoveries / n_discoveries
    } else {
      NA_real_
    }
  )
}

method_order = c(
  "MSstats+",
  "MSstats+ (pre-swap)",
  "MSstats",
  "MSstats (pre-swap)",
  "limpa",
  "MaxLFQ + limma",
  "msqrob2",
  "DEqMS"
)

comparison_table = rbindlist(list(
  spectronaut_metric(
    file.path("MSstats+", "MSstats+_model.csv"),
    "MSstats+",
    "pvalue",
    "log2FC"
  ),
  spectronaut_metric(
    file.path("MSstats+_preswap", "MSstats+_preswap_model.csv"),
    "MSstats+ (pre-swap)",
    "pvalue",
    "log2FC"
  ),
  spectronaut_metric(
    file.path("MSstats", "MSstats_model.csv"),
    "MSstats",
    "pvalue",
    "log2FC"
  ),
  spectronaut_metric(
    file.path("MSstats_preswap", "MSstats_preswap_model.csv"),
    "MSstats (pre-swap)",
    "pvalue",
    "log2FC"
  ),
  spectronaut_metric(
    file.path("limpa", "limpa_model.csv"),
    "limpa",
    "pvalue",
    "logFC"
  ),
  spectronaut_metric(
    file.path("limma", "limma_model.csv"),
    "MaxLFQ + limma",
    "pvalue",
    "logFC"
  ),
  spectronaut_metric(
    file.path("msqrob2", "msqrob2_model.csv"),
    "msqrob2",
    "pval",
    "logFC"
  ),
  spectronaut_metric(
    file.path("DEqMS", "deqms_model.csv"),
    "DEqMS",
    "pvalue",
    "logFC"
  )
))

comparison_table[, Method := factor(Method, levels = method_order)]
setorder(comparison_table, Method)

rounded = copy(comparison_table)
rounded[, c("CSF_TPR", "CSF_PPV") := lapply(.SD, round, 3),
        .SDcols = c("CSF_TPR", "CSF_PPV")]

fwrite(rounded, "CSF_Spectronaut_comparison_table.csv")
writeLines(
  capture.output(print(rounded)),
  "CSF_Spectronaut_comparison_table.txt"
)
print(rounded)
