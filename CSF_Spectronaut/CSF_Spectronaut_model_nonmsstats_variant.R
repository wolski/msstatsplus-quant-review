args = commandArgs(trailingOnly = TRUE)

if (length(args) != 1 || !args[[1]] %in% c("log2", "vsn")) {
  stop("Usage: Rscript CSF_Spectronaut_model_nonmsstats_variant.R <log2|vsn>",
       call. = FALSE)
}

variant = args[[1]]

required_packages = c(
  "data.table",
  "tidyverse",
  "QFeatures",
  "msqrob2",
  "iq",
  "limma",
  "limpa",
  "DEqMS"
)

if (variant == "vsn") {
  required_packages = c(required_packages, "vsn")
}

missing_packages = required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

library(data.table)
library(tidyverse)
library(QFeatures)
library(msqrob2)
library(iq)
library(limma)
library(limpa)
library(DEqMS)

source("../benchmark_experiments_functions.R")

normalize_matrix = function(input_matrix, variant) {
  input_matrix = as.matrix(input_matrix)
  input_dimnames = dimnames(input_matrix)
  storage.mode(input_matrix) = "double"

  if (variant == "log2") {
    normalized_matrix = log2(input_matrix)
  } else {
    normalized_matrix = vsn::justvsn(input_matrix)
  }

  dimnames(normalized_matrix) = input_dimnames
  normalized_matrix
}

summarize_deqms_no_ref_col = function(dat, group_col = 2) {
  dat.ratio = dat
  dat.ratio[, 3:ncol(dat)] = dat.ratio[, 3:ncol(dat)] -
    matrixStats::rowMedians(as.matrix(dat.ratio[, 3:ncol(dat)]), na.rm = TRUE)
  dat.summary = plyr::ddply(
    dat.ratio,
    colnames(dat)[group_col],
    function(x) {
      matrixStats::colMedians(as.matrix(x[, 3:ncol(dat)]), na.rm = TRUE)
    }
  )
  colnames(dat.summary)[2:ncol(dat.summary)] = colnames(dat)[3:ncol(dat)]
  dat.new = dat.summary[, -1]
  rownames(dat.new) = dat.summary[, 1]
  dat.new
}

spectronaut_metric = function(file, method, source, p_col, fc_col) {
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
    Source = source,
    CSF_TPR = sum(significant & label == "Positive", na.rm = TRUE) /
      sum(finite_fc & label == "Positive", na.rm = TRUE),
    CSF_PPV = if (n_discoveries > 0) {
      1 - false_discoveries / n_discoveries
    } else {
      NA_real_
    }
  )
}

write_comparison_table = function(output_dir, variant_label) {
  method_order = c(
    "MSstats+",
    "MSstats",
    "limpa",
    "MaxLFQ + limma",
    "msqrob2",
    "DEqMS"
  )

  comparison_table = rbindlist(list(
    spectronaut_metric(
      file.path("MSstats+", "MSstats+_model.csv"),
      "MSstats+",
      "baseline unchanged",
      "pvalue",
      "log2FC"
    ),
    spectronaut_metric(
      file.path("MSstats", "MSstats_model.csv"),
      "MSstats",
      "baseline unchanged",
      "pvalue",
      "log2FC"
    ),
    spectronaut_metric(
      file.path(output_dir, "limpa_model.csv"),
      "limpa",
      variant_label,
      "pvalue",
      "logFC"
    ),
    spectronaut_metric(
      file.path(output_dir, "limma_model.csv"),
      "MaxLFQ + limma",
      variant_label,
      "pvalue",
      "logFC"
    ),
    spectronaut_metric(
      file.path(output_dir, "msqrob2_model.csv"),
      "msqrob2",
      variant_label,
      "pval",
      "logFC"
    ),
    spectronaut_metric(
      file.path(output_dir, "deqms_model.csv"),
      "DEqMS",
      variant_label,
      "pvalue",
      "logFC"
    )
  ))

  comparison_table[, Method := factor(Method, levels = method_order)]
  setorder(comparison_table, Method)

  rounded = copy(comparison_table)
  rounded[, c("CSF_TPR", "CSF_PPV") := lapply(.SD, round, 3),
          .SDcols = c("CSF_TPR", "CSF_PPV")]

  fwrite(rounded, file.path(output_dir, "CSF_Spectronaut_comparison_table.csv"))
  writeLines(
    capture.output(print(rounded)),
    file.path(output_dir, "CSF_Spectronaut_comparison_table.txt")
  )

  rounded
}

write_readme = function(output_dir, variant, variant_label) {
  readme_lines = c(
    paste0("# CSF Spectronaut non-MSstats variant: ", variant),
    "",
    paste0("Variant label: ", variant_label),
    "",
    "Input:",
    "- Reads raw, filtered, benchmark-swapped non-MSstats inputs from `../nonmsstats_preprocessed/` relative to this folder.",
    "",
    "Methods rerun:",
    "- msqrob2",
    "- MaxLFQ + limma",
    "- limpa",
    "- DEqMS",
    "",
    "Methods intentionally not rerun:",
    "- MSstats+ and base MSstats. Their scores are read from the baseline model CSVs and included as unchanged rows.",
    "- mapDIA. It is external and is not included in this variant table.",
    "",
    "Output files:",
    "- `msqrob2_model.csv`",
    "- `limma_model.csv`",
    "- `limpa_model.csv`",
    "- `deqms_model.csv`",
    "- `CSF_Spectronaut_comparison_table.csv`",
    "- `CSF_Spectronaut_comparison_table.txt`",
    "",
    "Dependency versions:",
    paste0("- R: ", getRversion()),
    paste0("- msqrob2: ", as.character(utils::packageVersion("msqrob2"))),
    paste0("- limma: ", as.character(utils::packageVersion("limma"))),
    paste0("- limpa: ", as.character(utils::packageVersion("limpa"))),
    paste0("- DEqMS: ", as.character(utils::packageVersion("DEqMS")))
  )

  if (variant == "vsn") {
    readme_lines = c(
      readme_lines,
      paste0("- vsn: ", as.character(utils::packageVersion("vsn")))
    )
  }

  writeLines(readme_lines, file.path(output_dir, "README.md"))
}

preprocess_dir = file.path(".", "nonmsstats_preprocessed")

required_inputs = file.path(
  preprocess_dir,
  c(
    "annotation.csv",
    "truth_labels.csv",
    "msqrob_feature_peakarea_wide.csv",
    "limma_feature_peakarea_long.csv",
    "feature_peakarea_wide.csv"
  )
)

if (!all(file.exists(required_inputs))) {
  stop(
    "Missing preprocessed input(s). Run CSF_Spectronaut_preprocess_nonmsstats.R first.",
    call. = FALSE
  )
}

output_dir = if (variant == "log2") {
  file.path(".", "V1_log2")
} else {
  file.path(".", "v2_vsn")
}

variant_label = if (variant == "log2") {
  "log2 rerun"
} else {
  "VSN normalized rerun"
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
unlink(
  file.path(
    output_dir,
    c(
      "msqrob_obj.rda",
      "msqrob2_model.csv",
      "limma_model.csv",
      "limpa_model.csv",
      "deqms_model.csv",
      "CSF_Spectronaut_comparison_table.csv",
      "CSF_Spectronaut_comparison_table.txt",
      "README.md"
    )
  ),
  force = TRUE
)

annotation = fread(file.path(preprocess_dir, "annotation.csv"))
truth_labels = fread(file.path(preprocess_dir, "truth_labels.csv"))
true_positives = truth_labels[Label == "Positive", Protein]

## MSqRob2 ---------------------------------------------------------------------
msqrob_wide = fread(file.path(preprocess_dir, "msqrob_feature_peakarea_wide.csv"))
colnames(msqrob_wide)[3:length(colnames(msqrob_wide))] = paste(
  "F.PeakArea",
  colnames(msqrob_wide)[3:length(colnames(msqrob_wide))],
  sep = "_"
)
ecols = grep("F.PeakArea", colnames(msqrob_wide))

pe = readQFeatures(
  msqrob_wide,
  fnames = 2,
  quantCols = ecols,
  name = "peptideRaw"
)

colData(pe)$condition = annotation$Condition[
  match(
    str_remove(colnames(msqrob_wide)[grep("F.PeakArea", colnames(msqrob_wide))],
               "F.PeakArea_"),
    annotation$R.FileName
  )
] %>%
  as.factor()

rowData(pe[["peptideRaw"]])$nNonZero = rowSums(
  assay(pe[["peptideRaw"]]) > 0,
  na.rm = TRUE
)
pe = zeroIsNA(pe, "peptideRaw")

peptide_assay = if (variant == "log2") {
  pe = logTransform(pe, base = 2, i = "peptideRaw", name = "peptideNorm")
  "peptideNorm"
} else {
  pe = logTransform(pe, base = 2, i = "peptideRaw", name = "peptideNorm")
  assay(pe[["peptideNorm"]]) = normalize_matrix(
    assay(pe[["peptideRaw"]]),
    variant = variant
  )
  "peptideNorm"
}

protein_filter = rowData(
  pe[[peptide_assay]]
)$ProteinName %in% smallestUniqueGroups(
  rowData(pe[[peptide_assay]])$ProteinName
)
pe = pe[protein_filter,]

pe = filterFeatures(pe, ~ nNonZero >= 2)

pe = QFeatures::impute(
  pe,
  i = peptide_assay,
  name = "peptideImp",
  method = "QRILC"
)

pe = aggregateFeatures(
  pe,
  i = "peptideImp",
  fcol = "ProteinName",
  name = "protein"
)

pe = msqrob(object = pe, i = "protein", formula = ~ condition)
L = makeContrast(
  "conditionCondition2=0",
  parameterNames = c("conditionCondition2", "conditionCondition1")
)
pe = hypothesisTest(object = pe, i = "protein", contrast = L)
save(pe, file = file.path(output_dir, "msqrob_obj.rda"))

msqrob2_model = rowData(pe[["protein"]])$`conditionCondition2`
msqrob2_model$Label = ifelse(
  rownames(msqrob2_model) %in% true_positives,
  "Positive",
  "Negative"
)
msqrob2_model$Protein = rownames(msqrob2_model)
fwrite(msqrob2_model, file = file.path(output_dir, "msqrob2_model.csv"))
print(paste("msqrob2", variant, "finished"))

## MaxLFQ + limma --------------------------------------------------------------
limma_input = fread(file.path(preprocess_dir, "limma_feature_peakarea_long.csv"))

maxlfq_input = preprocess(
  limma_input,
  primary_id = "PG.ProteinGroups",
  secondary_id = c("Feature"),
  sample_id = "R.FileName",
  intensity_col = "F.PeakArea",
  median_normalization = FALSE,
  log2_intensity_cutoff = 0,
  pdf_out = file.path(output_dir, "qc-plots.pdf"),
  pdf_width = 12,
  pdf_height = 8,
  intensity_col_sep = NULL,
  intensity_col_id = NULL,
  na_string = "0"
)

maxlfq_summarized = fast_MaxLFQ(maxlfq_input)
maxlfq_summarized = maxlfq_summarized$estimate

if (variant == "vsn") {
  maxlfq_summarized = normalize_matrix(maxlfq_summarized, variant = variant)
}

class = annotation$Condition[
  match(colnames(maxlfq_summarized), annotation$R.FileName)
] %>%
  as.factor()
design = model.matrix(~ 0 + class)

fit1 = lmFit(maxlfq_summarized, design = design)
cont = makeContrasts(classCondition2 - classCondition1, levels = design)
fit2 = contrasts.fit(fit1, contrasts = cont)
fit3 = eBayes(fit2)

limma_model = data.frame(
  "Protein" = rownames(fit3$coefficients),
  "logFC" = as.numeric(fit3$coefficients),
  "SE" = as.numeric(sqrt(fit3$s2.post) * fit3$stdev.unscaled),
  "DF" = fit3$df.total,
  "pvalue" = as.numeric(fit3$p.value)
)

rownames(limma_model) = NULL
limma_model$adj.pvalue = p.adjust(limma_model$pvalue, method = "BH")
limma_model$Label = ifelse(
  limma_model$Protein %in% true_positives,
  "Positive",
  "Negative"
)
fwrite(limma_model, file = file.path(output_dir, "limma_model.csv"))
print(paste("limma", variant, "finished"))

## Limpa -----------------------------------------------------------------------
feature_peakarea = fread(file.path(preprocess_dir, "feature_peakarea_wide.csv"))
mapper = as.data.frame(feature_peakarea[, .(PG.ProteinGroups, Feature)])
row.names(mapper) = mapper$Feature
limpa_dt = as.data.frame(feature_peakarea)
row.names(limpa_dt) = limpa_dt$Feature
limpa_dt = limpa_dt[, !colnames(limpa_dt) %in% c("PG.ProteinGroups", "Feature")]
limpa_dt = normalize_matrix(limpa_dt, variant = variant)

targets = as.data.frame(annotation)[c("R.FileName", "Condition")]
row.names(targets) = targets$R.FileName
limpa_elist = new(
  "EList",
  list(E = limpa_dt, genes = mapper, targets = targets)
)

pdf(file.path(output_dir, "limpa_dpc_plots.pdf"), width = 12, height = 8)
dpcfit = dpc(limpa_elist)
plotDPC(dpcfit)

y.protein = dpcQuant(
  limpa_elist,
  "PG.ProteinGroups",
  dpc = dpcfit
)

class = annotation$Condition[
  match(colnames(y.protein$E), annotation$R.FileName)
] %>%
  as.factor()
design = model.matrix(~ 0 + class)

fit = dpcDE(y.protein, design, plot = TRUE)
dev.off()

cont = makeContrasts(classCondition2 - classCondition1, levels = design)
fit = contrasts.fit(fit, contrasts = cont)
fit = eBayes(fit)

limpa_model = data.frame(
  "Protein" = rownames(fit$coefficients),
  "logFC" = as.numeric(fit$coefficients),
  "SE" = as.numeric(sqrt(fit$s2.post) * fit$stdev.unscaled),
  "DF" = fit$df.total,
  "pvalue" = as.numeric(fit$p.value)
)

rownames(limpa_model) = NULL
limpa_model$adj.pvalue = p.adjust(limpa_model$pvalue, method = "BH")
limpa_model$Label = ifelse(
  limpa_model$Protein %in% true_positives,
  "Positive",
  "Negative"
)
fwrite(limpa_model, file = file.path(output_dir, "limpa_model.csv"))
print(paste("limpa", variant, "finished"))

## DEqMS -----------------------------------------------------------------------
deqms_input = as.data.frame(feature_peakarea)
deqms_input[, 3:ncol(deqms_input)] = normalize_matrix(
  deqms_input[, 3:ncol(deqms_input)],
  variant = variant
)

deqms_summarized = summarize_deqms_no_ref_col(deqms_input, group_col = 1)

pep_count = deqms_input %>%
  group_by(PG.ProteinGroups) %>%
  summarise(count = n_distinct(Feature), .groups = "drop") %>%
  as.data.frame()

pep_count$count = pep_count$count + 1
pep_count = as.data.frame(pep_count)
rownames(pep_count) = pep_count$PG.ProteinGroups
pep_count$PG.ProteinGroups = NULL

class = annotation$Condition[
  match(colnames(deqms_summarized), annotation$R.FileName)
] %>%
  as.factor()
design = model.matrix(~ 0 + class)

fit1 = lmFit(deqms_summarized, design = design)
cont = makeContrasts(classCondition2 - classCondition1, levels = design)
fit2 = contrasts.fit(fit1, contrasts = cont)
fit3 = eBayes(fit2)

fit3$count = pep_count[rownames(fit3$coefficients), "count"]
fit4 = spectraCounteBayes(fit3)

deqms_model = data.frame(
  "Protein" = rownames(fit4$coefficients),
  "logFC" = as.numeric(fit4$coefficients),
  "SE" = as.numeric(sqrt(fit4$s2.post) * fit4$stdev.unscaled),
  "DF" = fit4$df.total,
  "pvalue" = as.numeric(fit4$p.value)
)
rownames(deqms_model) = NULL
deqms_model$adj.pvalue = p.adjust(deqms_model$pvalue, method = "BH")
deqms_model$Label = ifelse(
  deqms_model$Protein %in% true_positives,
  "Positive",
  "Negative"
)
fwrite(deqms_model, file = file.path(output_dir, "deqms_model.csv"))
print(paste("DEqMS", variant, "finished"))

write_readme(output_dir, variant, variant_label)
comparison_table = write_comparison_table(output_dir, variant_label)
print(comparison_table)
