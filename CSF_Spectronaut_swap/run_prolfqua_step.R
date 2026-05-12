## prolfqua step — standalone, defensive against dplyr namespace conflicts.
##
## Sourced by the main processing script; can also be run directly with
## VARIANT and OUT_SUFFIX env vars to produce a missing prolfqua_model.csv
## without re-running the rest of the pipeline.

run_prolfqua_step = function(merged_input, annotation, all_proteins, no_swap,
                              apply_vsn, vsn_func) {
  # Build the long-form precursor table without going through the dplyr pipe,
  # to avoid namespace collisions with prolfqua::rename / select.
  prolfqua_input = prepare_data_for_limma(merged_input, all_proteins, no_swap)
  prolfqua_input = as.data.frame(prolfqua_input)

  needed = c("PG.ProteinGroups", "EG.PrecursorId", "Feature", "R.FileName",
              "Condition", "F.PeakArea")
  missing_cols = setdiff(needed, colnames(prolfqua_input))
  if (length(missing_cols) > 0) {
    stop("prolfqua: missing columns in prepare_data_for_limma output: ",
         paste(missing_cols, collapse = ", "))
  }

  prolfqua_input = prolfqua_input[is.finite(prolfqua_input$F.PeakArea) &
                                    prolfqua_input$F.PeakArea > 0, , drop = FALSE]

  # Feature from prepare_data_for_limma is fragment-level
  # (paste(EG.PrecursorId, F.FrgIon, F.Charge)). The full hierarchy is
  # protein -> precursor -> fragment; model at protein (depth = 1).
  config = prolfqua::AnalysisConfiguration$new()
  config$file_name = "R.FileName"
  config$factors["group_"] = "Condition"
  config$hierarchy[["protein_Id"]]   = "PG.ProteinGroups"
  config$hierarchy[["precursor_Id"]] = "EG.PrecursorId"
  # F.FrgLossType (noloss / water loss / ammonia loss) distinguishes rows
  # that share the same (precursor, ion, charge); without it setup_analysis
  # collapses to the duplicate-summary table.
  config$hierarchy[["fragment_Id"]]  = c("Feature", "F.FrgLossType")
  config$hierarchy_depth = 1
  config$set_response("F.PeakArea")
  adata = prolfqua::setup_analysis(prolfqua_input, config)
  lfqdata = prolfqua::LFQData$new(adata, config)

  tr = lfqdata$get_Transformer()
  if (apply_vsn) {
    tr$intensity_matrix(.func = vsn_func)
  } else {
    tr$log2()
    tr$robscale()
  }
  lfqdata_trans = tr$lfq

  agg = lfqdata_trans$get_Aggregator("medpolish")
  agg$aggregate()
  lfq_protein = agg$lfq_agg

  contr_spec = c("Condition2_vs_Condition1" =
                   "group_Condition2 - group_Condition1")
  fa = prolfqua::ContrastsLMImputeFacade$new(lfq_protein, "~ group_",
                                              contr_spec)
  res = fa$get_contrasts()

  data.frame(
    Protein    = res$protein_Id,
    logFC      = res$diff,
    SE         = res$std.error,
    DF         = res$df,
    pvalue     = res$p.value,
    adj.pvalue = res$FDR,
    stringsAsFactors = FALSE
  )
}
