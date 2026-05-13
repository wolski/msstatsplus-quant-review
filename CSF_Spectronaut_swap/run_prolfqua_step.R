## prolfqua step — single source of truth for the prolfqua DE block.
##
## Sourced by run_nonmsstats.R. Returns a list with
##   $model              (data.frame: Protein, logFC, SE, DF, pvalue, adj.pvalue)
##   $preprocess_seconds (wall-clock: prep + normalization + aggregation)
##   $model_seconds      (wall-clock: ContrastsLMImputeFacade + get_contrasts)
##
## Uses base-R indexing to avoid the prolfqua/dplyr `rename` and `select`
## namespace collisions that bit earlier iterations of this code.

run_prolfqua_step = function(merged_input, annotation, all_proteins, no_swap,
                              normalization = "none",
                              vsn_func = NULL,
                              quantile_func = NULL,
                              median_func = NULL) {
  stopifnot(normalization %in% c("none", "vsn", "quantile", "median"))
  t_pre = proc.time()[3]
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

  # prolfquapp canonical aggregation:
  #   precursors -> natural log -> medpolish -> exp -> protein-level (linear)
  # Then normalize at PROTEIN level (vsn or log2+robscale), matching
  # prolfquapp::transform_lfqdata. Aggregating in log space then undoing the
  # log gives a stable additive decomposition; normalization at the protein
  # scale lets vsn (or robscale) see one value per protein per run.
  tr_log = lfqdata$get_Transformer()$intensity_array(log)
  agg = tr_log$lfq$get_Aggregator("medpolish")
  agg$aggregate()
  lfq_protein_log = agg$lfq_agg
  tr_inv = lfq_protein_log$get_Transformer()$intensity_array(exp, force = TRUE)
  lfq_protein = tr_inv$lfq
  lfq_protein$is_transformed(FALSE)

  tr_norm = lfq_protein$get_Transformer()
  if (normalization == "vsn") {
    stopifnot(!is.null(vsn_func))
    tr_norm$intensity_matrix(.func = vsn_func)
  } else if (normalization == "quantile") {
    # log2 first, then quantile on the log2-scale protein matrix. The
    # transformer flags log2'd data as "transformed" and refuses to run
    # intensity_matrix on it without force = TRUE -- without that, the
    # quantile step silently no-ops.
    stopifnot(!is.null(quantile_func))
    tr_norm$log2()
    tr_norm$intensity_matrix(.func = quantile_func, force = TRUE)
  } else if (normalization == "median") {
    stopifnot(!is.null(median_func))
    tr_norm$log2()
    tr_norm$intensity_matrix(.func = median_func, force = TRUE)
  } else {
    # "none" -> prolfquapp default: log2 + robscale at protein level.
    tr_norm$log2()
    tr_norm$robscale()
  }
  lfq_protein = tr_norm$lfq

  preprocess_seconds = as.numeric(proc.time()[3] - t_pre)
  t_mod = proc.time()[3]

  contr_spec = c("Condition2_vs_Condition1" =
                   "group_Condition2 - group_Condition1")
  fa = prolfqua::ContrastsLMImputeFacade$new(lfq_protein, "~ group_",
                                              contr_spec)
  res = fa$get_contrasts()
  model_seconds = as.numeric(proc.time()[3] - t_mod)

  list(
    model = data.frame(
      Protein    = res$protein_Id,
      logFC      = res$diff,
      SE         = res$std.error,
      DF         = res$df,
      pvalue     = res$p.value,
      adj.pvalue = res$FDR,
      stringsAsFactors = FALSE
    ),
    preprocess_seconds = preprocess_seconds,
    model_seconds      = model_seconds
  )
}
