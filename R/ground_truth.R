## ground_truth.R - TP/TN constructors and the label_proteins consumer.
## Three benchmark types: protein swap, sample swap, species mix.

truth_protein_swap <- function(swap_list_path) {
  d <- data.table::fread(swap_list_path)
  stopifnot(all(c("Protein", "Label") %in% names(d)))
  d <- d[, .(Protein, Label)]
  d$Label <- ifelse(d$Label == "Positive", "Positive", "Negative")
  as.data.frame(d)
}

truth_sample_swap <- function(ground_truth_tsv) {
  d <- data.table::fread(ground_truth_tsv)
  stopifnot(all(c("Protein", "Label") %in% names(d)))
  as.data.frame(d[, .(Protein, Label)])
}

truth_mix_of_proteome <- function(idmapping_path,
                                   negative_organism = "Homo sapiens") {
  d <- data.table::fread(idmapping_path)
  protein_col <- intersect(c("Entry", "Protein", "ProteinId"), names(d))[1]
  organism_col <- intersect(c("Organism", "OS"), names(d))[1]
  stopifnot(!is.na(protein_col), !is.na(organism_col))
  out <- data.frame(
    Protein = d[[protein_col]],
    Label = ifelse(grepl(negative_organism, d[[organism_col]], fixed = TRUE),
                   "Negative", "Positive"),
    stringsAsFactors = FALSE
  )
  out
}

## label_proteins: left-join Label onto a per-package model output.
## model_df must have a `Protein` column. Unmatched proteins get Label = NA;
## downstream code filters by Label %in% c("Positive","Negative").
label_proteins <- function(model_df, truth_df) {
  stopifnot("Protein" %in% names(model_df), all(c("Protein", "Label") %in% names(truth_df)))
  # Authors' scripts (CSF_Spectronaut/run_msstats.R etc.) already write a
  # Label column into their model CSVs. Drop it before merge so the result
  # has a single canonical Label (not Label.x/Label.y suffixed by merge).
  if ("Label" %in% names(model_df)) model_df$Label <- NULL
  merge(model_df, truth_df[, c("Protein", "Label")], by = "Protein",
        all.x = TRUE, sort = FALSE)
}
