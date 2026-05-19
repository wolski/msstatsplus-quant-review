## build_sample_swap_list.R - regenerate the CSF sample-swap ground-truth list
## on the full 3,041-protein universe (the one used by CSF_Spectronaut_protein_swap),
## with a 10 % Positive / 90 % Negative stratified split.
##
## Why: the authors' sample-swap benchmark uses a narrower 1,820-protein list
## (Spectronaut-quantified subset). This script reproduces the same 10/90 design
## on the full 3,041-protein universe so the sample_swap and protein_swap
## benchmarks operate on the same protein set.
##
## Input:  CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv (3,041 rows;
##         the Protein column is used, the Label column is replaced).
## Output: CSF_Spectronaut_sample_swap/CSF_protein_swap_list.csv (Protein, Label).
##
## Usage:  Rscript R/build_sample_swap_list.R

suppressPackageStartupMessages(library(data.table))

source_path <- "CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv"
out_dir     <- "CSF_Spectronaut_sample_swap"
out_path    <- file.path(out_dir, "CSF_protein_swap_list.csv")

src <- fread(source_path)
stopifnot("Protein" %in% names(src))
proteins <- unique(src$Protein)
n <- length(proteins)

set.seed(123)
n_pos <- round(0.10 * n)
pos_idx <- sort(sample.int(n, n_pos))
labels  <- rep("Negative", n)
labels[pos_idx] <- "Positive"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out <- data.table(Protein = proteins, Label = labels)
fwrite(out, out_path)

cat(sprintf("[build_sample_swap_list] wrote %s: %d proteins (%d Positive / %d Negative)\n",
            out_path, nrow(out), sum(out$Label == "Positive"),
            sum(out$Label == "Negative")))
