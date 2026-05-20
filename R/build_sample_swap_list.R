## build_sample_swap_list.R - build the sample-swap ground truth.
##
## Emulates the authors' CSF sample-swap design: 10 % of "quantified
## proteins" keep their labels (Positive), 90 % get their condition
## labels permuted (Negative). Protein universe = proteins quantified
## in non-blank runs with at least 2 unique precursors (matching the
## protein-swap script's filter in src/swap_spectronaut_report.py
## compute_protein_stats — gives both benchmarks the same 2244-protein
## universe and lets sample-swap and protein-swap PPV/TPR be compared
## directly).
##
## Input:  CSF_Spectronaut/Report.tsv  (Spectronaut DIA report)
## Output: CSF_Spectronaut_sample_swap/CSF_protein_swap_list.csv
##           columns: Protein, Label  (Label in {"Positive","Negative"})
##
## Usage:  Rscript R/build_sample_swap_list.R
##
## No dependency on the protein-swap folder.

suppressPackageStartupMessages(library(data.table))

report_path <- "CSF_Spectronaut/Report.tsv"
out_dir     <- "CSF_Spectronaut_sample_swap"
out_path    <- file.path(out_dir, "CSF_protein_swap_list.csv")

if (!file.exists(report_path)) {
  stop("Not found: ", report_path, " (run `make symlinks` first)")
}

cat(sprintf("[build_sample_swap_list] reading %s ...\n", report_path))
report <- fread(report_path, sep = "\t",
                select = c("R.Condition", "PG.ProteinGroups",
                           "EG.PrecursorId", "F.PeakArea"))
report <- report[tolower(R.Condition) != "blank"]
report <- report[is.finite(F.PeakArea) & F.PeakArea > 0]
report <- report[nzchar(PG.ProteinGroups)]

# Per-protein precursor count, then filter n_precursors >= 2 (same as the
# protein-swap script's compute_protein_stats filter).
prot_stats <- unique(report[, .(PG.ProteinGroups, EG.PrecursorId)])[
  , .(n_precursors = .N), by = PG.ProteinGroups
]
prot_stats <- prot_stats[n_precursors >= 2L]

proteins <- sort(prot_stats$PG.ProteinGroups)
n <- length(proteins)
cat(sprintf(
  "[build_sample_swap_list] %d proteins with n_precursors >= 2 in non-blank runs\n",
  n))

set.seed(123)
n_pos <- round(0.10 * n)
pos_idx <- sort(sample.int(n, n_pos))
labels  <- rep("Negative", n)
labels[pos_idx] <- "Positive"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out <- data.table(Protein = proteins, Label = labels)
fwrite(out, out_path)

cat(sprintf("[build_sample_swap_list] wrote %s: %d proteins (%d Positive / %d Negative)\n",
            out_path, nrow(out),
            sum(out$Label == "Positive"),
            sum(out$Label == "Negative")))
