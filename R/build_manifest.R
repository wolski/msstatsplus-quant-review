## build_manifest.R - emit quant/manifest.csv enumerating every cell in the
## §0.2 grid: <csffolder, dataset, normalization, package>.
## Run with: Rscript R/build_manifest.R

suppressPackageStartupMessages(library(data.table))

datasets_per_folder <- list(
  CSF_Spectronaut                = c("all_data", "good_data"),
  CSF_Spectronaut_protein_swap   = c("all_data", "good_data", "small", "small_good_data"),
  CSF_Spectronaut_sample_swap    = c("all_data", "good_data", "small_good_data"),
  Mix_of_Proteome                = c("all_data")
)
normalizations <- c("log2", "median", "quantile")
packages <- c("MSstats", "MSstats+", "MaxLFQ_limma", "DEqMS",
              "msqrob2", "prolfqua", "limpa")

rows <- list()
for (csf in names(datasets_per_folder)) {
  for (dset in datasets_per_folder[[csf]]) {
    for (norm in normalizations) {
      for (pkg in packages) {
        rows[[length(rows) + 1L]] <- data.table(
          csffolder     = csf,
          dataset       = dset,
          normalization = norm,
          package       = pkg
        )
      }
    }
  }
}
manifest <- rbindlist(rows)
fwrite(manifest, file = "manifest.csv")
cat(sprintf("[manifest] wrote %d rows to manifest.csv\n", nrow(manifest)))
