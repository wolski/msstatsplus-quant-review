# MSstats+ benchmark reanalysis

This repository contains a review and reanalysis of the quantitative proteomics
benchmark distributed with the MSstats+ manuscript and MassIVE reanalysis
`RMSV000000701.3`.

The repository is intended to make the review workflow reproducible: it records
the input layout expected from the MassIVE archive, the scripts used to build
benchmark subsets and synthetic swap ground truth, and the Makefile targets used
to run the method comparisons and render the review vignette.

## Related paper and data archive

- bioRxiv preprint: [Accounting for longitudinal peak quality metrics with
  MSstats+ enhances differential analysis in proteomic experiments with
  data-independent acquisition](https://doi.org/10.1101/2025.09.11.675573),
  Devon Kohler, Eralp Dogu, Mrittika Bhattacharya, Ozge Karayel, Manuel Magana,
  Anthony Wu, Veronica G. Anania, and Olga Vitek.
- ProteomeXchange dataset: [PXD066486](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD066486).
- MassIVE dataset: [MSV000098622](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=fdf984642108451cbdece7303e47f2d2).
- MassIVE FTP archive: `ftp://massive-ftp.ucsd.edu/v10/MSV000098622/`.
- MassIVE reanalysis page verified during this review:
  [RMSV000000701.2](https://massive.ucsd.edu/ProteoSAFe/QueryMSV?id=RMSV000000701.2).

The local working directory is named `RMSV000000701.3-rerun/quant` because it
tracks the third received/reviewed revision of the quantification reanalysis.
At the time this README was written, I could verify the public
`RMSV000000701.2` MassIVE page, but not a public `RMSV000000701.3` URL.

## Repository structure

- `Makefile` is the top-level orchestration file. It includes all folder
  Makefiles and exposes aggregate targets such as `prep`, `cells`,
  `diagnostics`, and `review`.
- `mk/common.mk` defines shared variables, shared R helper dependencies, and
  the CSF raw-input symlink rules used by multiple benchmark folders.
- `CSF_Spectronaut/` reproduces the authors' CSF Spectronaut benchmark branch.
  It uses log2 normalization only and runs the original per-folder R scripts.
- `CSF_Spectronaut_protein_swap/` builds a synthetic protein-swap benchmark
  from the CSF Spectronaut report, then evaluates all configured methods across
  subsets and normalization choices.
- `CSF_Spectronaut_sample_swap/` builds a sample-swap benchmark directly from
  the CSF Spectronaut report and evaluates the same method grid.
- `Mix_of_Proteome/` runs the controlled mixture-of-proteomes benchmark.
- `R/` contains the shared modelling code for MSstats, MSstats+, MaxLFQ+limma,
  DEqMS, msqrob2, prolfqua, limpa, normalization, timing, plotting, and
  comparison-table helpers.
- `src/` contains Python utilities for generating synthetic swapped reports and
  subset `Report.tsv`/`annotation.csv` directories.
- `vignettes/` contains the Quarto review and diagnostics reports.
- `results/` contains narrative notes and run summaries used during review.
- `TODO/` contains working notes and discrepancy reports; it is not part of the
  executable pipeline.

## Makefile workflow

Run commands from the repository root, i.e. from `quant/`.

```bash
make symlinks
make prep
make cells
make diagnostics
make review
```

The main targets are:

- `make symlinks`: creates local symlinks from canonical short names such as
  `Report.tsv` and `annotation.csv` to the raw archive filenames.
- `make prep`: creates synthetic swap reports and subset directories. This is
  where `all_data`, `good_data`, `small`, and `small_good_data` inputs are
  generated.
- `make cells`: runs model-fitting cells. Each cell corresponds to a benchmark
  folder, subset, normalization, and method block.
- `make diagnostics`: renders per-cell diagnostics once cell outputs exist.
- `make review`: strict review render. This depends on all expected cell stamps.
- `make review-best-effort`: renders the review against whatever cell outputs
  are present, but still requires the prep outputs used as truth tables.

The Makefile is deliberately non-recursive: the top-level `Makefile` includes
the folder Makefiles, so `make -jN` can schedule independent cell blocks across
folders. The comments at the top of `Makefile` describe practical parallelism
settings. The clean targets are layered:

- `make clean-models` removes model outputs and diagnostics.
- `make clean-subsets` also removes generated subset directories and prep
  stamps.
- `make clean-prep` also removes generated swap reports and root-level swap
  truth files.

Raw Spectronaut report files from the archive are not removed by the clean
targets.

## Reproducibility notes

The Python helper environment is described by `pyproject.toml`; it requires
Python 3.10 or newer and `polars`. The R workflow expects the packages used by
the model adapters in `R/`, including MSstats, MSstatsConvert, limma, DEqMS,
msqrob2, prolfqua, limpa, data.table, and Quarto for report rendering.

Large raw reports, generated subset reports, model outputs, logs, and local prep
stamps are generated artifacts. They should be restored from MassIVE or produced
with `make prep`/`make cells`, not committed as source files.
