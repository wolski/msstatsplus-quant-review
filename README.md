## MSstats+ benchmark reanalysis

This repository contains a review and reanalysis of the quantitative proteomics
benchmark distributed with the MSstats+ manuscript and MassIVE reanalysis
`RMSV000000701.3`.

Rendered review vignettes: <https://wolski.github.io/msstatsplus-quant-review/>.

The repository is intended to make the review workflow reproducible: it records
the input layout expected from the MassIVE archive, the scripts used to build
benchmark subsets and synthetic swap ground truth, and the Makefile targets used
to run the method comparisons and render the review vignette.

### Related paper and data archive

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

### Repository structure

- `Makefile` is the top-level orchestration file. It includes all folder
  Makefiles and exposes aggregate targets (`prep`, `cells`, `review`,
  and the separately-invoked `diagnostics`, `render-vignettes`, `gh-pages`).
  Run `make help` for the full list.
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
  Currently excluded from the active pipeline pending DEqMS/msqrob2 tweaks
  for its 3 vs 3 design; see the comment block at the top of `Makefile`.
- `R/` contains the shared modelling code for MSstats, MSstats+, MaxLFQ+limma,
  DEqMS, msqrob2, prolfqua, limpa, normalization, timing, plotting, and
  comparison-table helpers.
- `src/` contains Python utilities for generating synthetic swapped reports
  and subset `Report.tsv`/`annotation.csv` directories, plus the bash scripts
  `render_vignettes.sh` and `publish_gh_pages.sh` that drive the GitHub
  Pages workflow.
- `vignettes/` contains the Quarto review, diagnostics, and `index.qmd`
  landing page for the published site.
- `results/` contains narrative notes and run summaries used during review.
- `TODO/` contains working notes and discrepancy reports; it is not part of
  the executable pipeline. `TODO/TODO_ghpages.md` documents the publishing
  pipeline design.

### Makefile workflow

Run commands from the repository root, i.e. from `quant/`. Run `make help`
for the categorized list of available targets (pipeline, per-norm shortcuts,
GitHub Pages, clean).

The Makefile is deliberately non-recursive: the top-level `Makefile` includes
the folder Makefiles, so `make -jN` can schedule independent cell blocks across
folders. Each cell-block internally uses ~8 cores via MSstats, so practical
choices are `-j 2` on a laptop (~16 cores) and `-j 4` on a workstation
(~32 cores).

Raw Spectronaut report files from the archive are not removed by any of the
`clean-*` targets.

#### Full rebuild from scratch

The `make all` step takes hours, so detach it from the terminal and
capture everything to `make.log`:

```bash
make clean-prep

# Long-running step: detached from terminal, stdout + stderr captured.
nohup make -j5 -k all > make.log 2>&1 &

# Watch progress live; Ctrl-C to detach (the build keeps running).
tail -f make.log

# Once `make all` has finished:
make diagnostics
make review-best-effort
```

`clean-prep` removes every generated artefact (subsets, model outputs,
synthetic swap reports, truth files, prep stamps) while keeping the raw
Spectronaut symlinks intact. `-k` keeps the build going past known-failing
cells (MSstats+ × quantile on some subsets, DEqMS × small × log2); without
`-k` the strict `review` step would block on the first missing stamp.

`nohup` keeps the job alive after logout / SSH disconnect; `> make.log`
redirects stdout; `2>&1` redirects stderr to the same file; `&` runs in
the background. `make.log` is git-ignored.

`make diagnostics` is intentionally not part of `all` and runs sequentially
(forces `-j1` internally), because every diagnostics target renders the
same `vignettes/diagnostics.qmd` and concurrent Quarto invocations race on
the shared `vignettes/.quarto/` scratch directory.

#### Checking on a backgrounded run

```bash
jobs                                  # if you stayed in the same shell
ps -fp $(pgrep -f 'make -j')          # find the make process(es)
tail -n 50 make.log                   # last 50 log lines (one-shot)
tail -f make.log                      # follow live
```

#### Publishing rendered vignettes to GitHub Pages

The site at <https://wolski.github.io/msstatsplus-quant-review/> is served
from the `gh-pages` branch of this repo. Rendering and publishing are
deliberately separate so the publish step never touches `main` and can be
re-run without re-rendering. From `quant/`:

```bash
make render-vignettes   # quarto-render the vignettes; outputs *.html in vignettes/
make gh-pages           # force-push the staged HTMLs to origin/gh-pages
```

`render-vignettes` is best-effort: a single failed `.qmd` does not abort the
others. `gh-pages` then verifies that every HTML it needs to publish exists
on disk and aborts with a clear error if any are missing — it does not
silently re-render. After the first push, enable Pages once via the repo
Settings: *Pages → Build and deployment → Source: Deploy from a branch →
Branch: `gh-pages` / `/ (root)`*.

See `TODO/TODO_ghpages.md` for the pipeline design rationale.

### Reproducibility notes

The Python helper environment is described by `pyproject.toml`; it requires
Python 3.10 or newer and `polars`. The R workflow expects the packages used by
the model adapters in `R/`, including MSstats, MSstatsConvert, limma, DEqMS,
msqrob2, prolfqua, limpa, data.table, and Quarto for report rendering.

Large raw reports, generated subset reports, model outputs, logs, and local prep
stamps are generated artifacts. They should be restored from MassIVE or produced
with `make prep`/`make cells`, not committed as source files.
