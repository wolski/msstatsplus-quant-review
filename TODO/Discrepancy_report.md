# Discrepancy report

This report records discrepancies found while reviewing `vignettes/review.qmd`
against `TODO/WhatWeStillNeed.md`, plus code and vignette issues observed in
the shared `quant` R workflow.

## Scope

Reviewed files included:

- `TODO/WhatWeStillNeed.md`
- `TODO/TODO_step_by_step.md`
- `vignettes/review.qmd`
- `vignettes/diagnostics.qmd`
- `vignettes/swap_visualization.qmd`
- Shared R workflow files under `R/`
- Sweep scripts `run_all.fish` and `run_folder_sweep.fish`

No files were edited during the review. This report is the first written
artifact created from those findings.

## Findings

### 1. Protein-swap ground-truth counts in `review.qmd` are wrong

`vignettes/review.qmd` describes the protein-swap ground truth as
`1,820 + 9 = 1,829` proteins and approximately 50 swap pairs.

Observed from the on-disk files:

- `CSF_Spectronaut_protein_swap/CSF_protein_swap_list.csv`: 3,041 proteins
- Positive labels: 224
- Negative labels: 2,817
- `*_swap_ground_truth.tsv`: 112 protein-swap pairs

This should be corrected because it affects the benchmark description and the
interpretation of positive/negative class sizes.

Relevant location:

- `vignettes/review.qmd`, around the protein-swap dataset section.

### 2. Swap visualizations do not match the requested subset scope

`WhatWeStillNeed.md` and `TODO_step_by_step.md` ask for heatmap and NA heatmap
visualizations on `all_data`, with density plots on `good_data`.

Current `review.qmd` uses `good_data` for:

- sample-swap heatmap
- sample-swap NA heatmap
- protein-swap heatmap
- protein-swap NA heatmap

There is also an internal contradiction in `review.qmd`: a comment says the
protein-swap heatmap and NA heatmap intentionally stay on `all_data`, but the
code paths use `good_data`.

This is either a code-path error or a narrative/TODO mismatch. It should be
resolved explicitly before relying on the figures.

### 3. Parameterized vignette defaults are not runnable as checked

`vignettes/review.qmd` rendered successfully from the quant root:

```sh
cd RMSV000000701.3-rerun/quant
quarto render vignettes/review.qmd --to html
```

The parameterized vignettes failed when rendered the same way:

```sh
quarto render vignettes/diagnostics.qmd --to html
quarto render vignettes/swap_visualization.qmd --to html
```

Both fail because `base_dir: "."` resolves relative to the vignette execution
context, so paths such as `R/paths.R` and `R/figures.R` are not found.

Likely fixes to consider:

- Change the default `base_dir` to `..`, or
- Add robust path resolution like `review.qmd` uses.

### 4. `run_all.fish` skip detection does not match package output filenames

`run_all.fish` builds the expected model filename from the package name. This
does not match several canonical output names, for example:

- `MaxLFQ_limma` writes `limma_model.csv`
- `DEqMS` writes `deqms_model.csv`
- `msqrob2` writes `msqrob2_model.csv`

As a result, existing model outputs may not be detected correctly, and reruns
may happen unexpectedly. `run_folder_sweep.fish` uses a broader
`*_model.csv` check and is less exposed to this specific issue.

### 5. `run_nonmsstats_block.R` can hide package failures from automation

`R/run_nonmsstats_block.R` wraps each non-MSstats package run in `tryCatch` and
continues after errors. This is useful for long sweeps, but the script can still
exit successfully when one or more requested packages failed.

Downstream automation must parse logs to notice these failures. For benchmark
reproducibility, consider making the script exit non-zero when any requested
package fails, while still running the remaining packages first.

## Alignment with `WhatWeStillNeed.md`

The current `review.qmd` broadly covers the requested narrative:

- Introduction about Vitek/MSstats and benchmark claims
- Table 1 reproduction and threshold discussion
- CSF and K562 benchmark construction overview
- Sample-swap critique
- Protein-swap alternative
- Mix_of_Proteome comparison
- Cross-folder synthesis and conclusions

The main remaining discrepancies are:

- Incorrect protein-swap counts and pair count
- Heatmap/NA-heatmap subset mismatch (`good_data` vs requested `all_data`)
- Render failures for the reusable parameterized vignettes

## Suggested priority

1. Correct the protein-swap counts and pair count in `review.qmd`.
2. Decide whether heatmap and NA heatmap figures should follow the TODO
   request (`all_data`) or the current review code (`good_data`), then make the
   code and prose consistent.
3. Fix default path handling in `diagnostics.qmd` and
   `swap_visualization.qmd`.
4. Fix `run_all.fish` skip detection to use the canonical model filenames or a
   shared package-to-filename mapping.
5. Decide whether `run_nonmsstats_block.R` should exit non-zero after any
   package failure.
