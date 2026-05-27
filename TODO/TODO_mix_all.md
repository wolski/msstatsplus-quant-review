# Mix_of_Proteome — multi-condition, all-pairwise contrasts

Plan for re-enabling `Mix_of_Proteome/` as a standalone pipeline that fits
**one model per package over all 6 conditions** and computes **all 15
pairwise contrasts**. Authoritative reference for the canonical analysis is
[`../RMSV000000701.3/quant/Mix_of_Proteome/Mixture_of_proteomes_processing.R`](../../RMSV000000701.3/quant/Mix_of_Proteome/Mixture_of_proteomes_processing.R)
in the read-only reference deposit.

## Design principle

Mix gets its own runner script and stays fully isolated from the swap-pipeline
adapters. **Nothing in `R/models_*.R` or `R/run_*_block.R` changes.** The
existing `mix` alias at the top-level `Makefile` continues to recurse into
`Mix_of_Proteome/Makefile` (committed in `7c76cc0`); only the body of that
folder Makefile and the new R sources change.

## Discovery (already done, recorded here so the next session does not
have to re-investigate)

- `Mix_of_Proteome/Mix_of_Proteome_annotation.csv` — full 6-condition design
  (E5H50Y45 … E45H50Y5), 18 runs, 3 BioReplicates per condition.
- `Mix_of_Proteome/Mix_of_Proteome_annotation_contrast.csv` — collapsed to
  2 conditions (`Condition1`, `Condition2`), 6 runs labelled `Good`. This
  is what `annotation.csv` is currently symlinked to and is what blocks the
  multi-condition design.
- The current rerun adapters at `R/models_maxlfq_limma.R:39`,
  `R/models_deqms.R:66`, `R/models_limpa.R:57` all hard-code
  `limma::makeContrasts(classCondition2 - classCondition1, …)`. They cannot
  do multi-condition without modification.
- The reference processing script fits MSstats+, MSstats, msqrob2,
  limma + MaxLFQ, limpa, DEqMS, mapDIA (prep only). prolfqua is **not** in
  the original — it would be an addition specific to the rerun.
- Ground-truth comes from organism membership (Human = null; E. coli, Yeast
  = differentially abundant), and the **expected log2FC per contrast** is
  computed from the percentage mix design in
  `Mix_of_Proteome/Mix_of_Proteome_design.md` and in the reference script at
  [`Mixture_of_proteomes_processing.R:946-951`](../../RMSV000000701.3/quant/Mix_of_Proteome/Mixture_of_proteomes_processing.R#L946-L951)
  (15-element `human_ratios`, `ecoli_ratios`, `yeast_ratios` vectors).
- One previously-cited blocker (Makefile header: "DEqMS spectraCounteBayes
  too small for 3 vs 3") **does not apply** in the all-pairwise design:
  6 conditions × 3 reps = 18 runs gives DEqMS enough degrees of freedom.

## Files to add / change (all inside `quant/`)

### New files

1. **`R/mix_contrasts.R`** — single source of truth for the 15-contrast
   spec.
   - `comparisons` (15 human-readable labels, e.g. `"E10H50Y40 vs E20H50Y30"`),
   - `limma_comparisons` (15 strings in `classXXX - classYYY` form for limma /
     limpa / DEqMS / MaxLFQ),
   - `msqrob_contrasts` (15 strings in the asymmetric `condition…` /
     `condition… - condition…` form msqrob2 needs — first 5 are direct level
     coefficients vs the implicit reference, last 10 are differences),
   - `human_ratios`, `ecoli_ratios`, `yeast_ratios` (length-15 numeric vectors,
     expected per-contrast log2FC by organism),
   - `comparison_map` data.frame joining all the above on `Label`.
   - The msqrob2 sign-flip for the first 5 contrasts (reference deposit
     [`L1041-L1047`](../../RMSV000000701.3/quant/Mix_of_Proteome/Mixture_of_proteomes_processing.R#L1041-L1047))
     belongs here too, exposed as a helper function.

2. **`R/run_mix_processing.R`** — ported from the reference deposit's
   `Mixture_of_proteomes_processing.R`, ≈ 600 lines, with these adjustments:
   - CLI signature: `Rscript R/run_mix_processing.R <subset_dir> <normalization>`
     (mirrors `R/run_cell.R` and `R/run_nonmsstats_block.R`).
   - Replace the `data_folder = ""` constant with an argv-driven `subset_dir`.
   - Each package writes to `<subset_dir>/<norm>/swap/<pkg>/<pkg>_model.csv`,
     same path convention as the swap pipeline (so downstream comparison
     helpers can reuse `R/comparison_table.R` patterns).
   - Use `R/preprocess.R` normalization helpers so the `<normalization>` arg
     (log2 / median / quantile) actually changes the inputs. The reference
     script is log2-only.
   - Wrap each package in `tryCatch`; treat per-package failures as
     non-fatal *for the standalone Mix pipeline* (Mix is a research target,
     not a release-gate). Print a summary at the end.
   - Add prolfqua as a 7th adapter mirroring its package-native multi-contrast
     API. **If it doesn't fit cleanly in one session, commit without it and
     follow up.**
   - Drop mapDIA (the reference script only prepped the input but never
     invoked the binary).
   - Parameterize `numberOfCores` via env var
     (`MIX_CORES`, default `parallel::detectCores() / 2`).

### Changed files

3. **`Mix_of_Proteome/Makefile`** — rewrite. Today it has 3 cell rules
   driving `R/run_cell.R` + `R/run_nonmsstats_block.R`. After:
   ```make
   include mk/common.mk

   MIX_DEPS = R/run_mix_processing.R R/mix_contrasts.R \
              R/preprocess.R R/paths.R

   # Symlinks: switch annotation.csv to the 6-condition version.
   symlinks-mix: Mix_of_Proteome/Report.tsv Mix_of_Proteome/annotation.csv
   Mix_of_Proteome/annotation.csv:
   	ln -s Mix_of_Proteome_annotation.csv $@

   # Prep: same as today but driven off the 6-condition annotation.
   prep-mix: Mix_of_Proteome/.prep.stamp
   Mix_of_Proteome/.prep.stamp: …  # build_subsets.py with --subsets all_data

   # One cell per normalization. Single processing script fits all packages.
   cells-mix-all_data-<norm>: Mix_of_Proteome/all_data/<norm>/swap/.stamp
   Mix_of_Proteome/all_data/<norm>/swap/.stamp: $(MIX_DEPS) Mix_of_Proteome/.prep.stamp
   	Rscript R/run_mix_processing.R Mix_of_Proteome/all_data <norm>
   	touch $@
   ```
   - Bundle / msstats stamp split is dropped (one script, one stamp per cell).
   - Diagnostics rules are removed (see note below).

4. **`Makefile` (top-level)** — drop the `mix-diagnostics` alias and its
   `make help` line. The swap-style `diagnostics.qmd` is not Mix-aware
   (organism-based truth and per-contrast errors are a different shape),
   so leaving the alias around would be misleading. Reintroduce once
   `vignettes/diagnostics_mix.qmd` exists. The `mix` alias stays.

5. **`Mix_of_Proteome/Mix_of_Proteome_annotation.csv`** — unchanged on disk;
   `annotation.csv` symlink is repointed by the new symlinks rule.

### Explicitly NOT in scope (defer to a follow-up)

- **`vignettes/diagnostics_mix.qmd`** — new Mix-aware diagnostics rendering
  the per-contrast / per-organism error and AUROC plots from the reference
  analysis script (`Mixture_of_proteomes_analysis.R`).
- **Generalizing the swap adapters** to accept a `contrasts` parameter.
  Possible cleanup later; out of scope here because of the risk to the
  swap pipeline.
- **Adding mapDIA** back. The original never ran it from R.

## Risks

- **prolfqua port** — package-specific work to emit per-contrast long rows.
  Fallback: commit without prolfqua.
- **DEqMS spectraCounteBayes** — the prior comment about insufficient
  degrees of freedom was for a 3 vs 3 contrast; should be fine in 6 × 3 =
  18 runs. Verify on first run.
- **Runtime** — the reference script asks for `numberOfCores = 12` in
  MSstats. Parameterize via env var.
- **One-shot stamp** — a per-package failure aborts only its own try-block,
  not the whole stamp. The stamp is touched at the end regardless. That
  matches Mix's research-target status (versus the swap pipeline's
  release-gate posture).

## Order of operations once approved

1. `R/mix_contrasts.R` (small, pure data).
2. `R/run_mix_processing.R` (the heavy port + parameterize).
3. Rewrite `Mix_of_Proteome/Makefile`.
4. Repoint `annotation.csv` symlink.
5. Drop the `mix-diagnostics` alias from the top-level Makefile + `make help`.
6. Commit as one logical change.
7. User runs `make mix` (likely after the live `make all` build frees up
   cores).
8. Iterate on prolfqua / diagnostics in a follow-up.

## State at the time this plan was written (2026-05-21)

- The current `make -j5 -k all` rebuild is still in flight on the user's box
  (initial sample-swap cells + downstream).
- Pre-requisites already in place: `mix` and `mix-diagnostics` aliases in
  the top-level Makefile (`7c76cc0`), `msqrob2` enabled in `$(SCRIPTS)` and
  in the swap-folder bundles (`ceb6ae7`).
- The plan is **not implemented yet**. User wants to first run what is
  currently committed before deciding whether and how to start on this.
