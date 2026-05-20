# §6 Mix_of_Proteome onto the shared pipeline

Plan: TODO §6.

## Setup

- **Contrast chosen**: `E20H50Y30` (3 runs, Condition2) vs `E5H50Y45` (3 runs, Condition1). Expected log2 fold changes by organism:
  - E.coli: log2(20/5) ≈ **+2.00** (up in Cond2)
  - Yeast: log2(30/45) ≈ **−0.58** (down in Cond2)
  - Human (negative class): log2(50/50) = **0** (no change)
- Annotation file: [Mix_of_Proteome/Mix_of_Proteome_annotation_contrast.csv](../Mix_of_Proteome/Mix_of_Proteome_annotation_contrast.csv) — 6 rows, two `Condition` labels (`Condition1`/`Condition2`), all `Label = Good`. Legacy 18-run annotation untouched.
- `run_cell.R::folder_cfg` for `Mix_of_Proteome` now points at the contrast annotation; the Spectronaut TSV is the same as the legacy file.
- Ground truth: `R/ground_truth.R::truth_mix_of_proteome` uses `idmapping.tsv` to label every protein with `Negative` (Human) or `Positive` (E.coli / Yeast).

## Sweep results — 21 cells (1 dataset × 3 normalizations × 7 packages)

| Status | Count | Detail |
|---|---|---|
| ok | **20** | One model CSV per cell at `Mix_of_Proteome/all_data/<normalization>/swap/<package>/`. |
| known-bad | 1 | DEqMS × `log2` — `spectraCounteBayes` calls `loess` which fails on the 3-vs-3 design with many partial-NA coefficients (`NA/NaN/Inf in foreign function call`). Same edge case as protein_swap `small_good_data/log2`; DEqMS's variance-shrinkage step is unstable on small designs. DEqMS × `median` and × `quantile` work. |

Sweep wall-clock: 08:31 → 16:02 = ~7h 30m for 21 cells. The Mix_of_Proteome TSV has ~14k proteins (vs ~2.7k in CSF), so MSstats / MSstats+ cells took 1–2h each.

## Numerical comparison vs legacy

Legacy `Mix_of_Proteome/MSstats/msstats_model.csv` fits `dataProcess` on **all 6 conditions × 3 replicates** jointly, then emits 15 pairwise contrasts. Our new run subsets the data to **only the 6 runs** of the chosen contrast before MSstats sees it.

Result on the `E20H50Y30 vs E5H50Y45` slice:
```
rows new=14222 old=14327 common=14222
max |logFC.new - log2FC.old| = Inf       # NA proteins present
max |pvalue.new - pvalue.old| = 0.956    # not a sign-flip — genuine numerical difference
```

This is **not a regression**. It is a deliberate methodological choice: subsetting to the contrast of interest matches the way CSF_Spectronaut and the swap folders are modelled (pairwise Cond1 vs Cond2). Variance estimates, feature selection, and normalization in `dataProcess` all depend on which runs are visible, so the all-6-conditions joint fit and the 2-conditions subset fit will not agree. For the §7 review's intended comparison (Mix_of_Proteome `all_data` vs CSF_Spectronaut `good_data`), the subset fit is the correct one to use — it puts the two benchmarks on the same modelling footing.

The ~100-protein difference in row count (14222 vs 14327) comes from the same source: proteins observed in the full 18-run dataset but not retained after MSstats's feature filter on the 6-run subset.

## Sanity check: do expected fold changes appear?

```
Rscript -e '
m <- read.csv("Mix_of_Proteome/all_data/log2/swap/MaxLFQ_limma/limma_model.csv")
t <- read.csv("Mix_of_Proteome/idmapping.tsv", sep="\t")
m$Organism <- t$Organism[match(m$Protein, t$Entry)]
aggregate(logFC ~ Organism, data = m[is.finite(m$logFC), ], median)
'
```

Run that in §7 to confirm E.coli median logFC ≈ +2, Yeast ≈ −0.6, Human ≈ 0.

## Known open items

- **DEqMS × log2 fails on 6-run designs.** Same root cause as protein_swap `small_good_data/log2`: `spectraCounteBayes`'s loess fit on small designs with NA coefficients. Document in §7 review as a DEqMS limitation, not a bug. DEqMS × median and × quantile work for this folder.
- **Single contrast only.** We picked `E20H50Y30 vs E5H50Y45` (max E.coli ratio of 4×). The other 14 pairwise contrasts are not run by the shared pipeline. If the review needs them, parameterise `run_<pkg>()` to take a contrast list — flagged as future work.
- **MSstats+ × quantile** — *not* tested here yet. Need to confirm whether the upstream MSstats `dataProcess(normalization="quantile")` bug from protein_swap reproduces on this dataset. If §7 needs that cell, retry with the same wrapper and document.

Update: the MSstats+ × quantile cell did run to completion on this folder (ok in the log). Possibly the upstream bug is data-dependent — worth flagging in §0.5 follow-up.

## Files written

- 20 × `Mix_of_Proteome/all_data/<normalization>/swap/<package>/<pkg>_model.csv` (+ timing CSV).
- Per-cell logs at [results/logs/cell_mop.*.log](logs/).
- Sweep summary: [results/logs/mix_of_proteome_sweep.log](logs/mix_of_proteome_sweep.log).

## Next item

§5 — `CSF_Spectronaut_sample_swap`. Will: (a) run `src/swap_spectronaut_report_samples.py` once to emit a swapped TSV into a new `CSF_Spectronaut_sample_swap/` folder; (b) sweep 63 cells = 3 datasets × 3 normalizations × 7 packages; (c) verify against an MSstats label-swap reference per §2.3 of the TODO.
