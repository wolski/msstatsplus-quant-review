# Protein-swap good_data redo (subset definition fix)

## What was wrong

`R/dataset_subsets.R::subset_good_data` filtered by `Label == "Good"` for
every folder. For `CSF_Spectronaut_protein_swap` that yields 29 runs (neat +
1:2). Per WhatWeStillNeed.md §0.2 and the design rationale, the correct
subset for protein_swap is **neat-only (15 runs: 8 Cond1 + 7 Cond2)** —
because Cond1/Cond2 in this folder is the G1/G2 random split *within each
dilution*, and including the 1:2 runs adds a 2× dilution signal unrelated
to the protein-swap that contaminates within-condition variance.

The 42 cell outputs under `CSF_Spectronaut_protein_swap/good_data/` and
`small_good_data/` were therefore all computed against the wrong subset.

## Fix

`R/dataset_subsets.R::subset_good_data` is now folder-aware:

```r
if (identical(csffolder, "CSF_Spectronaut_protein_swap")) {
  keep <- grepl("NeatCSF", annotation$R.FileName)
  annotation[keep, , drop = FALSE]
} else {
  annotation[annotation$Label == "Good", , drop = FALSE]
}
```

`subset_small_good_data` already wraps `subset_good_data`, so it picks up
the fix automatically.

`vignettes/review.qmd` §7.6 effect-size/SD/mean-variance chunk switched from
`all_data/log2` rda → `good_data/log2` rda. The heatmap and NA-heatmap
chunks intentionally remain on `all_data` (they illustrate the dataset-level
swap pattern, which is more visible across the full 40-run set).

## Sweep results

42 cells re-run via [run_folder_sweep.fish](../run_folder_sweep.fish):

| Status | Count | Detail |
|---|---|---|
| ok   | **40** | All 5 non-MSstats packages × 6 (dataset, norm) blocks = 30, plus 10 of 12 MSstats/MSstats+ cells. |
| fail | 2 | **MSstats+ × good_data × quantile** and **MSstats+ × small_good_data × quantile** — the upstream `MSstats::dataProcess(normalization="quantile")` bug, observed previously on `all_data × quantile` for MSstats+. Now confirmed for all three protein_swap subsets. |

Wall-clock: 23:17 → 00:19 = **62 min**. Previous (wrong-subset) sweep took
~3h41m without `n_cores=8` + bundle; the new infrastructure paid off.

Notable: **DEqMS × small_good_data × log2 worked this time** — the previous
sweep failed it with the `spectraCounteBayes` loess error. The neat-only
4v4 design has a different missingness pattern than the previous wrong
neat+1/2 4v4 (probably fewer all-NA rows), enough for loess to converge.

## Table 1 — protein_swap good_data × log2 (neat-only, 15 runs)

| Package | N | TPR (p<0.05) | PPV (p<0.05) |
|---|---:|---:|---:|
| MSstats | 2673 | 0.873 | 0.576 |
| MSstats+ | 2671 | 0.863 | 0.570 |
| MaxLFQ_limma | 2673 | 0.896 | 0.667 |
| DEqMS | 2673 | 0.868 | 0.637 |
| msqrob2 | 2673 | 0.906 | 0.561 |
| prolfqua | 2673 | 0.877 | 0.592 |
| limpa | 2673 | 0.882 | 0.512 |

**Methods cluster tightly** (TPR 0.86–0.91; PPV 0.51–0.67). MSstats and
MSstats+ no longer outperform the moderated methods — they're in the middle
of the pack. This matches WhatWeStillNeed.md's stated expectation:

> I expect to see a similar performance profile for the models for both
> datasets [protein_swap good_data and Mix_of_Proteome], contrast it with
> the performance profile of the swap of sample dataset good_data.

The story now visible across the three CSF subsets:

- **sample_swap good_data**: MSstats+ TPR=0.983 / PPV=1.000 vs limpa
  TPR=0.703 / PPV=1.000 — large spread, MSstats+ wins by a wide margin.
  The variance-asymmetry biases moderated methods.
- **protein_swap good_data**: MSstats+ TPR=0.863 / PPV=0.570 vs limpa
  TPR=0.882 / PPV=0.512 — tight cluster, MSstats+ no longer leads. The
  variance-asymmetry runs the other direction (Positives inflated, not
  Negatives), and the small N (15 runs) reduces all PPVs.
- **Mix_of_Proteome all_data**: also tight cluster (similar to
  protein_swap good_data).

The variance-asymmetry story is what creates MSstats+'s win in the
sample_swap benchmark, not method quality.

## Table 1 — protein_swap small_good_data × log2 (4v4 from neat-only)

| Package | N | TPR (p<0.05) | PPV (p<0.05) |
|---|---:|---:|---:|
| MSstats | 2664 | 0.755 | 0.586 |
| MSstats+ | 2664 | 0.759 | 0.528 |
| MaxLFQ_limma | 2673 | 0.802 | 0.515 |
| DEqMS | 2673 | 0.731 | 0.517 |
| msqrob2 | 2673 | 0.783 | 0.462 |
| prolfqua | 2673 | 0.755 | 0.516 |
| limpa | 2673 | 0.802 | 0.463 |

Still tight clustering at 4v4. TPR drops uniformly by ~0.10 vs the full
good_data; PPV drops less. Moderated methods do not pull ahead at small N
for protein_swap (unlike sample_swap, where they get penalised by the
opposite variance asymmetry).

## What needs the user's attention next

- `review.qmd` §7.6 / §7.7 narrative around protein_swap currently still
  reads as if the comparison was on neat+1/2 (the old wrong subset). The
  underlying figures and tables are correct now, but a paragraph or two of
  narrative explaining the **neat-only** good_data choice and the
  expected/observed tight clustering would tighten the story.
- MSstats+ × quantile fails on every protein_swap subset (`all_data`,
  `good_data`, `small_good_data`). Other folders have the same bug only on
  some subsets. Worth a stand-alone note in the methods section: an
  upstream MSstats limitation, not a benchmark artefact.

## Next item

§8 — Zenodo + GitHub publication.
