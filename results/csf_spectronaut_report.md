# §3 CSF_Spectronaut — Table 1 reproduction (all_data × log2 × 7 packages)

Plan: TODO §3, narrowed to *only* `all_data × log2 × 7 packages` per user direction. The full §3.3–§3.5 grid (2 datasets × 3 normalizations) is not needed — the sample_swap folder (§5) already covers the variance-shrinkage story for `good_data` and the alternative normalizations.

## What was done

- **No rerun.** The legacy CSF_Spectronaut author scripts had already written `MSstats/`, `MSstats+/`, `limma/`, `DEqMS/`, `msqrob2/`, `prolfqua/`, `limpa/` under `CSF_Spectronaut/all_dilutions/V1_log2/`. These are the author-faithful outputs of the original benchmark.
- **Path reshape only.** Copied each legacy package directory into the canonical §0.1 layout, with the one rename (`limma → MaxLFQ_limma`):

```
CSF_Spectronaut/all_dilutions/V1_log2/<legacy>/   →   CSF_Spectronaut/all_data/log2/swap/<canonical>/
  MSstats        MSstats
  MSstats+       MSstats+
  limma          MaxLFQ_limma
  DEqMS          DEqMS
  msqrob2        msqrob2
  prolfqua       prolfqua
  limpa          limpa
```

The two preswap variants (`MSstats_preswap/`, `MSstats+_preswap/`) are not part of the canonical layout — they live alongside in `all_dilutions/V1_log2/` as alternative methodological checks. Per §1 plan, the noswap branch is disabled in the new layout.

## Table 1 — Manuscript reproduction

P-value threshold = 0.05; ground truth from [CSF_protein_swap_list.csv](../CSF_Spectronaut/CSF_protein_swap_list.csv) (182 Positives / 1638 Negatives).

| Package | N matched | TPR | PPV | Design doc target (§7) |
|---|---:|---:|---:|---|
| **MSstats** | 1786 | **0.777** | **1.000** | 0.777 / 1.000 ✓ |
| **MSstats+** | 1786 | **0.983** | **0.994** | 0.983 / 0.994 ✓ |
| MaxLFQ_limma | 1786 | 0.880 | 1.000 | — |
| DEqMS | 1786 | 0.857 | 0.993 | — |
| msqrob2 | 1786 | 0.966 | 1.000 | — |
| prolfqua | 1786 | 0.971 | 0.994 | — |
| limpa | 1786 | 0.714 | 1.000 | — |

**MSstats and MSstats+ match the design doc's three-decimal targets exactly.** Other five packages weren't tabulated in the design doc; their numbers here are the reference values for the §7 review.

## Cross-check vs sample_swap (§5)

Both folders run on the same underlying CSF data; the swap stage differs:

- **CSF_Spectronaut** (this report): in-script post-summarization label swap (`run_msstats.R` → `swap_condition_labels` on `ProteinLevelData`).
- **CSF_Spectronaut_sample_swap** (§5): pre-`SpectronauttoMSstatsFormat` TSV-level rewrite.

| Package | CSF (post-summarization) | sample_swap (TSV-level) | Δ |
|---|---|---|---|
| MSstats | 0.777 / 1.000 | 0.777 / 1.000 | 0 |
| MSstats+ | 0.983 / 0.994 | 0.983 / 1.000 | +1 row Negative-near-threshold |
| MaxLFQ_limma | 0.880 / 1.000 | 0.903 / 0.994 | +0.023 TPR |
| DEqMS | 0.857 / 0.993 | 0.863 / 0.993 | +0.006 TPR |
| msqrob2 | 0.966 / 1.000 | 0.977 / 0.872 | **−0.128 PPV** |
| prolfqua | 0.971 / 0.994 | 0.971 / 0.988 | −0.006 PPV |
| limpa | 0.714 / 1.000 | 0.703 / 1.000 | −0.011 TPR |

The two swap stages produce essentially equivalent benchmarks except for **msqrob2's PPV** (1.000 vs 0.872). Hypothesis: msqrob2's mandatory peptide-level `QFeatures::normalize(method = "center.median")` (see [R/README_normalization.md](../R/README_normalization.md)) interacts with the swap stage — when intensities are pre-swapped at TSV level, the peptide-level centering averages over a heterogeneous (mixed-condition) pool for Negative proteins, biasing their estimated effect away from zero. Worth flagging in the §7 review as a real methodological note.

## Files (re-)written

- `CSF_Spectronaut/all_data/log2/swap/<package>/*.csv` for all 7 packages (copies of legacy outputs, no recompute).

## What's not in §3

- `good_data`, `median`, `quantile` cells: not run for CSF_Spectronaut. Per user direction the swap_swap folder already exercises those.
- `_preswap` outputs left at their legacy paths (alternative-stage comparator; not in the canonical layout).

## Next item

§7 — populate `vignettes/review.qmd` against the results now sitting in `results/` and the cell outputs across the four folders. This is the last analytical step before §8 (Zenodo + GitHub publication).
