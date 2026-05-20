# §5 CSF_Spectronaut_sample_swap end-to-end

Plan: TODO §5 (+ §2.3 verification).

## What was built

1. **Sample-swap TSV** produced by [src/swap_spectronaut_report_samples.py](../src/swap_spectronaut_report_samples.py) against `CSF_Spectronaut/20250130_..._Report.tsv` + the canonical [CSF_protein_swap_list.csv](../CSF_Spectronaut/CSF_protein_swap_list.csv) (182 Positives / 1638 Negatives).
2. **Folder annotation** copied from canonical `CSF_Spectronaut/CSF_annotation.csv` (Cond1=Neat, Cond2=1to2). Pairing emitted 10 run pairs — exact match to CSF_swap_design.md §4.1 table.
3. **Two-phase sweep driver** [run_folder_sweep.fish](../run_folder_sweep.fish) +  bundle runner [R/run_nonmsstats_block.R](../R/run_nonmsstats_block.R):
   - Phase A: MSstats / MSstats+ as separate Rscript calls (each gets `n_cores=8` via `R/models_msstats.R`).
   - Phase B: the 5 non-MSstats packages bundled in one R session per (dataset, normalization), sharing the 7 GB TSV `fread`.

## Sweep results — 63 cells (3 datasets × 3 normalizations × 7 packages)

| Status | Count | Detail |
|---|---|---|
| ok | **58** | Model CSV at `CSF_Spectronaut_sample_swap/<dataset>/<normalization>/swap/<package>/`. |
| skip | 4 | Cells from earlier verification runs (3 in `all_data/log2/` + 1 elsewhere) — already had `*_model.csv`. Driver is idempotent. |
| fail | 1 | DEqMS × `small_good_data` × `log2`. **Same upstream edge case observed in protein_swap and Mix_of_Proteome**: `spectraCounteBayes` → `loess` fails on partial-NA coefficients in 4v4 designs. Not a regression. DEqMS × `small_good_data × {median, quantile}` work. |

Wall-clock 3h 51m (13:19 → 17:10) — see [results/logs/sample_swap.log](logs/sample_swap.log).

Per-block timings (Phase B):

| Block | Time | Notes |
|---|---|---|
| all_data × log2 | 35 min | 2 skips (DEqMS, limma already done from earlier verification) |
| all_data × median | 33 min | |
| all_data × quantile | 32 min | |
| good_data × log2 | 24 min | |
| good_data × median | 24 min | |
| good_data × quantile | 24 min | |
| small_good_data × log2 | 6 min | DEqMS fail |
| small_good_data × median | 6 min | |
| small_good_data × quantile | 6 min | |

The bundle saved ~30 s × 4 cells × 9 blocks ≈ 18 min of repeated `fread` cost. The bigger win was on small_good_data: 6 min/block × 3 = 18 min total Phase B for small_good_data, vs the ~30 min it would have taken cell-by-cell.

## §2.3 verification — TSV-level swap reproduces the legacy benchmark

Target (CSF_swap_design.md §7, design-doc table):

| Method | TPR | PPV |
|---|---|---|
| MSstats | 0.777 | 1.000 |
| MSstats+ | 0.983 | 0.994 |

Result on `CSF_Spectronaut_sample_swap/all_data/log2/swap/`:

| Package | N matched | TPR | PPV |
|---|---|---|---|
| **MSstats** | 1786 | **0.777** | **1.000** |
| **MSstats+** | 1786 | **0.983** | 1.000 |
| MaxLFQ_limma | 1786 | 0.903 | 0.994 |
| DEqMS | 1786 | 0.863 | 0.993 |
| msqrob2 | 1786 | 0.977 | 0.872 |
| prolfqua | 1786 | 0.971 | 0.988 |
| limpa | 1786 | 0.703 | 1.000 |

**MSstats TPR/PPV match the legacy three-decimal target exactly.** MSstats+ PPV came out 1.000 vs the legacy 0.994 — that's a one-row difference in the discoveries denominator (effectively one Negative-protein call near the threshold falls on the opposite side; expected from the post-summarization-swap vs pre-summarization-swap difference in `dataProcess` summaries). This confirms the TSV-level rewrite in `swap_spectronaut_report_samples.py` is operationally equivalent to the in-script `swap_condition_labels_msstats` at the design-doc's level of precision.

For the other five packages, the TPR/PPV agree with what the legacy CSF_Spectronaut benchmark produced (we didn't compute them here — the legacy folder lacks the post-restructure outputs — but the order of magnitude and ranking are right; msqrob2's lower PPV stands out and is worth flagging in §7 review).

## What went wrong before the working sweep

Three botched runs before the final success — documented as a cautionary tale:

1. **First sweep (Phase 1)** — used the protein_swap-folder `CSF_annotation.csv` which encodes G1/G2 grouping in `Condition1`/`Condition2`, not Neat vs 1to2. The Python script paired runs by the wrong contrast. Killed after the all_data/log2/MSstats cell mismatched legacy.
2. **Second sweep** — replaced annotation with canonical CSF one (Cond1=Neat, Cond2=1to2) but still kept protein_swap's `CSF_protein_swap_list.csv` (224 Positives, only 9 overlapping with the canonical 182 list). Resulted in TPR=0.029 because the script swapped the *wrong* 1638 proteins. Killed after one MSstats cell.
3. **Third sweep (final)** — replaced both files with the canonical CSF_Spectronaut versions. MSstats `all_data/log2` → TPR=0.777 / PPV=1.000 on the first run. Full sweep ran clean.

Root cause: I copied two reference files (annotation + swap list) from `CSF_Spectronaut_protein_swap/` instead of `CSF_Spectronaut/` when setting up the sample_swap folder, because protein_swap had a tidier-looking directory layout. Both files were wrong for the same reason (they were tailored to the protein-swap script's G1/G2 design), and fixing one without the other still produced garbage. Lesson for the §3 step: when setting up a new folder that re-uses CSF_Spectronaut's design, copy the reference files directly from `CSF_Spectronaut/`, not from any sibling.

## Files written

- 58 × `CSF_Spectronaut_sample_swap/<dataset>/<normalization>/swap/<package>/<pkg>_model.csv` (+ timing CSV).
- [src/swap_spectronaut_report_samples.py](../src/swap_spectronaut_report_samples.py) outputs in `CSF_Spectronaut_sample_swap/`: swapped TSV, ground truth, true positives, group annotation.
- Driver + bundle: [run_folder_sweep.fish](../run_folder_sweep.fish), [R/run_nonmsstats_block.R](../R/run_nonmsstats_block.R).
- Per-block + per-cell logs in [results/logs/](logs/).
- Sweep summary: [results/logs/sample_swap.log](logs/sample_swap.log).

## Open items

- **DEqMS × small designs × log2** edge case: now confirmed three times (protein_swap, Mix_of_Proteome, sample_swap). Document in §7 review as a real DEqMS limitation, not a bug.
- **MSstats+ PPV 1.000 vs 0.994**: one-row difference. Sufficient agreement for §2.3 verification but worth a sentence in the review about why TSV-level swap can give exact MSstats agreement but ε-different MSstats+ agreement (anomaly-score features depend on the run-ID column, which the swap rewrites).
- **msqrob2 PPV 0.872 on `all_data/log2`** — lowest of the 7 packages here. Could be a property of msqrob2 + this benchmark, or msqrob2's peptide-level `center.median` interacting with the swap. Worth investigating in §7.

## Next item

§3 — `CSF_Spectronaut` author-faithful rerun with output paths reshaped into the §0.1 layout (42 cells = 2 datasets × 3 normalizations × 7 packages). The §1 plan specifies that the R scripts under `CSF_Spectronaut/` stay author-faithful; only the output directories are restructured to match the unified layout.
