- We do not generalize the modelling code in `CSF_Spectronaut` — it stays close in structure to the MSstats authors' code. This folder is mainly to replicate the authors' results.

- However, we structure its outputs so they are compatible with the new protein_swap and sample_swap folders.

Shared code that runs on more than one folder lives in `quant/vignettes/` (`.qmd`/`.Rmd`), `quant/R/` (R), and `quant/src/` (Python helpers, e.g. the swap scripts).

The folders `CSF_Spectronaut_protein_swap`, `CSF_Spectronaut_sample_swap` and `Mix_of_Proteome` are analysed by the shared code in `quant/R/` and `quant/vignettes/`.

For `CSF_Spectronaut_protein_swap`: do not compute the noswap — disable that branch for now. Re-enable at the end.

Folder structure (uniform across all four folders):
```
<csffolder>/<dataset>/<normalization>/swap/<modellingPackage>
```
`swap` is literal everywhere; `noswap` is disabled. Example:
```
CSF_Spectronaut_protein_swap/all_data/log2/swap/DEqMS
```

Datasets `<all_data>`, `<good_data>`, `<small>`, `<small_good_data>`:
- `<good_data>` is folder-dependent:
  - `CSF_Spectronaut_protein_swap`: **neat samples only** (TP/TN come from protein swap, not dilution).
  - `CSF_Spectronaut`, `CSF_Spectronaut_sample_swap`: neat + 1/2, balanced replicate counts.
- `<small_good_data>` = 4 vs 4 sampled from `<good_data>`.

Normalization methods: `log2`, `median`, `quantile`. Audit per-package whether "median" is centering-only or centering + z-scaling, and document.

Modelling packages (same set everywhere): MSstats, MSstats+, MaxLFQ+limma, DEqMS, msqrob2 (faster code), prolfqua, limpa.



Regarding what analysis:

# CSF_Spectronaut (`RMSV000000701.3-rerun/quant/CSF_Spectronaut`)

**Normalization scope:** CSF_Spectronaut is used **only to replicate the manuscript's Table 1**, so we use **only log2** (no median, no quantile). The 3-normalization scope (log2, median, quantile) applies only to the shared-pipeline folders: `CSF_Spectronaut_protein_swap`, `CSF_Spectronaut_sample_swap`, `Mix_of_Proteome`.

- `<all_data>` — replicate the authors' results across MSstats, MaxLFQ+limma, and the other packages. Exceptions are msqrob2 (faster code) and prolfqua (new).
- `<good_data>` — show that MSstats+ does not outperform non-MSstats packages because it judges sample quality, but because it does *not* moderate variance.
  - The authors introduce `Mix_of_Proteome` to show that MSstats+ works well on normal datasets. `<good_data>` should be a normal dataset and the modelling tools should behave similarly to how they behave on `Mix_of_Proteome`.

# CSF_Spectronaut_protein_swap (rename of `CSF_Spectronaut_swap`)

We rename the folder to make clear that we swap proteins in group 2, not samples.

Current protein-swap universe: the Spectronaut report has 3,041 non-blank
proteins with at least one precursor, but the protein-swap generator filters to
`n_precursors >= 2` before building the benchmark truth list. This leaves 2,244
proteins: 224 Positives in 112 matched pairs and 2,020 Negatives. The 797
excluded proteins are single-precursor proteins. We keep this conservative
definition because single-precursor proteins cannot participate in the paired
protein-swap construction and are less reliable benchmark Negatives.

- `<all_data>` (as before)
- `<small>` (as before)
- Trash `no_high_dilution`.
- Add `<good_data>`: neat samples only — show that for this dataset the tools behave similarly to `Mix_of_Proteome`.
- Add `<small_good_data>`: 4 vs 4 samples — show that p-value moderation matters.

# CSF_Spectronaut_sample_swap (new)

- In `CSF_Spectronaut`, the sample swap happens in the scripts, not in the data file — which makes it hard to run the same scripts on protein_swap and the original folder. If we start from a Spectronaut file with samples swapped for 90% of proteins (as in `CSF_Spectronaut`), we can run the same scripts as for `CSF_Spectronaut_protein_swap`.
- Use the same conservative `n_precursors >= 2` universe as protein-swap:
  2,244 proteins total, 224 Positives, and 2,020 Negatives. Do not use the full
  3,041 one-or-more-precursor universe for the rebuilt sample-swap benchmark.
- `<all_data>`
- `<good_data>` — all samples created by swapping between neat and 1/2.
- `<small_good_data>` — 4 vs 4 samples; here MSstats+ and MSstats perform even worse: small prior N makes variance moderation kick in harder, hurting the tools.

# Mix_of_Proteome

- Create `<all_data>` subfolder and analyse with the same scripts used for `CSF_Spectronaut_sample_swap` and `CSF_Spectronaut_protein_swap`.

# Documenting what Msstats did and what we did.

- Zenodo - create a zenodo repository for the download RMSV00000701.3 most importantly the quant folder. Point to the msstats+ bioarchive publication and it make it clear that it is to document that state of the scripts as submitted by the authors.

- My code, also Zenodo. However since we work only with the spectronaut data, this part of the quant folder is sufficient.  The easiest way to create a zenodo upload is by creating a github repo. We have git already in the quant. 


--- Review ---

The review should be generated from a qmd file. The qmd file is part of the newly created repository.

Structure:

## Introduction
- introduction about the authors of the article, foremost about olga vitek is and her contribution to the community, about msstats
- Msstats - short overview of msstats (when first puslished) and that it was repeatedly used as a reference in many software publications, e.g. msqrob2 and limpa or prolfqua and it was shown that the other software was performing better. Provide references to these publications.
- what the MSstatsPlus publication is about and what it claims


### Benchmark result Table 1 biorxv and resubmission 2.

- Focus on table 1 in the first draft of the manuscript and then on a quite different table 1 in the second draft of the manuscript. Discuss briefly how the tables were obtained. (p-value thresholds). Reference here the biorxv publication.

UPDATE: 
- Introduce the benchmark dataset, and how they were created both CSF and K562.
 - create a table summarizing the number of samples for each dilution. 
 - Create visualizations explaining the swap schema use figure 1([text](../CSF_Spectronaut/CSF_swap_visualization.qmd))
 - Summarize how many samples per group, and how many <good> and <bad> samples per group.



- Based on the folder CSF_Spectronaut: Present again Table 1 but from our computation. Also mention that numerical differences are because we wanted to speed up the analysis for msqrob2. Show Table 1 for both p-value and FDR filtering, and both for the entire dataset and for the <good_data>  (neat+1/2).
  - Contrast it with the mix_of_proteomes results specifically the part of he <good_data>. While for mix_of_proeomes all the tools perform silimar I expect to see a different picture of the <good_data> from CSF_Spectronaut.
  - Then show figues for the p-value distributions and FDR distributions for <good_data>
  (./quant/CSF_Spectronaut/all_dilutions/V1_log_diagnostics.qmd)
  - we discuss the the assumptions underlying DEA, that for H0 we expect an uniform distribution. 
  - That the p-value histograms obtained from the benchmark dataset is uniptical for <all_data>, but untipical p-value distributions are frequently idicative of outliers
  - Then show the p-value histograms for <good_data> where we expect H0 to be uniformly distributed, point out that it is not, which motivates inspecting the properties of the the swap_samples benchmark dataset.

### The Sample Swap Benchmark dataset

We show figure 1 (msstatsplus/RMSV000000701.3-rerun/quant/CSF_Spectronaut/CSF_swap_visualization.qmd) Swapping schema.

- Next based on the the CSF_spectronaut_sample_swap we visualize the dataset: 
  - we show (using prolfqua):
  - heatmap <all_data>
  - NA_heatmap <all_data> 
  - density plots of the <good_data> 
- Then all the figures we have in "Per-protein effect size and within-condition variance" see file (msstatsplus/RMSV000000701.3-rerun/quant/CSF_Spectronaut/CSF_swap_visualization.qmd)
  - We specifically discuss the diffference in variance for H0 and H1, and how this affects methods which perform variance shrinkage. 
  - Give some estimates how much the sample size impacts the effect of the shrinkage. few samples - then prior is strong, more samples prior is week. 
  - Show table 1 for the CSF_spectronaut_sample_swap <good_data> and <small_good_data> to illustrate the effect of the sample size.

### The Protein Swap Benchmark dataset


- Describe how we created CSF_Spectronaut_protein_swap, an alternative benchmark dataset where we swap the proteins instead of samples.

- State explicitly that the protein-swap truth list uses the conservative
  `n_precursors >= 2` universe: 2,244 proteins total, 224 Positives in 112
  matched pairs, and 2,020 Negatives. The rebuilt sample-swap truth list uses
  the same 2,244-protein universe, but without pairing Positives.

Crete a figure illustrating the protein swapping.

- Next based on the the CSF_spectronaut_protein_swap we visualize the dataset: 
  - we show (using prolfqua):
  - heatmap <all_data>
  - NA_heatmap <all_data> 
  - density plots of the <good_data> 
- Then all the figures we have in "Per-protein effect size and within-condition variance" see file (msstatsplus/RMSV000000701.3-rerun/quant/CSF_Spectronaut/CSF_swap_visualization.qmd) but for the CSF_spectronaut_proein_swap dataset

# Benchmark results for the Protein Swap dataset

- Compare results of `<good_data>` for protein_swap with `Mix_of_Proteome`. Expect a similar performance profile across models for both (contrast it with the performance profile of `<good_data>` from sample_swap).

- Show performance results for `<all_data>` and `<small>` for the protein_swap dataset and contrast with results from the sample_swap dataset.


# Conclusions

