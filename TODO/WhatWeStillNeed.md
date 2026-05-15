- We do not generalize the modelling code in CSF_spectronaut, it should stay close in structure to the code of the Msstats authors. This folder is mainly to replicate the results of the authors.


- However, we want to structure the outputs in such a way that they are compatible with the new protein_swap and sample_swap folders.


All code that can run on more then one folder, e.g. should go into subfolder quant vignettes (qmd and rmd files) and and R folder for the R files.

The folders CSF_spectronaut_protein_swap and CSF_spectronaut_protein_swap and Mixture_of_Protoeomes, should be analysed by the code from the R folder


Regarding the CSF_Specnaut_protein_swap folder -  do not compute the noswap please, lets disable that branch, for the moment. We will enable it in the end.

Folder structure:
The folder strucure should be 
small_set/V1_.../...

that is: 
<csffolder>/<dataset>/<normalization>/<swap/noswap>/<modellingPacakege>

for instance:
CSF_Spectronaut_protein_swap/all_dilutions/log2/noswap/DEqMS



Regarding what analysis:

# CSF_spectronaut (/Users/wolski/projects/reviews/msstatsplus/RMSV000000701.3-rerun/quant/CSF_Spectronaut)

- one analysis with the entire dataset <all_data> - to show that we replicate the results for MSstats, MaxLFQLimma, and the other packages. Exception is msqrob2 (faster code) and proflqua (new)

- one analysis with only the good samples <good_data> -> to show that MSstats+ does not perform better then other non msstats package because it analyses sample quality but because, it does not moderate variance.

  -  the reason is that they introduce the Mix_of_proteomes to show that msstats+ works good on normal datasets. the good_data should be a normal dataset and the behaviour of the modelling tools should be similar then on the mix of proetomes dataset.

# CSF_spectronaut_protein_swap (/Users/wolski/projects/reviews/msstatsplus/RMSV000000701.3-rerun/quant/CSF_Spectronaut_swap)

we rename this folder to make it clear that we swap the proteins in group 2 not the samples.

- We run <all_data> (as before)
- We run <small> (as before)
- We trash the no_high_dilution
- We add a dataset <good_data> with only the neat_samples again, to show that for this dataset the tools behave similarily as for the mix_of_proteomes data.
- And we also create <small_good_data> with 3 vs 3 samples to show, that p-value moderations matters


The normalization methods we try is log2, median_scaling and quantile

# CSF_spectronaut_sample_swap (new)

 - the problem with CSF_spectronaut is that the sample swap is happening in the scripts not in the file, this makes it difficult to run the same scripts on the protein_swap and the original folder. However if we start from a spectronaut file with swapped samples for 90% of the proteins as in the CSF_spectronaut, we can run the same scripts which we run CSF_spectronaut_protein_swap

- We run <all_data>
- We run <good_data> - all the samples created by swapping between neat and 1/2
- We run <small_good_data> - 3 vs 3 samples, and here we show, that MSstats+ and MSstats perform even better - because here because of small prior N the variance moderation kicks in much harder. Making the tools behave worse.


# Mix_of_proteomes

- We create subfolder <all_data> and analyse the data using the scripts which we used to model CSF_spectronaut_sample_swap, CSF_spectronaut_protein_swap

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


### Benchmark result Table 1 biorxv and resubmission2.

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

Crete a figure illustrating the protein swapping.

- Next based on the the CSF_spectronaut_protein_swap we visualize the dataset: 
  - we show (using prolfqua):
  - heatmap <all_data>
  - NA_heatmap <all_data> 
  - density plots of the <good_data> 
- Then all the figures we have in "Per-protein effect size and within-condition variance" see file (msstatsplus/RMSV000000701.3-rerun/quant/CSF_Spectronaut/CSF_swap_visualization.qmd) but for the CSF_spectronaut_proein_swap dataset

# Benchmark results for the Protein Swap dataset

- compare results of <all_good> for protein_swap with those of the mixture_of_proteomes. I expect to see a similar performance profile for the models for both datasets (contrast it with the preformance profile of the swap of sample dataset <all_good>). Include it here.

- Show the performance result for <all_data> and <small> for the protein_swap dataset and contrast them with the results obtained witht the sample_swap dataset.


# Conclusions



