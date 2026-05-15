# CSF benchmark: dataset design and the role of the run/label "swap"

This document explains how the CSF Spectronaut (and DIANN) benchmark builds a
ground truth of differentially abundant ("Positive") and non-differentially
abundant ("Negative") proteins, and what side effects the swap step has on the
benchmark.

The swap **removes** the dilution signal for the proteins that are supposed to
serve as the null class — it does not introduce one.

---

## 1. The raw experiment

The CSF dataset is a **dilution series** of human cerebrospinal fluid acquired
on a single LC-MS platform. Each MS run is a different dilution of the same
underlying CSF pool (`CSF_annotation.csv`):

| Sample type | Dilution | # runs | Benchmark assignment |
|---|---|---|---|
| Neat CSF | 1× | 20 | Condition1 (mostly) |
| 1:2 CSF | 0.5× | 16 | Condition2 (mostly) |
| 1:4, 1:8, 1:16, 1:32, 1:64 | various | 2 each | Split between Cond1/Cond2 (`Bad`) |
| Blank | — | 7 | Excluded |

The **canonical comparison** is **Neat (Cond1) vs 1:2 (Cond2)**. Runs are
*interleaved in injection order* (the `Order` column) — Neat and 1:2 alternate
roughly every 1-2 injections — which spreads any temporal drift evenly across
both conditions. The swap exploits this ordering.

Because both conditions are technical dilutions of the **same** biological
material, every detected protein should differ by ~2× (log2 ≈ 1) between
conditions. **In the raw data every protein is DE — there is no intrinsic
null class.**

---

## 2. Why a null class is needed

TPR and PPV require both true positives and true negatives. The raw dilution
series only provides positives (everything is DE). The swap manufactures the
missing null class out of the same real data, deterministically.

---

## 3. The protein swap list

[`CSF_protein_swap_list.csv`](CSF_protein_swap_list.csv) is a fixed list of
1820 proteins:

- **182 Positive** (~10%): designated **true DE** in the post-swap dataset.
- **1638 Negative** (~90%): designated **true non-DE** in the post-swap dataset.

The labels describe the **intended ground-truth status** after the swap, not
the operation. The swap is applied to the 90% Negative subset; the 10%
Positive subset is left alone.

---

## 4. The swap operation

The clearest implementation is `prepare_data_for_msqrob` in
[`benchmark_experiments_functions.R`](../benchmark_experiments_functions.R);
the new `swap_condition_labels_msstats` mirrors it.

### 4.1 Pairing runs

Sort runs by `Order`. Within each condition, pick every second run
(`seq(2, n, by = 2)`); pair them index-for-index across conditions. For this
dataset (annotated in [`CSF_annotation_with_Swap.csv`](CSF_annotation_with_Swap.csv)):

| Pair | Cond1 run | Cond2 run | Type |
|---|---|---|---|
| 1 | Seq2 (Neat) | Seq4 (1:2) | Good |
| 2 | Seq6 (Neat) | Seq9 (1:2) | Good |
| 3 | Seq10 (Neat) | Seq13 (1:2) | Good |
| 4 | Seq15 (Neat) | Seq17 (1:2) | Good |
| 5 | Seq18 (Neat) | Seq21 (1:2) | Good |
| 6 | Seq22 (Neat) | Seq25 (1:2) | Good |
| 7 | Seq26 (Neat) | Seq29 (1:2) | Good |
| 8 | Seq32 (1:2 in Cond1) | Seq31 (1:4) | Bad |
| 9 | Seq36 (1:8 in Cond1) | Seq35 (1:16) | Bad |
| 10 | Seq40 (1:32 in Cond1) | Seq39 (1:64) | Bad |

10 of the 20 Cond1 runs and 10 of the 20 Cond2 runs are paired (50% in each
condition). The other half stays unpaired — and that is essential to the
mechanism (§4.3).

### 4.2 What gets exchanged

For each pair `(X ∈ Cond1, Y ∈ Cond2)`:

- For every **Negative** protein, intensity values measured in Run X and Run Y
  are exchanged — equivalently, the run identifier on each Negative
  measurement is replaced with its partner's.
- For every **Positive** protein, nothing changes.

### 4.3 Per-pair vs condition-level

Comparing the **two runs in a single pair** is misleading. Take Pair 1
(Seq2 ↔ Seq4) after the swap:

| Class | Seq2 (Cond1) holds | Seq4 (Cond2) holds | Apparent DE Seq2 vs Seq4 |
|---|---|---|---|
| Positive (10%) | original Seq2 (Neat) | original Seq4 (1:2) | ~2× — same as raw |
| Negative (90%) | Seq4's data (1:2) | Seq2's data (Neat) | ~2× — **reversed sign** |

A pairwise comparison would show 100% of proteins as DE. But the benchmark
pools all 20 Cond1 runs against all 20 Cond2 runs. For a Negative protein the
Cond1 bucket is then 10 unpaired runs (true Neat) + 10 paired runs holding
1:2 data; Cond2 mirrors that. Both buckets become 50/50 mixtures with
**identical means** → the protein is null at the condition level. For
Positive proteins, both buckets stay pure → the original ~2× difference
survives.

The 50% swap rate is what makes this work. 0% would leave DE intact for
everyone; 100% would just invert the labels with the same effect size. Only a
partial swap yields a true null.

---

## 5. Variance asymmetry and its impact on moderated methods

The swap creates the desired null *mean* for Negative proteins, but it also
creates unequal variance classes. Positive proteins remain pure dilution
comparisons; Negative proteins become 50/50 mixtures of two dilution levels
inside each condition.

In log2 units, with biological noise variance σ² and a 1-log2 dilution shift:

- **Positive protein:** all runs within a condition have the same dilution
  level, so variance is approximately σ².
- **Negative protein:** each condition contains a 50/50 mixture of Neat and
  1:2 runs, so variance is approximately σ² + 0.25
  (`0.5 * 0.5 * 1²` from the mixture term).

Empirically, the Negative-protein SD mode shifts from ~0.4 to ~0.7 in CSF
(variance 0.16 to 0.49). Thus the swap produces the intended null mean by
making the null class more heterogeneous.

This matters because limma, DEqMS, msqrob2 and limpa moderate the per-protein
variance, whereas MSstats and MSstats+ do not. Moderation shrinks each
protein's empirical variance `s_g²` toward a shared empirical-Bayes prior
`s_0²`. In this benchmark the prior is estimated from a deliberately
imbalanced mixture: ~90% high-variance Negative proteins and ~10% low-variance
Positive proteins.

```
s_0²  ≈  0.9 * (σ² + 0.25) + 0.1 * σ²  =  σ² + 0.225
```

The prior is therefore close to the Negative-protein variance and far above
the Positive-protein variance. Under the simple model, the prior is 0.225
above the Positive variance but only 0.025 below the Negative variance. The
minority Positive class is therefore pulled upward more strongly than the
majority Negative class is pulled downward.

| Class | Share | Mean effect | True variance | Moderation effect | Benchmark consequence |
|---|---:|---|---|---|---|
| Positive | ~10% | ~1 log2 | low | variance pulled up strongly | lower t-statistics, fewer TP, lower TPR |
| Negative | ~90% | ~0 log2 | high | variance pulled down weakly | possible FP if nulls cross cutoff |

For the current Spectronaut log2 table, the second effect is not large enough
to make limpa, limma or msqrob2 call Negative proteins significant at
`p < 0.05`; their PPV is therefore 1.000. This does not remove the design
concern, and we do not interpret this table as evidence of a PPV penalty for
those methods. The observed bias appears mainly as reduced sensitivity,
because the true Positive proteins pay the larger variance-shrinkage penalty.

The root cause is a violation of the exchangeability assumption behind a
single shared variance prior. The benchmark constructs a bimodal variance
distribution tied to the ground-truth label: low variance for Positives, high
variance for Negatives. A one-prior moderation model then shrinks the two
ground-truth classes toward the wrong common target. The resulting performance
gap is therefore partly a property of the benchmark design, not only of method
quality.

---

## 6. Cleaner alternative: pre-swap the source file

Today the swap is re-implemented inside each method's processing path (four
helpers in `benchmark_experiments_functions.R` plus the post-hoc
`ProteinLevelData` swap for MSstats/MSstats+). A cleaner design: build a
single `CSF_Spectronaut_swapped.tsv` once — for every Negative-protein row,
replace `R.FileName` with its partner's (`CSF_annotation_with_Swap.csv`
already has the pairing) — and have every method run its standard converter
on that file with no further manipulation. This kills the stage-of-swap
question (§7), removes four near-duplicate helpers, and makes the ground
truth inspectable on disk.

---

## 7. Where the swap is currently applied

| Method | Stage |
|---|---|
| msqrob2, MaxLFQ + limma, limpa, DEqMS | Raw evidence, *before* protein-level summarization |
| MSstats, MSstats+ (original) | `ProteinLevelData`, *after* `dataProcess` |
| MSstats, MSstats+ (pre-swap) — new | Long-format input, *before* `dataProcess` |

Pre-swap variants were added to test whether stage matters. Result on
Spectronaut at p < 0.05:

| Method | TPR | PPV |
|---|---|---|
| MSstats+ | 0.983 | 0.994 |
| MSstats+ (pre-swap) | 0.983 | 0.994 |
| MSstats | 0.777 | 1.000 |
| MSstats (pre-swap) | 0.777 | 1.000 |

Identical to three decimals. The post- vs pre-summarization choice is not a
confounder for this dataset.

---

## 8. How TPR/PPV are computed

After each method runs on the swapped dataset, every protein gets:

```r
Label = ifelse(Protein %in% true_positives, "Positive", "Negative")
```

[`CSF_Spectronaut_comparison_table.R`](CSF_Spectronaut_comparison_table.R)
then computes:

- `TPR = #(p < 0.05 ∧ finite logFC ∧ Label="Positive")
       / #(finite logFC ∧ Label="Positive")`
- `PPV = 1 − #(p < 0.05 ∧ Label="Negative") / #(p < 0.05)`

A method that recovers the dilution signal without over-calling the null
class scores high on both.

---

## 9. Caveats

1. **Single deterministic realization.** Pairing is fixed by
   `seq(2, n, by = 2)` and the Positive/Negative split is fixed in the swap
   list. No replicates, no confidence intervals. Two methods within ~1% of
   each other are indistinguishable from design noise. Re-pairing or
   re-sampling the protein split across seeds would give uncertainty bands.
2. **Origin of the 182 Positive accessions is not in this repo.** The swap
   list is checked in as data; no script regenerates it.
3. **TPR denominator drops non-finite logFC.** Proteins where one condition
   is entirely missing are excluded, inflating TPR vs. using the full 182.
4. **Raw `p < 0.05`, no FDR.** Apples-to-apples across methods, but worth
   stating in the manuscript.
5. **Bad-pair contamination (§4.1, pairs 8-10).** Three of ten pairs swap
   between mismatched dilution factors (e.g. 1:8 ↔ 1:16). For Negative
   proteins these contribute extra heterogeneity to the null class — a minor
   but real source of variance inflation on top of §5.
