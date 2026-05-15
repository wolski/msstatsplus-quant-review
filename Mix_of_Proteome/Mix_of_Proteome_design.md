# Mix of Proteome benchmark: dataset design and ground-truth construction

This document explains how the Mix_of_Proteome (Olsen Astral, 2025-04-22
re-analysis of a 2023 acquisition) benchmark builds its ground truth of
differentially abundant ("Positive") and non-differentially abundant
("Negative") proteins. Companion notebook:
[`Mix_of_Proteome_visualization.qmd`](Mix_of_Proteome_visualization.qmd).

The key contrast with the CSF and K562 dilution-series benchmarks is that
**Mix_of_Proteome does not need a swap**. The ground truth is built into the
mixing design: every comparison contains a constant-abundance organism (Human)
that defines the null class and two organisms (E. coli, Yeast) whose
abundances change by a known ratio between conditions and define the positive
class.

---

## 1. The raw experiment

A single 200 ng tryptic peptide mix is acquired in DIA on a single LC-MS
platform (Astral, 30 min, 2 Th × 3.5 ms, 180 K). The mix is composed of three
organism digests in deliberately varied proportions
([`Mix_of_Proteome_annotation.csv`](Mix_of_Proteome_annotation.csv)):

| Condition   | E. coli (%) | Human (%) | Yeast (%) | BioReplicates | Runs (Order) |
|---|---|---|---|---|---|
| `E5H50Y45`  |  5 | 50 | 45 | 3 | 3, 4, 5 |
| `E10H50Y40` | 10 | 50 | 40 | 3 | 1, 2, 7 |
| `E20H50Y30` | 20 | 50 | 30 | 3 | 8, 9, 10 |
| `E30H50Y20` | 30 | 50 | 20 | 3 | 11, 12, 13 |
| `E40H50Y10` | 40 | 50 | 10 | 3 | 14, 15, 16 |
| `E45H50Y5`  | 45 | 50 |  5 | 3 | 6, 17, 18 |

**18 runs total, 6 conditions × 3 replicates, 15 pairwise comparisons.** The
Human fraction is held at 50% in every condition; only E. coli and Yeast vary.

---

## 2. Why no swap is needed

CSF/K562 are dilution series of a single proteome, so every detected protein
differs by the dilution factor between conditions and there is no intrinsic
null class. They have to manufacture one with a partial run/label swap.

Mix_of_Proteome carries the null and positive classes physically:

- **Human (50% in every condition)** → expected log2FC = 0 in every pairwise
  comparison → **null class** (true negatives).
- **E. coli (5–45%)** → expected log2FC = `log2(c1 / c2)` with the
  Condition-coded percentage → **positive class** (true positives).
- **Yeast (45–5%)** → mirrors E. coli inversely → **positive class**.

Organism membership is the ground truth, looked up from
[`idmapping.tsv`](idmapping.tsv). No swap, no protein swap list, no run
re-labelling.

---

## 3. Expected log2 fold changes per comparison

For each ordered pairwise contrast `c1 vs c2`, the expected log2FC is `0` for
Human, `log2(E_c1 / E_c2)` for E. coli, and `log2(Y_c1 / Y_c2)` for Yeast
([`Mixture_of_proteomes_processing.R:947-1017`](Mixture_of_proteomes_processing.R#L947-L1017)).
Numerically:

| Contrast (c1 vs c2)       | Human | E. coli | Yeast |
|---|---:|---:|---:|
| E10H50Y40 vs E20H50Y30 | 0 | log2(10/20) = -1.000 | log2(40/30) =  0.415 |
| E10H50Y40 vs E30H50Y20 | 0 | log2(10/30) = -1.585 | log2(40/20) =  1.000 |
| E10H50Y40 vs E40H50Y10 | 0 | log2(10/40) = -2.000 | log2(40/10) =  2.000 |
| E10H50Y40 vs E45H50Y5  | 0 | log2(10/45) = -2.170 | log2(40/5)  =  3.000 |
| E10H50Y40 vs E5H50Y45  | 0 | log2(10/5)  =  1.000 | log2(40/45) = -0.170 |
| E20H50Y30 vs E30H50Y20 | 0 | log2(20/30) = -0.585 | log2(30/20) =  0.585 |
| E20H50Y30 vs E40H50Y10 | 0 | log2(20/40) = -1.000 | log2(30/10) =  1.585 |
| E20H50Y30 vs E45H50Y5  | 0 | log2(20/45) = -1.170 | log2(30/5)  =  2.585 |
| E20H50Y30 vs E5H50Y45  | 0 | log2(20/5)  =  2.000 | log2(30/45) = -0.585 |
| E30H50Y20 vs E40H50Y10 | 0 | log2(30/40) = -0.415 | log2(20/10) =  1.000 |
| E30H50Y20 vs E45H50Y5  | 0 | log2(30/45) = -0.585 | log2(20/5)  =  2.000 |
| E30H50Y20 vs E5H50Y45  | 0 | log2(30/5)  =  2.585 | log2(20/45) = -1.170 |
| E40H50Y10 vs E45H50Y5  | 0 | log2(40/45) = -0.170 | log2(10/5)  =  1.000 |
| E40H50Y10 vs E5H50Y45  | 0 | log2(40/5)  =  3.000 | log2(10/45) = -2.170 |
| E45H50Y5  vs E5H50Y45  | 0 | log2(45/5)  =  3.170 | log2(5/45)  = -3.170 |

Effect-size magnitudes range from **0.17 to 3.17 log2 units**. Comparisons
between adjacent conditions (e.g. E10 vs E20) are hard (small fold change);
extreme comparisons (e.g. E5 vs E45) are easy.

---

## 4. How TPR / FDR are computed

After each method runs on the full dataset, every protein is mapped to its
organism via `idmapping.tsv` and labelled:

```r
truth = ifelse(Organism == "Homo sapiens (Human)", "Negative", "Positive")
```

[`Mixture_of_proteomes_processing.R:669-731`](Mixture_of_proteomes_processing.R#L669-L731)
defines:

- `FDR  = #(adj.p < 0.05 ∧ Human) / #(adj.p < 0.05)`     (per contrast, then pooled)
- `TPR  = #(adj.p < 0.05 ∧ non-Human) / #(non-Human)`    (per contrast, then pooled)

i.e. *any* significant Human protein is a false positive, and *any* significant
E. coli or Yeast protein is a true positive, regardless of effect-size match.
The `error = log2FC − expected` column (lines 1072-1076) is reserved for
calibration plots, not for FDR scoring.

ROC/AUROC is computed against the same binary truth at sliding `adj.p`
thresholds (lines 1118-1218).

---

## 5. Variance considerations — why moderation is not penalised here

Empirical-Bayes variance moderation (limma, DEqMS, msqrob2 internally,
limpa/dpc) shrinks per-protein variances toward a common prior fit across the
protein set. In CSF/K562 the swap creates **bimodal** within-condition
variance (low for Positives, high for Negatives) which biases the prior; here
the situation is different and largely benign:

- All 18 runs are biological replicates of the *same* total digest mass on the
  *same* instrument. Within-condition variance is dominated by injection /
  measurement noise, which has no reason to differ between Human and
  E. coli or Yeast peptides at matched intensity.
- The mixture composition changes the **mean** abundance per organism per
  condition (driving the fold-change signal), but the **within-condition
  spread** at a given mean is set by the ionisation/quantitation noise floor,
  not by which organism the peptide came from.
- The empirical-Bayes prior therefore sees an approximately unimodal variance
  distribution. Moderation does what it is meant to do — stabilise small-N
  variance estimates — without a class-asymmetric bias.

Caveat: at the **low end** of each organism's dynamic range (E. coli at 5 %,
Yeast at 5 %) abundance approaches the limit of quantification and missingness
rises. That inflates the per-protein SD for those organisms in those
conditions — but only at the low-abundance tail, not as a class-wide shift.
Section 4 of the companion notebook visualises this directly with
mean-variance trends per organism.

The take-away: unlike CSF/K562, **Mix_of_Proteome is a benchmark on which
variance moderation is not at a structural disadvantage**. Differences between
moderating and non-moderating methods on this dataset reflect the methods,
not the design.

---

## 6. (No alternative pre-swap variant — N/A)

There is no swap, so the §6/§7 of the CSF/K562 docs has no analogue.

---

## 7. Caveats

1. **Constant-Human assumption.** The benchmark assumes Human peptide signal
   is identical across conditions. In a digest where total peptide mass is
   constant but the *non*-Human fractions change, the *amount* of Human
   peptide loaded is constant in absolute terms only if the input was mixed
   on a mass basis — the protocol assumes this. Any normalisation choice
   (none, median, equalised-Human) interacts with the FDR scoring of Human
   proteins. Two of the methods in [Mixture_of_proteomes_processing.R](Mixture_of_proteomes_processing.R)
   run with `normalization = FALSE` (MSstats, MSstats+, lines 56, 85), the
   others use their defaults; this is a known confounder and is documented in
   the comparison table.

2. **No replicates of the design.** Three BioReplicates per condition is the
   only randomness. Two methods within ~1 % AUROC of each other are within
   sampling noise.

3. **Adjacent-condition comparisons are hard.** E10 vs E20 has E. coli log2FC
   of -1 but Yeast log2FC of only +0.41. Methods can score well on extreme
   contrasts (E5 vs E45) and poorly on near contrasts; reporting per-contrast
   AUROC matters more here than on CSF/K562 where every comparison is the
   same canonical Neat-vs-1:2.

4. **Missingness is class-asymmetric at the extremes.** At E. coli 5% and
   Yeast 5%, peptide drop-out is higher than at 50% Human. Methods that
   discard rows with non-finite logFC (e.g. the limma table loop in
   [Mixture_of_proteomes_processing.R:274-285](Mixture_of_proteomes_processing.R#L274-L285))
   silently bias the TPR denominator toward proteins detectable in both
   conditions. This is consistent across methods but inflates absolute TPR.

5. **Three-organism FDR pooling.** `FDR = FP_Human / discoveries` pools
   E. coli and Yeast hits in the denominator. A method biased toward calling
   only E. coli (or only Yeast) gets the same FDR as one calling both, even
   though calibration on the un-called organism may be poor. Per-organism
   AUROC and the calibration plot (`error = log2FC − expected` per organism)
   are the better diagnostics.

---

## 8. Summary

| Aspect | CSF / K562 | Mix_of_Proteome |
|---|---|---|
| Source material | One proteome, dilution series | Three proteomes, varied proportions |
| Null class | Manufactured by 50% run swap on 90 % of proteins | Built-in: Human (50 % in every condition) |
| Positive class | Manufactured by leaving 10 % of proteins un-swapped | Built-in: E. coli + Yeast |
| Effect size | One value (~1 log2) for every Positive | Known per organism, per comparison; range 0.17–3.17 log2 |
| Variance moderation penalty | Yes — bimodal variance from the swap | No — approximately unimodal |
| Number of comparisons | 1 canonical (Cond1 vs Cond2) | 15 pairwise |
