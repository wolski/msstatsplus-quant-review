# TODO Review 2: Assessment And Rewrite Plan

This document records the assessment of `vignettes/review2.qmd` and the
recommended next steps for turning it into an effective reviewer response.

## Executive Recommendation

Use two documents:

1. A short main review, about 1-2 pages, written for the editor and authors.
2. A longer reproducible supplement based on the current `review2.qmd` /
   `review2.pdf`, containing the tables, diagnostics, sensitivity analyses, and
   output inventory.

COMMENT: leave review2.qmd untouched, copy the relevant (all figures tables figure discussions, technical conclusions) to  review_supplement.qmd


The current review is scientifically serious and mostly well grounded, but it
is too close to a parallel research article. The strongest path is to make the
main review narrow, procedural, and hard to dismiss, while keeping the full
reanalysis available as supporting material.

## Core Finding To Preserve

The central critique is not that MSstats+ is implemented incorrectly. The
central critique is that the K562/CSF sample-swap benchmark construction can
couple ground-truth class to within-condition variance.

The benchmark starts from a dilution experiment where essentially all proteins
are differentially abundant. The manuscript creates a null class by reassigning
condition labels for most proteins. In the CSF sample-swap construction,
Positive proteins retain the cleaner dilution comparison, while Negative
proteins become within-condition mixtures of dilution levels.

This matters for methods that use empirical-Bayes or related variance
moderation. A shared variance prior estimated from a class-imbalanced mixture
can be influenced by the high-variance Negative class and can reduce
sensitivity for the lower-variance Positive class. Therefore Table 1 supports a
narrower conclusion than the manuscript currently suggests: MSstats+ performs
well on this specific sample-swap stress-test construction, but Table 1 alone
does not establish broad superiority over established differential-abundance
workflows.

Recommended short wording:

> My concern is not that MSstats+ is implemented incorrectly. My concern is
> that the sample-swap benchmark induces a label-dependent variance structure:
> the synthetic null proteins are constructed as within-condition mixtures of
> dilution levels, whereas the positive proteins retain the original dilution
> comparison. This makes the benchmark partly a test of how each method responds
> to this artificial variance asymmetry.

And:

> Therefore, Table 1 supports a narrower conclusion than currently stated:
> MSstats+ performs well on this specific sample-swap benchmark with induced
> anomalous runs. It does not by itself establish general superiority over
> variance-moderated workflows in standard differential-abundance settings.

## What The Main Review Should Contain

The main review should be short and decision-relevant:

1. Acknowledge that MSstats+ may be valuable and that the benchmark is useful as
   a quality-degradation stress test.
2. State that the benchmark construction, not the software implementation, is
   the main concern.
3. Explain that the benchmark has no natural null because all proteins are
   differential before label reassignment.
4. Explain that label reassignment creates Negative proteins that are
   within-condition mixtures of dilution levels, while Positives keep the
   original dilution effect.
5. Explain why this can disadvantage moderated methods and complicate the
   interpretation of Table 1.
6. Ask the authors to narrow the benchmark claims or add sensitivity analyses
   where truth labels are not coupled to variance structure.
7. Ask that nominal p-value and adjusted-p-value operating points be clearly
   separated in the main benchmark interpretation.
8. Ask that benchmark generation be separated from tool-specific analysis code,
   ideally by writing a fixed benchmark report, annotation, and truth table to
   disk before running methods.

## What Should Move To The Supplement

Keep these in the long supplement, not in the main review:

- Full reproduced CSF Spectronaut benchmark tables.
- P-value distribution panels.
- Per-protein effect size and within-condition variance figures.
- Rebuilt sample-swap results.
- Protein-swap benchmark construction and result tables.
- Output availability and incomplete-cell inventory.
- Runtime, memory, and prolfqua inclusion comments.
- Any rows where pipelines differ from the manuscript, especially `msqrob2`.

The protein-swap and rebuilt sample-swap analyses should be framed as
sensitivity analyses, not as definitive replacement benchmarks. Their purpose is
to show that method rankings are sensitive to benchmark construction.

Recommended wording:

> These analyses are not proposed as definitive alternative benchmarks. They
> are sensitivity checks showing that the performance ranking is not stable once
> the sample-swap-induced H0/H1 variance asymmetry is removed or made explicit.

## Risks In The Current Review

The current `review2.qmd` has several vulnerabilities:

- It is long enough to read like a research article rather than a peer review.
- It includes incomplete output cells, which can weaken the final message.
- It relies heavily on CSF Spectronaut; the manuscript also includes K562,
  DIA-NN, mixture-of-proteomes, and clinical examples.
- The protein-swap benchmark changes the data-generating mechanism and protein
  universe, so authors can argue that it tests a different question.
- Some `msqrob2` rows are not directly comparable because the rerun uses a
  different or faster pipeline.
- Strong causal statements about moderation need either direct evidence or
  softer wording.

## Tone And Wording Changes

Avoid:

- "biased benchmark" as the lead phrase.
- "the authors are wrong."
- "the benchmark is invalid."
- "MSstats+ does not outperform other methods."
- "the protein-swap benchmark proves the original conclusion false."

Prefer:

- "the benchmark construction induces label-dependent variance structure."
- "this complicates interpretation of Table 1."
- "the evidence supports performance on this sample-swap construction, but not
  the broader claim without additional sensitivity analyses."
- "method rankings are sensitive to benchmark construction."

Specific edits to consider:

- Replace "This is the correct object for criticizing the manuscript benchmark"
  with "This object matches the authors' benchmark code path and is therefore
  the relevant object for evaluating the Table 1 CSF result."
- Replace "The two trends show the moderation problem directly" with "The two
  trends show the variance structure that can disadvantage methods using a
  shared moderation trend."
- Replace "the apparent advantage of MSstats+ is not a general property" with
  "the CSF Table 1 advantage alone does not establish a general property."

## Likely Author Pushback And Preemption

Likely pushback:

- The benchmark intentionally models longitudinal quality degradation.
- Permutation preserves run-level quality distributions and is a standard way
  to create negatives.
- MSstats+ also performs well in the mixture-of-proteomes benchmark.
- The protein-swap benchmark is artificial and tests a different question.
- Some reanalysis pipelines differ from the authors' submitted pipelines.
- Nominal p-values were intentionally used and adjusted results are in the
  supplement.

Preemptive framing:

- Agree that the K562/CSF benchmark is valuable as a quality-degradation stress
  test.
- State that preserving run-level structure is not the issue; the issue is that
  in a dilution series, the permutation also creates class-dependent
  within-condition mixture variance.
- Acknowledge that the mixture benchmark supports the usefulness of MSstats+,
  while maintaining that Table 1 should not be overgeneralized.
- Present protein-swap as a diagnostic sensitivity analysis, not as a required
  replacement benchmark.
- Keep non-identical pipelines and incomplete results out of the main argument.
- Accept nominal p-values as a diagnostic, but ask that discovery-oriented
  claims also foreground adjusted-p-value results.

## Additional Evidence That Would Strengthen The Review

Highest-value additions:

- Quantitative summaries of within-condition SD by class before and after swap:
  median, IQR, effect size, and possibly bootstrap intervals.
- The same variance diagnostic for K562, because the manuscript's benchmark
  claim includes K562 as well as CSF.
- Direct method-level evidence for moderation effects: raw variance, moderated
  variance, t-statistic, and p-value by class where available.
- A minimal reproducible script showing the exact sample-swap operation and the
  resulting class-specific variance shift from the authors' submitted files.
- A clean sensitivity comparison where preprocessing, normalization, protein
  universe, and method pipelines are as aligned as possible.

## Credit And Acknowledgement

Do not ask the authors directly for acknowledgement or authorship as part of the
peer review. That would be procedurally awkward and could weaken the review.

If the authors adopt the alternative benchmark, code, generated data, or figures,
then attribution is appropriate. The safer route is to raise this with the
editor, not the authors.

Suggested wording to the editor:

> This review includes a substantial independent reproducible reanalysis. If the
> authors choose to incorporate the alternative benchmark, generated data, code,
> or figures from this review into the manuscript, I ask the journal to advise on
> appropriate citation or acknowledgement.

Do not make acknowledgement a condition of acceptance. Make benchmark validity
and claim scope the condition.

## Created Documents

- `vignettes/review_main.qmd`: short main review, intended for the editor and
  authors.
- `vignettes/review_supplement.qmd`: technical supplement copied from
  `review2.qmd` and reframed as supporting material. The original
  `review2.qmd` was left untouched.

## Action Checklist

- [x] Draft a 1-2 page main review focused on the label-dependent variance
      critique.
- [x] Move full tables, figures, inventories, and sensitivity analyses into a
      supplement.
- [x] Remove incomplete or partial cells from the main argument.
- [x] Mark protein-swap and rebuilt sample-swap as sensitivity analyses.
- [x] Soften causal language unless direct moderation evidence is added.
- [ ] Add quantitative before/after variance summaries if time permits.
- [ ] Add K562 variance diagnostics if time permits.
- [x] Keep `msqrob2` pipeline differences disclosed but noncentral.
- [x] Ask the editor, not the authors, about attribution if the reanalysis is
      reused in the manuscript.
