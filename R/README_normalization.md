# Normalization audit (§0.5)

Goal: settle what `log2`, `median`, and `quantile` actually mean across the 7 modelling packages used in the benchmark, so the shared `quant/R/normalize.R` (built in §1) can fix one definition per label.

Scope of this document: how each package currently normalizes in the existing scripts. The conclusion drives `quant/R/normalize.R`; this file does not modify any modelling code.

## TL;DR

- **Centering vs scaling**: all current "median" paths are **centering only** (subtract per-sample/per-run median). None of the 7 packages applies scaling (MAD / IQR / SD division) under the `median` label. The Explore-agent claim that MSstats `equalizeMedians` does centering + scaling is wrong; see source quote below.
- **Real inconsistency is the *level* at which normalization is applied**:
  - **MSstats / MSstats+** normalize at the **feature (precursor) level**, grouped by `RUN × FRACTION`, *before* protein summarization.
  - **MaxLFQ+limma, DEqMS, msqrob2, prolfqua, limpa** normalize the **protein-level matrix** *after* summarization (msqrob2 also does an earlier `center.median` at the peptide level via `QFeatures::normalize`).
- **Quantile**: MSstats uses its own `.normalizeQuantile` (operates on the feature × run wide table). All others use `limma::normalizeBetweenArrays(method="quantile")` on the protein-level log2 matrix.
- **Log2**: in all 7 packages, "log2" means "no inter-sample normalization, abundance on log2 scale". The log2 transform itself is applied internally by each package.

## Per-package call sites (current code)

| Package | Norm | Call site | Underlying function | What it does | Level | NA handling |
|---|---|---|---|---|---|---|
| MSstats | log2 | [run_msstats.R:24](../CSF_Spectronaut_swap/run_msstats.R#L24), [:56](../CSF_Spectronaut_swap/run_msstats.R#L56), [:100](../CSF_Spectronaut_swap/run_msstats.R#L100) | `dataProcess(normalization = FALSE)` | log2 only, no normalization | feature | MSstats internal |
| MSstats | median | [run_msstats.R:25](../CSF_Spectronaut_swap/run_msstats.R#L25), [:56](../CSF_Spectronaut_swap/run_msstats.R#L56), [:100](../CSF_Spectronaut_swap/run_msstats.R#L100) | `dataProcess(normalization = "equalizeMedians")` → `.normalizeMedian` | **centering only** (per-run median shift to common target) | feature, per RUN × FRACTION | NA-safe (`median(..., na.rm=TRUE)`) |
| MSstats | quantile | [run_msstats.R:11](../CSF_Spectronaut_swap/run_msstats.R#L11) | `dataProcess(normalization = "quantile")` → `.normalizeQuantile` | quantile normalization via wide-table per fraction | feature × run | NA-safe |
| MSstats+ | log2 / median / quantile | same as MSstats | same as MSstats | same as MSstats | same as MSstats | same as MSstats |
| MaxLFQ+limma | log2 | [run_nonmsstats.R:155-158](../CSF_Spectronaut_swap/run_nonmsstats.R#L155-L158) | `iq::maxLFQ` (already log2) | none | protein matrix | n/a |
| MaxLFQ+limma | median | [run_nonmsstats.R:163](../CSF_Spectronaut_swap/run_nonmsstats.R#L163) | `median_normalize_log2_matrix` | **centering only** (`sweep(m, 2, colMedians(m), "-")`) | protein matrix | NA-safe (`colMedians(..., na.rm=TRUE)`) |
| MaxLFQ+limma | quantile | [run_nonmsstats.R:161](../CSF_Spectronaut_swap/run_nonmsstats.R#L161) | `quantile_normalize_log2_matrix` → `limma::normalizeBetweenArrays(method="quantile")` | quantile | protein matrix | limma rank-based, NA-safe |
| DEqMS | log2 | implicit (matrix already log2) | none | none | protein matrix | n/a |
| DEqMS | median | [run_nonmsstats.R:291](../CSF_Spectronaut_swap/run_nonmsstats.R#L291) | `median_normalize_log2_matrix` | **centering only** | protein matrix | NA-safe |
| DEqMS | quantile | [run_nonmsstats.R:287](../CSF_Spectronaut_swap/run_nonmsstats.R#L287) | `quantile_normalize_log2_matrix` | quantile via limma | protein matrix | NA-safe |
| msqrob2 | log2 | [run_nonmsstats.R:69](../CSF_Spectronaut_swap/run_nonmsstats.R#L69) | `QFeatures::normalize(method="center.median")` at peptide, then summarize | peptide-level **centering only** | peptide → protein matrix | QFeatures internal |
| msqrob2 | median | [run_nonmsstats.R:86](../CSF_Spectronaut_swap/run_nonmsstats.R#L86) (after peptide center.median + protein summary) | `median_normalize_log2_matrix` | **centering only** (on top of peptide-level center.median) | protein matrix | NA-safe |
| msqrob2 | quantile | [run_nonmsstats.R:84](../CSF_Spectronaut_swap/run_nonmsstats.R#L84) | `quantile_normalize_log2_matrix` | quantile (on top of peptide-level center.median) | protein matrix | NA-safe |
| prolfqua | log2 | [run_prolfqua_step.R:82](../CSF_Spectronaut_swap/run_prolfqua_step.R#L82) | `tr_norm$log2()` only | none | protein matrix | prolfqua internal |
| prolfqua | median | [run_prolfqua_step.R:76-79](../CSF_Spectronaut_swap/run_prolfqua_step.R#L76-L79) | `tr_norm$log2()` then `intensity_matrix(.func = median_normalize_log2_matrix)` | **centering only** | protein matrix | NA-safe |
| prolfqua | quantile | [run_prolfqua_step.R:68-75](../CSF_Spectronaut_swap/run_prolfqua_step.R#L68-L75) | `tr_norm$log2()` then `intensity_matrix(.func = quantile_normalize_log2_matrix)` | quantile via limma | protein matrix | NA-safe |
| limpa | log2 | [run_nonmsstats.R:200-207](../CSF_Spectronaut_swap/run_nonmsstats.R#L200-L207) | matrix log2 only | none | protein matrix | n/a |
| limpa | median | [run_nonmsstats.R:211](../CSF_Spectronaut_swap/run_nonmsstats.R#L211) | `median_normalize_log2_matrix(log2(...))` | **centering only** | protein matrix | NA-safe |
| limpa | quantile | [run_nonmsstats.R:207](../CSF_Spectronaut_swap/run_nonmsstats.R#L207) | `quantile_normalize_log2_matrix(log2(...))` | quantile via limma | protein matrix | NA-safe |

## Reference: MSstats `equalizeMedians` source

Quoted from `MSstats` 4.18.1 (installed at `~/Library/R/4.5-arm64/MSstats`), dumped via `Rscript -e 'print(asNamespace("MSstats")$.normalizeMedian)'`:

```r
.normalizeMedian <- function (input) {
    if (length(unique(input$LABEL)) == 1L) label <- "L" else label <- "H"
    input[, ABUNDANCE_RUN := .getMedian(.SD, label),
          by = c("RUN", "FRACTION"),
          .SDcols = c("ABUNDANCE", "LABEL")]
    input[, ABUNDANCE_FRACTION := median(ABUNDANCE_RUN, na.rm = TRUE),
          by = "FRACTION"]
    input[, ABUNDANCE := ABUNDANCE - ABUNDANCE_RUN + ABUNDANCE_FRACTION]
    input <- input[, !(colnames(input) %in% c("ABUNDANCE_RUN", "ABUNDANCE_FRACTION")),
                   with = FALSE]
    input
}
```

Reading: `ABUNDANCE` is already log2 (set upstream by `dataProcess`). For each `(RUN, FRACTION)` group, compute the run-level median of the reference-label feature abundances. For each `FRACTION`, take the median of those run-medians as the common target. Shift each run so its median lands on that target: `ABUNDANCE := ABUNDANCE − ABUNDANCE_RUN + ABUNDANCE_FRACTION`.

This is **centering only** — a pure additive shift on the log2 scale, no scaling factor. The in-repo comment at [run_step_common.R:73](../CSF_Spectronaut_swap/run_step_common.R#L73) ("Equivalent in spirit to MSstats's equalizeMedians") is correct on the arithmetic axis; the only deviation is the **level** at which the operation is applied.

## Decision for the shared `quant/R/normalize.R`

Adopt the following canonical definitions for the three labels used across `<csffolder>/<dataset>/<normalization>/swap/<modellingPackage>`:

- **`log2`**: log2 only, no inter-sample normalization.
- **`median`**: per-sample median centering on a log2-scale matrix, NA-safe (`sweep(m, 2, colMedians(m, na.rm=TRUE), "-")`).
- **`quantile`**: `limma::normalizeBetweenArrays(method = "quantile")` on a log2-scale matrix.

Level: all three are applied at the **protein-level matrix** for the non-MSstats packages (where they already are today). For MSstats / MSstats+, we keep using `dataProcess`'s built-in routines because (a) they operate at the feature level by design and the package's downstream model assumes that, and (b) the arithmetic matches (centering only / quantile). The divergence is in the level, not the operation.

Per-package adjustments to enforce these definitions:

| Package | `log2` | `median` | `quantile` |
|---|---|---|---|
| MSstats / MSstats+ | `dataProcess(normalization = FALSE)` | `dataProcess(normalization = "equalizeMedians")` | `dataProcess(normalization = "quantile")` |
| MaxLFQ+limma, DEqMS, msqrob2, prolfqua, limpa | matrix log2 only | `median_normalize_log2_matrix(protein_log2_matrix)` | `quantile_normalize_log2_matrix(protein_log2_matrix)` |

Caveat (document, do not "fix"): msqrob2 currently does a peptide-level `QFeatures::normalize(method = "center.median")` *before* protein summarization. This is an upstream centering step that happens regardless of the chosen normalization label. Keep it for `msqrob2` (it is part of msqrob2's documented pipeline); note it in the per-package row in the manuscript table.

## Open question (not exposed yet)

Should we add a fourth label, `median_scaled`, for centering + per-sample MAD scaling? Useful as a sensitivity analysis if the variance-shrinkage story in §7.5 turns out to depend on whether per-sample spread is equalized. Default: do not expose. Re-evaluate after the first end-to-end results of §5 / §6 are in.

## Empirical smoke check (optional, for re-verification)

If the canonical definitions above are ever in doubt, regenerate this block by running on any protein-level log2 matrix `M`:

```r
m1 <- median_normalize_log2_matrix(M)
m2 <- quantile_normalize_log2_matrix(M)
# Both should leave protein-level structure interpretable on log2 scale.
all.equal(matrixStats::colMedians(m1, na.rm=TRUE),
          rep(0, ncol(m1)))            # TRUE: medians shifted to zero
stopifnot(diff(range(matrixStats::colMedians(m2, na.rm=TRUE))) < 1e-8)
# Per-sample SDs should NOT be equalized by either method:
matrixStats::colSds(m1, na.rm=TRUE)    # expect heterogeneous
matrixStats::colSds(m2, na.rm=TRUE)    # expect heterogeneous (quantile equalizes the full distribution, not the SD per se)
```

If `colSds(m1, ...)` ever comes back equal across columns under the shared helper, the helper has been changed to scaling — investigate before proceeding with §1.
