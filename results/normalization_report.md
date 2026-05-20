# Normalization across the 7 modelling packages — short report

Short summary of how `log2`, `median`, and `quantile` are realised in each package's current implementation. Full audit with call-site references: [../R/README_normalization.md](../R/README_normalization.md).

## Headline

| | log2 | median | quantile |
|---|---|---|---|
| What it does | log2 only, no inter-sample normalization | **per-sample median centering** (centering only — no scaling) | quantile normalization of the per-sample distributions |
| Operation on log2 scale | `M` | `M − colMedian(M)` (NA-safe) | `limma::normalizeBetweenArrays(method = "quantile")` |

All three labels are **centering or distribution-matching only**. None of the seven packages applies per-sample scaling (MAD / IQR / SD division) under the `median` label. Verified directly from `MSstats:::.normalizeMedian` source (v4.18.1): `ABUNDANCE := ABUNDANCE − ABUNDANCE_RUN + ABUNDANCE_FRACTION` — a pure additive shift.

## What differs between packages

Not the arithmetic, but the **level** at which normalization is applied:

| Package | Level of normalization |
|---|---|
| MSstats, MSstats+ | feature / precursor, grouped by RUN × FRACTION (inside `dataProcess`) |
| MaxLFQ + limma | protein matrix, after MaxLFQ summarization |
| DEqMS | protein matrix, after summarization |
| msqrob2 | protein matrix; **plus** an upstream peptide-level `center.median` via `QFeatures::normalize` that runs regardless of label |
| prolfqua | protein matrix, after medpolish summarization |
| limpa | protein matrix |

Implication for the benchmark: the `median` label is comparable across packages on the arithmetic axis (centering only, no hidden scaling). It is **not** comparable on the level axis — MSstats normalizes feature intensities before protein rollup, the others normalize the protein matrix after rollup. Either approach is defensible; we keep each package's native level rather than forcing one or the other.

## What we will fix for the shared `quant/R/normalize.R`

- **`log2`**: log2 only, no normalization.
- **`median`**: NA-safe `sweep(m, 2, colMedians(m, na.rm = TRUE), "-")` on log2 protein matrix; for MSstats / MSstats+, use `dataProcess(normalization = "equalizeMedians")`.
- **`quantile`**: `limma::normalizeBetweenArrays(method = "quantile")` on log2 protein matrix; for MSstats / MSstats+, use `dataProcess(normalization = "quantile")`.

## Open question (deferred)

Whether to add a fourth label `median_scaled` (centering + per-sample MAD scaling) as a sensitivity analysis for the variance-shrinkage story. Default: do not expose. Revisit after the first end-to-end results land.
