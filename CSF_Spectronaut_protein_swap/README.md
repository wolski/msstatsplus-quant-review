# CSF Spectronaut protein-swap benchmark

This folder contains the protein-swap benchmark derived from the CSF
Spectronaut report. The benchmark swaps precursor intensities rank-for-rank
between matched protein pairs in a randomly selected `G2` half of the runs, then
tests methods on the `G2 - G1` contrast.

Run commands from the repository root (`quant/`), not from this folder.

```bash
make prep-protein-swap
make cells-protein-swap
make diag-protein-swap
```

`make prep-protein-swap` generates:

- `CSF_Spectronaut_protein_swap/Report.tsv`
- `CSF_Spectronaut_protein_swap/annotation.csv`
- `CSF_Spectronaut_protein_swap/*_swap_ground_truth.tsv`
- `CSF_Spectronaut_protein_swap/*_swap_true_positives.tsv`
- `CSF_Spectronaut_protein_swap/*_swap_group_annotation.csv`
- subset directories: `all_data`, `good_data`, `small`, `small_good_data`

`make cells-protein-swap` runs the shared modelling pipeline in `../R/` across
the configured subset, normalization, and method grid. Outputs are written under:

```text
CSF_Spectronaut_protein_swap/<subset>/<normalization>/swap/<method>/
```

The folder-specific Makefile is included by the top-level `Makefile`, but can
also be invoked directly from `quant/`:

```bash
make -f CSF_Spectronaut_protein_swap/Makefile prep-protein-swap
```

Generated reports, subset directories, model outputs, diagnostics, and prep
stamps are not source files. Recreate them with the Makefile targets above.
