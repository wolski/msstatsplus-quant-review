# §2 Sample-swap Python script — what was built

Plan: [/Users/wolski/.claude/plans/jolly-petting-giraffe.md](/Users/wolski/.claude/plans/jolly-petting-giraffe.md) — §2 of the parent TODO. Note: that plan file holds only the §1 plan; the §2 work spec is in [TODO/TODO_step_by_step.md](../TODO/TODO_step_by_step.md) §2.

## What changed

- Moved [src/swap_spectronaut_report.py](../src/swap_spectronaut_report.py) (was `quant/swap_spectronaut_report.py`) → `quant/src/`. Updated callers in [CSF_Spectronaut_swap/README.md:16](../CSF_Spectronaut_swap/README.md#L16). The other references (qmd author byline, scaffolding/TODO docs) describe the script and need no path edit.
- Created [src/swap_spectronaut_report_samples.py](../src/swap_spectronaut_report_samples.py) — the new TSV-level sample-swap script.

## Algorithm

Mirrors `swap_condition_labels_msstats()` in [benchmark_experiments_functions.R:63-96](../benchmark_experiments_functions.R#L63-L96), with the cleaner "rewrite on disk" approach from §6 of [CSF_Spectronaut/CSF_swap_design.md](../CSF_Spectronaut/CSF_swap_design.md):

1. Read TSV + annotation + protein swap list.
2. Build the partner mapping: within each of the two conditions (default `Condition1`/`Condition2`), sort runs by `Order`, pick every-second run (R's `seq(2, n, by = 2)`, 0-based indices `1, 3, 5, ...`), pair index-for-index across conditions. Pair count = `min(len(a), len(b))`.
3. For each TSV row whose `PG.ProteinGroups` is a **Negative** protein, replace `R.FileName` with its partner's. Intensities themselves are NOT touched — swapping the run identifier is operationally equivalent, because downstream code joins by `R.FileName` and reads the experimental `Condition` from the annotation.
4. Positive proteins, unpaired runs, and blanks pass through untouched.

## CLI surface

```
python src/swap_spectronaut_report_samples.py \
  --report             <Spectronaut TSV> \
  --annotation         <annotation CSV with R.FileName, Condition, Order> \
  --protein-swap-list  <CSV with Protein, Label columns> \
  --out-dir            <output dir> \
  [--cond-a Condition1] [--cond-b Condition2] [--blank-condition Blank]
```

## Outputs (next to `<stem>`)

| File | Schema |
|---|---|
| `<stem>_sample_swap.tsv` | Schema-identical Spectronaut TSV; only `R.FileName` cells for Negative-protein rows are rewritten. |
| `<stem>_sample_swap_ground_truth.tsv` | `Protein, Label` (full list, Positives + Negatives). |
| `<stem>_sample_swap_true_positives.tsv` | Positive proteins only. |
| `<stem>_sample_swap_group_annotation.csv` | `pair_id, cond_a, run_a, cond_b, run_b` — one row per pair. |

## Verification

- Both Python scripts parse cleanly (`ast.parse`).
- Smoke test on synthetic 8-row TSV (1 Positive, 1 Negative protein, 4 runs, 2 conditions):
  - `[pairs] 1 run pair (Condition1 <-> Condition2)` — `Seq2 <-> Seq4` (the every-second-run picks per the R helper).
  - POS_A rows: all 4 untouched.
  - NEG_B rows: `Seq2` rewritten to `Seq4` and vice versa; `Seq1`/`Seq3` (unpaired) untouched.

Confirms the R helper's behaviour 1:1: `r1 = Seq2, r2 = Seq4; partner[r1] = r2; partner[r2] = r1`; only Negative-protein rows get rewritten.

## End-to-end smoke test against the real TSV (deferred)

Per the §2.3 spec we should rerun MSstats on the original TSV (with the in-script label swap) and MSstats no-swap on the new TSV and confirm protein-level results agree. This requires the full ~7.4 GB Spectronaut TSV and a full MSstats run; deferred to §5 when `CSF_Spectronaut_sample_swap/` is wired up end-to-end with the shared pipeline from §1.

## Known design choices to revisit

1. **`R.Condition` is not rewritten** — only `R.FileName`. `R.Condition` in the TSV is the raw dilution label (`neat`, `1:2`, ...), used downstream only to drop blanks. The design factor used by every model is the annotation's `Condition` column, read after a fresh merge by `R.FileName`. Rewriting `R.FileName` is enough for the downstream pipeline to see the swap; rewriting `R.Condition` would only confuse anyone inspecting the TSV directly.
2. **Deterministic pairing.** No `--seed`. The R helper is deterministic too (`seq(2, n, by = 2)` and pairs in Order). Re-pairing or re-sampling the protein split across seeds would give uncertainty bands — flagged as future work in caveat §9.1 of CSF_swap_design.md.
3. **Two-condition only.** The script's pairing is between exactly two condition labels. For `<good_data>` (neat + 1/2 balanced) this is exactly the intended scope; for `<all_data>` the unpaired runs (1:4, 1:8, ...) are left untouched, which means their Negative-protein measurements stay attached to their original runs. Downstream subsetting is then responsible for either including or excluding those runs.

## Next item

§3 — CSF_Spectronaut output restructure: keep the author-faithful R scripts; reshape only the output paths into the §0.1 layout.
