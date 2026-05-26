"""TSV-level sample swap for Spectronaut reports.

Implements the cleaner alternative described in §6 of
quant/CSF_Spectronaut/CSF_swap_design.md: for every row whose protein is on
the Negative list, replace R.FileName with its paired partner's R.FileName.
The intensity columns are NOT touched - swapping the run identifier is
equivalent to swapping intensities between paired runs, because downstream
code joins the TSV with the annotation by R.FileName and reads the
experimental Condition from the annotation.

Pairing mirrors swap_condition_labels_msstats() in
benchmark_experiments_functions.R: sort runs by Order, pick every second
run within each condition (R: seq(2, n, by = 2); 0-based: indices 1, 3, 5,
...), pair index-for-index across the two conditions.

Outputs (next to <stem>):
  <stem>_sample_swap.tsv              schema-identical Spectronaut TSV
  <stem>_sample_swap_ground_truth.tsv per-protein Label (Positive/Negative)
  <stem>_sample_swap_true_positives.tsv  Positive proteins only
  <stem>_sample_swap_group_annotation.csv pairing details (cond, run -> partner)
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import cyclopts
import polars as pl

app = cyclopts.App(name="swap-spectronaut-report-samples", help=__doc__,
                    help_on_error=True)

COL_RUN = "R.FileName"
COL_COND_RAW = "R.Condition"
COL_PG = "PG.ProteinGroups"


# CLI moved below to @app.default on main().


def build_pairing(annotation: pl.DataFrame, cond_a: str, cond_b: str,
                  blank: str) -> tuple[dict[str, str], pl.DataFrame]:
    """Build run -> partner-run map for the two named conditions.

    Returns (partner_map, pairing_table). Runs whose Condition is neither
    cond_a nor cond_b (or which fall outside the every-second-run subset) are
    absent from partner_map; those rows in the TSV will be left untouched.
    """
    a = annotation.filter(pl.col("Condition") == cond_a).sort("Order")
    b = annotation.filter(pl.col("Condition") == cond_b).sort("Order")
    # seq(2, n, by = 2) in R (1-based) -> 0-based indices 1, 3, 5, ...
    a_runs = a[COL_RUN].to_list()[1::2]
    b_runs = b[COL_RUN].to_list()[1::2]
    n_pairs = min(len(a_runs), len(b_runs))
    a_runs = a_runs[:n_pairs]
    b_runs = b_runs[:n_pairs]

    partner: dict[str, str] = {}
    for ra, rb in zip(a_runs, b_runs):
        partner[ra] = rb
        partner[rb] = ra

    pairing_rows = []
    for i, (ra, rb) in enumerate(zip(a_runs, b_runs), start=1):
        pairing_rows.append({"pair_id": i, "cond_a": cond_a, "run_a": ra,
                             "cond_b": cond_b, "run_b": rb})
    pairing_tbl = pl.DataFrame(pairing_rows) if pairing_rows else pl.DataFrame(
        {"pair_id": [], "cond_a": [], "run_a": [], "cond_b": [], "run_b": []}
    )
    return partner, pairing_tbl


def apply_run_swap(df: pl.DataFrame, partner: dict[str, str],
                    negative_proteins: set[str]) -> pl.DataFrame:
    """For rows whose PG.ProteinGroups is in negative_proteins, rewrite
    R.FileName to its partner. Rows with no partner (e.g. unpaired runs,
    blanks) pass through untouched."""
    if not partner or not negative_proteins:
        return df

    map_df = pl.DataFrame(
        {"_self": list(partner.keys()), "_partner": list(partner.values())}
    )
    df = (
        df.with_columns(
            pl.col(COL_PG).is_in(list(negative_proteins)).alias("_is_neg")
        )
        .join(map_df, left_on=COL_RUN, right_on="_self", how="left")
        .with_columns(
            pl.when(pl.col("_is_neg") & pl.col("_partner").is_not_null())
              .then(pl.col("_partner"))
              .otherwise(pl.col(COL_RUN))
              .alias(COL_RUN)
        )
        .drop(["_is_neg", "_partner"])
    )
    return df


@app.default
def main(
    *,
    report: Path,
    annotation: Path,
    protein_swap_list: Path,
    out_dir: Path,
    cond_a: str = "Condition1",
    cond_b: str = "Condition2",
    blank_condition: str = "Blank",
) -> int:
    """TSV-level sample swap for Spectronaut reports.

    Parameters
    ----------
    report
        Spectronaut TSV (pre-swap).
    annotation
        Annotation CSV with R.FileName, Condition, Order columns.
    protein_swap_list
        CSV with Protein, Label columns; Label in {Positive, Negative}.
    out_dir
        Output directory; created if missing.
    cond_a, cond_b
        Annotation Condition values to pair across.
    blank_condition
        Annotation Condition value identifying blank runs (excluded
        from pairing).
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[load] report: {report}", file=sys.stderr)
    df = pl.read_csv(report, separator="\t", null_values=["NaN", ""],
                     infer_schema_length=20000, low_memory=False)
    print(f"[load] {df.height:,} rows x {df.width} cols", file=sys.stderr)

    annotation_df = pl.read_csv(annotation)
    swap_list = pl.read_csv(protein_swap_list)
    assert {"Protein", "Label"}.issubset(swap_list.columns), \
        "protein-swap-list must have Protein,Label columns"

    negatives = set(swap_list.filter(pl.col("Label") == "Negative")["Protein"].to_list())
    positives = set(swap_list.filter(pl.col("Label") == "Positive")["Protein"].to_list())
    print(f"[truth] {len(positives)} positives, {len(negatives)} negatives",
          file=sys.stderr)

    partner, pairing_tbl = build_pairing(annotation_df, cond_a, cond_b,
                                          blank_condition)
    print(f"[pairs] {len(partner)//2} run pairs ({cond_a} <-> {cond_b})",
          file=sys.stderr)

    original_cols = df.columns
    df_swapped = apply_run_swap(df, partner, negatives).select(original_cols)

    stem = report.stem
    # Canonical Makefile contract: write Report.tsv + annotation.csv directly.
    out_report = out_dir / "Report.tsv"
    out_annotation = out_dir / "annotation.csv"
    out_gt = out_dir / f"{stem}_sample_swap_ground_truth.tsv"
    out_tp = out_dir / f"{stem}_sample_swap_true_positives.tsv"
    out_pairs = out_dir / f"{stem}_sample_swap_group_annotation.csv"

    print(f"[write] {out_report}", file=sys.stderr)
    df_swapped.write_csv(out_report, separator="\t", null_value="NaN",
                          include_header=True, quote_style="never",
                          line_terminator="\r\n")

    print(f"[write] {out_annotation}", file=sys.stderr)
    shutil.copy(annotation, out_annotation)

    print(f"[write] {out_gt}", file=sys.stderr)
    swap_list.write_csv(out_gt, separator="\t")

    print(f"[write] {out_tp}", file=sys.stderr)
    (swap_list.filter(pl.col("Label") == "Positive")
              .write_csv(out_tp, separator="\t"))

    print(f"[write] {out_pairs}", file=sys.stderr)
    pairing_tbl.write_csv(out_pairs)

    return 0


if __name__ == "__main__":
    app()
