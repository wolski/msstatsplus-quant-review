"""build_subsets.py — slice a Spectronaut report + annotation into subset
directories. Each subset gets its own Report.tsv + annotation.csv inside
<out_dir>/<subset>/.

Subset rules:
  all_data         — every non-blank run
  good_data        — selected by --good-rule:
                       label_good : annotation Label == "Good"
                       neat_only  : R.FileName matches "NeatCSF"
  small            — 5 runs per Condition, stratified sample (seed=123)
  small_good_data  — 4 runs per Condition, stratified sample from good_data
                     (seed=123)

The annotation CSV is the source of truth for which runs to keep. The
Spectronaut TSV is filtered to rows whose R.FileName appears in the
filtered annotation. Blanks are dropped first (annotation Condition or
TSV R.Condition == "blank", case-insensitive).

Usage:
  python src/build_subsets.py \
    --report   <folder>/Report.tsv \
    --annotation <folder>/annotation.csv \
    --out-dir  <folder> \
    --subsets  all_data good_data small_good_data \
    --good-rule label_good
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Literal

import cyclopts
import polars as pl

Subset = Literal["all_data", "good_data", "small", "small_good_data"]
GoodRule = Literal["label_good", "neat_only"]

app = cyclopts.App(name="build-subsets", help=__doc__)


def drop_blanks(annotation: pl.DataFrame, report: pl.DataFrame) -> tuple[pl.DataFrame, pl.DataFrame]:
    ann = annotation.filter(pl.col("Condition").str.to_lowercase() != "blank")
    rep = report.filter(pl.col("R.Condition").str.to_lowercase() != "blank")
    return ann, rep


def filter_good(annotation: pl.DataFrame, good_rule: str) -> pl.DataFrame:
    if good_rule == "label_good":
        if "Label" not in annotation.columns:
            raise ValueError("good-rule=label_good requires a 'Label' column in annotation.")
        return annotation.filter(pl.col("Label") == "Good")
    if good_rule == "neat_only":
        return annotation.filter(pl.col("R.FileName").str.contains("NeatCSF"))
    raise ValueError(f"Unknown good_rule: {good_rule}")


def stratified_per_condition(annotation: pl.DataFrame, n: int, seed: int) -> pl.DataFrame:
    out = []
    for cond, group in annotation.group_by("Condition", maintain_order=True):
        take = min(n, group.height)
        sampled = group.sample(n=take, seed=seed, with_replacement=False)
        out.append(sampled)
    return pl.concat(out) if out else annotation.clear()


def write_subset(annotation: pl.DataFrame, report: pl.DataFrame,
                 out_dir: Path, subset: str) -> None:
    keep = set(annotation["R.FileName"].to_list())
    subset_report = report.filter(pl.col("R.FileName").is_in(keep))
    sub_dir = out_dir / subset
    sub_dir.mkdir(parents=True, exist_ok=True)
    annotation.write_csv(sub_dir / "annotation.csv")
    subset_report.write_csv(sub_dir / "Report.tsv", separator="\t")
    print(f"[build_subsets] {subset}: {annotation.height} runs, "
          f"{subset_report.height} report rows -> {sub_dir}/",
          file=sys.stderr)


@app.default
def main(
    *,
    report: Path,
    annotation: Path,
    out_dir: Path,
    subsets: list[Subset],
    good_rule: GoodRule | None = None,
    seed: int = 123,
    small_n: int = 5,
    small_good_n: int = 4,
) -> int:
    """Slice a Spectronaut report + annotation into subset directories.

    Parameters
    ----------
    report
        Spectronaut Report.tsv to read.
    annotation
        Annotation CSV with R.FileName + Condition (+ optional Label).
    out_dir
        Output directory; subset folders created under it.
    subsets
        One or more subsets to materialise. Pass multiple times,
        e.g. `--subsets all_data --subsets good_data`.
    good_rule
        Required when `--subsets` contains `good_data` or
        `small_good_data`. `label_good` keeps Label=="Good";
        `neat_only` keeps R.FileName matching "NeatCSF".
    seed
        RNG seed for stratified-sample subsets.
    small_n
        Runs per Condition in the `small` subset.
    small_good_n
        Runs per Condition in the `small_good_data` subset.
    """
    needs_good_rule = any(s in subsets for s in ("good_data", "small_good_data"))
    if needs_good_rule and good_rule is None:
        print("error: --good-rule is required when --subsets contains "
              "good_data or small_good_data", file=sys.stderr)
        return 2

    print(f"[build_subsets] reading {report} ...", file=sys.stderr)
    report_df = pl.read_csv(report, separator="\t",
                             infer_schema_length=10000, ignore_errors=True)
    print(f"[build_subsets] reading {annotation} ...", file=sys.stderr)
    annotation_df = pl.read_csv(annotation)

    annotation_df, report_df = drop_blanks(annotation_df, report_df)

    if "all_data" in subsets:
        write_subset(annotation_df, report_df, out_dir, "all_data")

    good = None
    if "good_data" in subsets or "small_good_data" in subsets:
        good = filter_good(annotation_df, good_rule)

    if "good_data" in subsets:
        write_subset(good, report_df, out_dir, "good_data")

    if "small" in subsets:
        small = stratified_per_condition(annotation_df, small_n, seed)
        write_subset(small, report_df, out_dir, "small")

    if "small_good_data" in subsets:
        sgd = stratified_per_condition(good, small_good_n, seed)
        write_subset(sgd, report_df, out_dir, "small_good_data")

    return 0


if __name__ == "__main__":
    app()
