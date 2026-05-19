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

import argparse
import sys
from pathlib import Path

import polars as pl

KNOWN_SUBSETS = ("all_data", "good_data", "small", "small_good_data")
KNOWN_GOOD_RULES = ("label_good", "neat_only")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--report", required=True, type=Path)
    p.add_argument("--annotation", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--subsets", required=True, nargs="+", choices=KNOWN_SUBSETS)
    p.add_argument("--good-rule", choices=KNOWN_GOOD_RULES,
                   help="Required if 'good_data' or 'small_good_data' is in --subsets.")
    p.add_argument("--seed", type=int, default=123)
    p.add_argument("--small-n", type=int, default=5,
                   help="Runs per Condition for 'small' subset.")
    p.add_argument("--small-good-n", type=int, default=4,
                   help="Runs per Condition for 'small_good_data' subset.")
    return p.parse_args()


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


def main() -> int:
    args = parse_args()

    needs_good_rule = any(s in args.subsets for s in ("good_data", "small_good_data"))
    if needs_good_rule and args.good_rule is None:
        print("error: --good-rule is required when --subsets contains good_data or small_good_data",
              file=sys.stderr)
        return 2

    print(f"[build_subsets] reading {args.report} ...", file=sys.stderr)
    report = pl.read_csv(args.report, separator="\t",
                         infer_schema_length=10000, ignore_errors=True)
    print(f"[build_subsets] reading {args.annotation} ...", file=sys.stderr)
    annotation = pl.read_csv(args.annotation)

    annotation, report = drop_blanks(annotation, report)

    if "all_data" in args.subsets:
        write_subset(annotation, report, args.out_dir, "all_data")

    good = None
    if "good_data" in args.subsets or "small_good_data" in args.subsets:
        good = filter_good(annotation, args.good_rule)

    if "good_data" in args.subsets:
        write_subset(good, report, args.out_dir, "good_data")

    if "small" in args.subsets:
        small = stratified_per_condition(annotation, args.small_n, args.seed)
        write_subset(small, report, args.out_dir, "small")

    if "small_good_data" in args.subsets:
        sgd = stratified_per_condition(good, args.small_good_n, args.seed)
        write_subset(sgd, report, args.out_dir, "small_good_data")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
