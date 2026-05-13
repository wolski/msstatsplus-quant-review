"""Create a small subset of the swapped Spectronaut report + matching
annotation. Used to benchmark how methods (especially variance-moderation
ones) behave under small sample sizes.

Picks the runs explicitly so the resulting design is balanced
(8 G1 + 8 G2) and spans neat / 1to2 / 1to4 / 1to8 / 1to16 dilutions.

Usage:
  .venv/bin/python make_small_dataset.py
"""
from __future__ import annotations
from pathlib import Path
import polars as pl

HERE = Path(__file__).resolve().parent
SWAP_REPORT = HERE / "20250130_163144_CSF dilutions Jan 2025 no normalization_Report.tsv"
ORIG_REPORT = HERE / ".." / "CSF_Spectronaut" / SWAP_REPORT.name
FULL_ANNOTATION = HERE / "CSF_annotation.csv"
OUT_DIR = HERE / "small"
OUT_DIR.mkdir(exist_ok=True)
OUT_SWAP   = OUT_DIR / SWAP_REPORT.name
OUT_ORIG   = OUT_DIR / ("orig_" + SWAP_REPORT.name)
OUT_ANNOTATION = OUT_DIR / "CSF_annotation.csv"

# 6 neat (3 G1 + 3 G2), 4 1to2 good (2/2), 6 bad across 1to4/1to8/1to16 (1/1 each).
KEEP_RUNS = [
    # Neat — G1 (Condition1)
    "20250123_Tulum_NeatCSF-DD_Seq5",
    "20250123_Tulum_NeatCSF-DD_Seq11",
    "20250123_Tulum_NeatCSF-DD_Seq15",
    # Neat — G2 (Condition2)
    "20250123_Tulum_NeatCSF-DD_Seq1",
    "20250123_Tulum_NeatCSF-DD_Seq2",
    "20250123_Tulum_NeatCSF-DD_Seq10",
    # 1to2 good — G1
    "20250123_Tulum_1to2CSF-DD_Seq3",
    "20250123_Tulum_1to2CSF-DD_Seq7",
    # 1to2 good — G2
    "20250123_Tulum_1to2CSF-DD_Seq4",
    "20250123_Tulum_1to2CSF-DD_Seq12",
    # 1to4 bad — G1 + G2
    "20250123_Tulum_1to4CSF-DD_Seq31",
    "20250123_Tulum_1to4CSF-DD_Seq34",
    # 1to8 bad — G1 + G2
    "20250123_Tulum_1to8CSF-DD_Seq33",
    "20250123_Tulum_1to8CSF-DD_Seq36",
    # 1to16 bad — G1 + G2
    "20250123_Tulum_1to16CSF-DD_Seq35",
    "20250123_Tulum_1to16CSF-DD_Seq38",
]
assert len(KEEP_RUNS) == 16

ann = pl.read_csv(FULL_ANNOTATION)
missing = set(KEEP_RUNS) - set(ann["R.FileName"].to_list())
if missing:
    raise SystemExit(f"runs missing from annotation: {sorted(missing)}")

ann_subset = ann.filter(pl.col("R.FileName").is_in(KEEP_RUNS))
print(f"[ann] {ann_subset.height} runs:")
print(ann_subset.group_by("Condition").len())
ann_subset.write_csv(OUT_ANNOTATION)
print(f"[ann] wrote {OUT_ANNOTATION}")

def subset(src: Path, dst: Path, label: str) -> None:
    print(f"[{label}] subsetting {src.name}")
    df = (
        pl.scan_csv(src, separator="\t",
                    infer_schema_length=20000, null_values=["NaN", ""])
          .filter(pl.col("R.FileName").is_in(KEEP_RUNS))
          .collect()
    )
    print(f"[{label}] {df.height:,} rows kept")
    df.write_csv(dst, separator="\t",
                 include_header=True, quote_style="never",
                 null_value="NaN", line_terminator="\r\n")
    print(f"[{label}] wrote {dst}")

subset(SWAP_REPORT, OUT_SWAP, "swap")
subset(ORIG_REPORT, OUT_ORIG, "orig")
