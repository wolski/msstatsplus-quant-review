"""Shared test fixtures and helpers for src/swap_spectronaut_report.py.

The module under test lives in ../src/ which is not a Python package, so we
prepend it to sys.path before any test imports it.
"""

from __future__ import annotations

import sys
from pathlib import Path

import polars as pl
import pytest

# Make src/ importable as a flat namespace.
_HERE = Path(__file__).resolve().parent
_SRC = _HERE.parent / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))


# ---------------------------------------------------------------------------
# Synthetic protein stats frames (for build_pairs tests)
# ---------------------------------------------------------------------------

def make_prot_stats(rows: list[dict]) -> pl.DataFrame:
    """Build a prot_stats-shaped frame from a list of dicts.

    Each row needs: PG.ProteinGroups, prot_mean_log2, n_precursors, n_peptides.
    """
    return pl.DataFrame(rows, schema={
        "PG.ProteinGroups": pl.String,
        "prot_mean_log2":   pl.Float64,
        "n_precursors":     pl.Int32,
        "n_peptides":       pl.Int32,
    })


@pytest.fixture
def simple_prot_stats() -> pl.DataFrame:
    """6 proteins, two n_precursors bins (3 and 7).

    Bin n=3 has 4 proteins (P1..P4) at log2 mean 6, 7.0, 8.5, 10. So pairs:
      P1↔P3 (log2fc 2.5), P1↔P4 (log2fc 4.0),
      P2↔P3 (log2fc 1.5), P2↔P4 (log2fc 3.0),
      etc.
    Bin n=7 has 2 proteins (Q1, Q2) at log2 mean 5, 7.5 — only one possible pair.
    """
    return make_prot_stats([
        {"PG.ProteinGroups": "P1", "prot_mean_log2":  6.0, "n_precursors": 3, "n_peptides": 2},
        {"PG.ProteinGroups": "P2", "prot_mean_log2":  7.0, "n_precursors": 3, "n_peptides": 2},
        {"PG.ProteinGroups": "P3", "prot_mean_log2":  8.5, "n_precursors": 3, "n_peptides": 2},
        {"PG.ProteinGroups": "P4", "prot_mean_log2": 10.0, "n_precursors": 3, "n_peptides": 2},
        {"PG.ProteinGroups": "Q1", "prot_mean_log2":  5.0, "n_precursors": 7, "n_peptides": 6},
        {"PG.ProteinGroups": "Q2", "prot_mean_log2":  7.5, "n_precursors": 7, "n_peptides": 6},
    ])


@pytest.fixture
def spill_prot_stats() -> pl.DataFrame:
    """Designed to trigger spill-over: bin n=7 has two proteins too close in
    abundance to satisfy min_log2fc=1.0, so its quota rolls down to n=3.
    """
    return make_prot_stats([
        # n=7 bin: only one possible pair, with log2fc 0.3 < min_log2fc=1.0 ⇒
        # no valid candidate, quota rolls to n=3.
        {"PG.ProteinGroups": "Q1", "prot_mean_log2": 6.0, "n_precursors": 7, "n_peptides": 6},
        {"PG.ProteinGroups": "Q2", "prot_mean_log2": 6.3, "n_precursors": 7, "n_peptides": 6},
        # n=3 bin: lots of candidates.
        {"PG.ProteinGroups": "P1", "prot_mean_log2":  4.0, "n_precursors": 3, "n_peptides": 2},
        {"PG.ProteinGroups": "P2", "prot_mean_log2":  5.0, "n_precursors": 3, "n_peptides": 2},
        {"PG.ProteinGroups": "P3", "prot_mean_log2":  6.0, "n_precursors": 3, "n_peptides": 2},
        {"PG.ProteinGroups": "P4", "prot_mean_log2":  7.0, "n_precursors": 3, "n_peptides": 2},
        {"PG.ProteinGroups": "P5", "prot_mean_log2":  8.0, "n_precursors": 3, "n_peptides": 2},
        {"PG.ProteinGroups": "P6", "prot_mean_log2":  9.0, "n_precursors": 3, "n_peptides": 2},
    ])


# ---------------------------------------------------------------------------
# Synthetic Spectronaut report (for compute_protein_stats / apply_swap)
# ---------------------------------------------------------------------------

def make_report_rows(spec: list[dict]) -> pl.DataFrame:
    """Materialise a list of dicts into a Spectronaut-shaped frame.

    Each row needs the columns referenced by the script:
    R.FileName, R.Condition, PG.ProteinGroups, PEP.StrippedSequence,
    EG.PrecursorId, FG.Quantity, F.PeakArea, F.PeakHeight,
    F.NormalizedPeakArea, F.NormalizedPeakHeight, FG.MS2RawQuantity,
    PEP.Quantity, PG.Quantity.
    """
    return pl.DataFrame(spec)


@pytest.fixture
def tiny_report() -> pl.DataFrame:
    """4 proteins × 2 precursors × 2 fragments × 4 runs (3 cond + 1 blank).

    All four proteins have n_precursors = 2 so they live in the same bin.
    """
    runs = [
        ("R1", "cond1"),
        ("R2", "cond1"),
        ("R3", "cond2"),
        ("Rb", "blank"),
    ]
    proteins = [
        # (PG, base_intensity)
        ("PH1", 10000.0),
        ("PH2",  9000.0),
        ("PL1",   500.0),
        ("PL2",   400.0),
    ]
    rows = []
    for pg, base in proteins:
        for prec_i in (1, 2):
            prec = f"_{pg}p{prec_i}_.2"
            pep  = f"{pg}_seq{prec_i}"
            for frag_i in (1, 2):
                for run, cond in runs:
                    # Inject a tiny per-run perturbation so SD is non-trivial.
                    perturb = 1.0 + 0.05 * (hash(run + prec) % 5)
                    fg_q = base * perturb if cond != "blank" else None
                    rows.append({
                        "R.FileName":           run,
                        "R.Condition":          cond,
                        "PG.ProteinGroups":     pg,
                        "PEP.StrippedSequence": pep,
                        "EG.PrecursorId":       prec,
                        "FG.Quantity":          fg_q,
                        "F.PeakArea":           (fg_q or 0.0) / 10.0 * frag_i,
                        "F.PeakHeight":         (fg_q or 0.0) / 20.0 * frag_i,
                        "F.NormalizedPeakArea":   (fg_q or 0.0) / 10.0 * frag_i,
                        "F.NormalizedPeakHeight": (fg_q or 0.0) / 20.0 * frag_i,
                        "FG.MS2RawQuantity":    fg_q,
                        "PEP.Quantity":         fg_q,
                        "PG.Quantity":          (fg_q or 0.0) * 2,
                    })
    return pl.DataFrame(rows)
