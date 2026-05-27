"""Tests for the row-moving (identity-relabel) swap in apply_swap.

These build small, fully-controlled reports + an explicit `groups` frame and
`prec_match`, then call apply_swap directly so the assertions don't depend on
the pair-selection RNG.

Key properties under test:
  * row count is preserved exactly (nothing dropped, nothing inserted),
  * G1 rows are byte-identical before/after,
  * in G2, protein A ends up with protein B's intensities and vice versa,
  * with full detection the per-protein mean_diff is exactly antisymmetric:
    mean_diff(HIGH) == -mean_diff(LOW),
  * an NA hole in one partner fills the other partner symmetrically rather
    than dropping data,
  * the multiset of FG.Quantity within (pair x G2 run) is invariant.
"""

import math

import polars as pl
import pytest

from swap_spectronaut_report import apply_swap


# ---------------------------------------------------------------------------
# Controlled report builder
# ---------------------------------------------------------------------------

# A = HIGH abundance, B = LOW abundance. 3 precursors each, one fragment row
# per precursor per run (so FG.Quantity == the precursor-level value and the
# per-protein mean is a clean mean over 3 precursors).
PREC_SPECS = [
    ("A", "a1", "pepA1", 1000.0),
    ("A", "a2", "pepA2", 800.0),
    ("A", "a3", "pepA3", 600.0),
    ("B", "b1", "pepB1", 100.0),
    ("B", "b2", "pepB2", 80.0),
    ("B", "b3", "pepB3", 60.0),
]
RUNS = ["R1", "R2", "R3", "R4"]          # R1,R2 -> G1 ; R3,R4 -> G2


def make_report(missing: set[tuple[str, str]] | None = None) -> pl.DataFrame:
    """Build a controlled Spectronaut-shaped report.

    `missing` is a set of (run, precursor) pairs to omit (simulating an
    undetected precursor in that run).
    """
    missing = missing or set()
    rows = []
    for pg, prec, pep, base in PREC_SPECS:
        for run in RUNS:
            if (run, prec) in missing:
                continue
            rows.append({
                "R.FileName": run,
                "R.Condition": "Neat",
                "PG.ProteinGroups": pg,
                "PEP.StrippedSequence": pep,
                "EG.PrecursorId": prec,
                "EG.ModifiedSequence": f"_{prec}_",
                "FG.Charge": 2,
                "FG.Quantity": base,
                "F.PeakArea": base / 2.0,
            })
    return pl.DataFrame(rows)


def make_groups() -> pl.DataFrame:
    return pl.DataFrame({
        "R.FileName": RUNS,
        "R.Condition": ["Neat"] * 4,
        "group": ["G1", "G1", "G2", "G2"],
    })


def make_prec_match() -> pl.DataFrame:
    # Rank-aligned: a1<->b1, a2<->b2, a3<->b3.
    return pl.DataFrame({
        "pair_id": [0, 0, 0],
        "EG.PrecursorId_hi": ["a1", "a2", "a3"],
        "EG.PrecursorId_lo": ["b1", "b2", "b3"],
    })


def _protein_run_mean_log2(df: pl.DataFrame) -> pl.DataFrame:
    return (
        df.group_by(["PG.ProteinGroups", "R.FileName"])
          .agg(log2I=pl.col("FG.Quantity").log(2).mean())
    )


def _mean_diff(df: pl.DataFrame, groups: pl.DataFrame) -> dict[str, float]:
    """Per-protein mean_diff = mean_G1 - mean_G2 of per-run log2 means."""
    pr = _protein_run_mean_log2(df).join(
        groups.select("R.FileName", "group"), on="R.FileName")
    agg = pr.group_by("PG.ProteinGroups", "group").agg(m=pl.col("log2I").mean())
    wide = agg.pivot(on="group", values="m", index="PG.ProteinGroups")
    return {
        row["PG.ProteinGroups"]: row["G1"] - row["G2"]
        for row in wide.iter_rows(named=True)
    }


# ---------------------------------------------------------------------------
# Tests — full detection
# ---------------------------------------------------------------------------

def test_row_count_preserved():
    df = make_report()
    out = apply_swap(df, make_groups(), make_prec_match())
    assert out.height == df.height, "relabel swap must not add or drop rows"


def test_columns_preserved():
    df = make_report()
    out = apply_swap(df, make_groups(), make_prec_match())
    assert set(out.columns) == set(df.columns)


def test_g1_untouched():
    df = make_report()
    out = apply_swap(df, make_groups(), make_prec_match())
    g1 = ["R1", "R2"]
    pre = df.filter(pl.col("R.FileName").is_in(g1)).sort(
        ["R.FileName", "PG.ProteinGroups", "EG.PrecursorId"])
    post = out.filter(pl.col("R.FileName").is_in(g1)).sort(
        ["R.FileName", "PG.ProteinGroups", "EG.PrecursorId"])
    assert pre.equals(post)


def test_g2_proteins_take_partner_intensities():
    df = make_report()
    out = apply_swap(df, make_groups(), make_prec_match())
    # In a G2 run, protein A should now carry B's intensities (100/80/60)
    # and protein B should carry A's (1000/800/600).
    g2 = out.filter(pl.col("R.FileName") == "R3")
    a_vals = sorted(g2.filter(pl.col("PG.ProteinGroups") == "A")["FG.Quantity"].to_list())
    b_vals = sorted(g2.filter(pl.col("PG.ProteinGroups") == "B")["FG.Quantity"].to_list())
    assert a_vals == [60.0, 80.0, 100.0]
    assert b_vals == [600.0, 800.0, 1000.0]


def test_mean_diff_is_antisymmetric_under_full_detection():
    df = make_report()
    out = apply_swap(df, make_groups(), make_prec_match())
    md = _mean_diff(out, make_groups())
    assert md["A"] == pytest.approx(-md["B"], abs=1e-9), \
        f"expected antisymmetric mean_diff, got A={md['A']}, B={md['B']}"
    # And the magnitude equals the true log2fc between A and B.
    expected = (
        (math.log2(1000) + math.log2(800) + math.log2(600)) / 3
        - (math.log2(100) + math.log2(80) + math.log2(60)) / 3
    )
    assert md["A"] == pytest.approx(expected, abs=1e-9)


def test_fg_multiset_preserved_within_pair_per_g2_run():
    df = make_report()
    out = apply_swap(df, make_groups(), make_prec_match())
    for run in ["R3", "R4"]:
        pre = sorted(df.filter(pl.col("R.FileName") == run)["FG.Quantity"].to_list())
        post = sorted(out.filter(pl.col("R.FileName") == run)["FG.Quantity"].to_list())
        assert pre == post, f"FG.Quantity multiset changed in {run}"


# ---------------------------------------------------------------------------
# Tests — NA hole handling
# ---------------------------------------------------------------------------

def test_na_hole_fills_partner_not_dropped():
    # b3 is undetected in G2 run R3. The relabel should MOVE a3's row to
    # protein B (B gains b3's slot with A's value) and protein A loses a3 in
    # R3 (no b3 to source from). No row is dropped relative to input.
    df = make_report(missing={("R3", "b3")})
    out = apply_swap(df, make_groups(), make_prec_match())

    assert out.height == df.height, "no rows should be dropped"

    r3 = out.filter(pl.col("R.FileName") == "R3")
    a_precs = set(r3.filter(pl.col("PG.ProteinGroups") == "A")["EG.PrecursorId"].to_list())
    b_precs = set(r3.filter(pl.col("PG.ProteinGroups") == "B")["EG.PrecursorId"].to_list())
    # A lost a3 in R3 (its source b3 was missing); B has all three including b3.
    assert a_precs == {"a1", "a2"}
    assert b_precs == {"b1", "b2", "b3"}
    # B's b3 carries A's a3 value (600) — the value moved, not dropped.
    b3_val = r3.filter((pl.col("PG.ProteinGroups") == "B") &
                       (pl.col("EG.PrecursorId") == "b3"))["FG.Quantity"].item()
    assert b3_val == 600.0


def test_na_hole_multiset_still_preserved():
    df = make_report(missing={("R3", "b3")})
    out = apply_swap(df, make_groups(), make_prec_match())
    for run in ["R3", "R4"]:
        pre = sorted(df.filter(pl.col("R.FileName") == run)["FG.Quantity"].to_list())
        post = sorted(out.filter(pl.col("R.FileName") == run)["FG.Quantity"].to_list())
        assert pre == post


# ---------------------------------------------------------------------------
# Tests — identity relabel specifics
# ---------------------------------------------------------------------------

def test_identity_columns_relabelled_together():
    # When a3's row becomes b3, ALL of b3's identity columns (peptide,
    # modified sequence) must move as a unit, not just the precursor id.
    df = make_report()
    out = apply_swap(df, make_groups(), make_prec_match())
    g2 = out.filter(pl.col("R.FileName") == "R3")
    # Find the row now labelled protein B, precursor b1 — it should carry
    # B-side identity throughout and A's a1 intensity (1000).
    row = g2.filter((pl.col("PG.ProteinGroups") == "B") &
                    (pl.col("EG.PrecursorId") == "b1"))
    assert row.height == 1
    r = row.to_dicts()[0]
    assert r["PEP.StrippedSequence"] == "pepB1"
    assert r["EG.ModifiedSequence"] == "_b1_"
    assert r["FG.Quantity"] == 1000.0     # A's a1 intensity moved here


def test_negatives_untouched():
    # Add a Negative protein C not in prec_match; it must be unchanged in G2.
    df = make_report()
    extra = pl.DataFrame([
        {"R.FileName": run, "R.Condition": "Neat", "PG.ProteinGroups": "C",
         "PEP.StrippedSequence": "pepC1", "EG.PrecursorId": "c1",
         "EG.ModifiedSequence": "_c1_", "FG.Charge": 2,
         "FG.Quantity": 250.0, "F.PeakArea": 125.0}
        for run in RUNS
    ])
    df = pl.concat([df, extra])
    out = apply_swap(df, make_groups(), make_prec_match())
    pre_c = df.filter(pl.col("PG.ProteinGroups") == "C").sort("R.FileName")
    post_c = out.filter(pl.col("PG.ProteinGroups") == "C").sort("R.FileName")
    assert pre_c.equals(post_c)
