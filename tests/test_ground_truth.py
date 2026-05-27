"""Tests for build_ground_truth (no drops in the relabel design)."""

import polars as pl
import pytest

from swap_spectronaut_report import build_ground_truth


def _pairs():
    # One pair: HIGH protein PA (mean log2 8), LOW protein PB (mean log2 6).
    return pl.DataFrame({
        "pair_id": [0],
        "PG.ProteinGroups_hi": ["PA"],
        "PG.ProteinGroups_lo": ["PB"],
        "prot_mean_log2_hi": [8.0],
        "prot_mean_log2_lo": [6.0],
        "n_precursors": [3],
        "n_peptides_hi": [3],
        "n_peptides_lo": [2],
        "log2fc": [2.0],
        "pep_count_diff": [1],
        "fc_above_min": [1.5],
    })


def _prec_match():
    return pl.DataFrame({
        "pair_id": [0, 0, 0],
        "EG.PrecursorId_hi": ["a1", "a2", "a3"],
        "EG.PrecursorId_lo": ["b1", "b2", "b3"],
    })


def _groups():
    return pl.DataFrame({
        "R.FileName": ["R1", "R2", "R3", "R4"],
        "R.Condition": ["Neat"] * 4,
        "group": ["G1", "G1", "G2", "G2"],
    })


def test_gt_columns_and_values():
    gt, tp, g_runs = build_ground_truth(_pairs(), _prec_match(), _groups())
    assert set(gt.columns) == {
        "pair_id", "PG_high", "PG_low", "n_precursors",
        "n_peptides_hi", "n_peptides_lo", "n_precursors_used", "expected_log2fc",
    }
    row = gt.to_dicts()[0]
    assert row["PG_high"] == "PA"
    assert row["PG_low"] == "PB"
    assert row["n_precursors"] == 3
    assert row["n_precursors_used"] == 3      # all precursors swapped, none dropped
    assert row["expected_log2fc"] == pytest.approx(2.0)


def test_tp_has_signed_roles():
    _, tp, _ = build_ground_truth(_pairs(), _prec_match(), _groups())
    # Two rows: the LOW protein is "up" (+log2fc), the HIGH protein is "down".
    by_prot = {r["PG.ProteinGroups"]: r for r in tp.to_dicts()}
    assert by_prot["PB"]["role"] == "up"
    assert by_prot["PB"]["expected_log2fc"] == pytest.approx(2.0)
    assert by_prot["PA"]["role"] == "down"
    assert by_prot["PA"]["expected_log2fc"] == pytest.approx(-2.0)
    assert by_prot["PA"]["partner_PG"] == "PB"
    assert by_prot["PB"]["partner_PG"] == "PA"


def test_g_runs_pivot():
    _, _, g_runs = build_ground_truth(_pairs(), _prec_match(), _groups())
    # One condition row with G1 and G2 run lists.
    assert g_runs.height == 1
    assert {"G1", "G2"}.issubset(set(g_runs.columns))
