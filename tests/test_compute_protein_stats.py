"""Tests for compute_protein_stats."""

import polars as pl
import pytest

from swap_spectronaut_report import compute_protein_stats


def test_returns_expected_columns(tiny_report):
    stats, prec_mean, pep_mean = compute_protein_stats(
        tiny_report, blank_cond="blank", min_precursors=2,
    )
    assert set(stats.columns) >= {
        "PG.ProteinGroups", "prot_mean_log2", "n_precursors", "n_peptides",
    }
    assert "EG.PrecursorId" in prec_mean.columns
    assert "PEP.StrippedSequence" in pep_mean.columns


def test_min_precursors_filters_out_singletons(tiny_report):
    # Each protein in tiny_report has 2 precursors. Bumping min_precursors
    # to 3 should leave the universe empty.
    stats, _, _ = compute_protein_stats(
        tiny_report, blank_cond="blank", min_precursors=3,
    )
    assert stats.height == 0


def test_min_precursors_2_keeps_all(tiny_report):
    stats, _, _ = compute_protein_stats(
        tiny_report, blank_cond="blank", min_precursors=2,
    )
    assert stats.height == 4
    assert set(stats["PG.ProteinGroups"].to_list()) == {"PH1", "PH2", "PL1", "PL2"}


def test_blank_condition_excluded(tiny_report):
    # tiny_report has FG.Quantity = NULL in the blank run. Even if it had
    # values, those should not contribute to prot_mean_log2.
    stats, _, _ = compute_protein_stats(
        tiny_report, blank_cond="blank", min_precursors=2,
    )
    # The PH proteins (~10000 intensity, log2 ≈ 13.3) should land above
    # the PL proteins (~500, log2 ≈ 8.9). If the blank null had leaked
    # into the mean it would either be NaN or drag the values down.
    means = dict(zip(stats["PG.ProteinGroups"], stats["prot_mean_log2"]))
    assert means["PH1"] > 12.0
    assert means["PL1"] < 10.0


def test_n_precursors_counts_distinct(tiny_report):
    stats, _, _ = compute_protein_stats(
        tiny_report, blank_cond="blank", min_precursors=2,
    )
    # Each protein has exactly 2 distinct EG.PrecursorId values.
    for n in stats["n_precursors"].to_list():
        assert n == 2
