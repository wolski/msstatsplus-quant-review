"""Tests for the --good-rule reference-run filter."""

import polars as pl
import pytest

from swap_spectronaut_report import (
    compute_protein_stats,
    filter_runs_by_good_rule,
)


def _ann(rows):
    return pl.DataFrame(rows, schema={
        "R.FileName": pl.String,
        "Condition":  pl.String,
        "Label":      pl.String,
    })


def test_label_good_keeps_good_only():
    ann = _ann([
        {"R.FileName": "r1", "Condition": "c1", "Label": "Good"},
        {"R.FileName": "r2", "Condition": "c1", "Label": "Bad"},
        {"R.FileName": "r3", "Condition": "c2", "Label": "Good"},
    ])
    out = filter_runs_by_good_rule(ann, "label_good")
    assert out["R.FileName"].to_list() == ["r1", "r3"]


def test_label_good_requires_label_column():
    ann = pl.DataFrame({"R.FileName": ["r1"], "Condition": ["c1"]})
    with pytest.raises(ValueError, match="Label"):
        filter_runs_by_good_rule(ann, "label_good")


def test_neat_only_matches_neatcsf_in_filename():
    ann = _ann([
        {"R.FileName": "20250123_NeatCSF-DD_Seq1", "Condition": "c1", "Label": "Good"},
        {"R.FileName": "20250123_1to2CSF-DD_Seq2", "Condition": "c2", "Label": "Good"},
        {"R.FileName": "20250123_NeatCSF-DD_Seq3", "Condition": "c1", "Label": "Bad"},
    ])
    out = filter_runs_by_good_rule(ann, "neat_only")
    assert out["R.FileName"].to_list() == [
        "20250123_NeatCSF-DD_Seq1", "20250123_NeatCSF-DD_Seq3",
    ]


def test_unknown_rule_raises():
    ann = _ann([{"R.FileName": "r1", "Condition": "c1", "Label": "Good"}])
    with pytest.raises(ValueError, match="Unknown good_rule"):
        filter_runs_by_good_rule(ann, "banana")


def test_reference_runs_restricts_protein_stats(tiny_report):
    # tiny_report has 4 proteins × 2 precursors × 2 frags × 4 runs (3 cond + blank).
    # Restrict the reference set to a single run; protein means should differ
    # from the unfiltered computation.
    stats_all, _, _ = compute_protein_stats(
        tiny_report, blank_cond="blank", min_precursors=2,
    )
    stats_one, _, _ = compute_protein_stats(
        tiny_report, blank_cond="blank", min_precursors=2,
        reference_runs=["R1"],
    )
    # Both should keep all 4 proteins.
    assert stats_all.height == 4
    assert stats_one.height == 4
    # The per-protein means must differ in general when the reference set shrinks
    # from 3 runs to 1 run.
    by_pg_all = dict(zip(stats_all["PG.ProteinGroups"], stats_all["prot_mean_log2"]))
    by_pg_one = dict(zip(stats_one["PG.ProteinGroups"], stats_one["prot_mean_log2"]))
    diffs = [abs(by_pg_all[pg] - by_pg_one[pg]) for pg in by_pg_all]
    assert max(diffs) > 1e-9, "expected per-protein abundance to shift when reference set shrinks"
