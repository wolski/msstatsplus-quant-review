"""Tests for the stratified-by-n_precursors build_pairs."""

import polars as pl
import pytest

from swap_spectronaut_report import build_pairs


def test_no_protein_used_twice(simple_prot_stats):
    pairs = build_pairs(simple_prot_stats, swap_fraction=0.5,
                        min_log2fc=0.5, seed=42)
    used = pairs["PG.ProteinGroups_hi"].to_list() + pairs["PG.ProteinGroups_lo"].to_list()
    assert len(used) == len(set(used)), \
        f"some protein appears twice across all pairs: {used}"


def test_all_pairs_meet_min_log2fc(simple_prot_stats):
    pairs = build_pairs(simple_prot_stats, swap_fraction=0.5,
                        min_log2fc=1.5, seed=42)
    assert (pairs["log2fc"] >= 1.5).all()


def test_higher_min_log2fc_reduces_candidate_pool(simple_prot_stats):
    pairs_low  = build_pairs(simple_prot_stats, swap_fraction=1.0,
                              min_log2fc=0.1, seed=42)
    pairs_high = build_pairs(simple_prot_stats, swap_fraction=1.0,
                              min_log2fc=3.0, seed=42)
    assert pairs_low.height >= pairs_high.height


def test_smallest_fc_above_min_preferred(simple_prot_stats):
    # In the n=3 bin (P1=6.0, P2=7.0, P3=8.5, P4=10.0) the smallest log2fc
    # >= 0.5 is P2↔P1 = 1.0. The matcher should pick it before any larger fc.
    pairs = build_pairs(simple_prot_stats, swap_fraction=0.16,
                        min_log2fc=0.5, seed=42)
    n3_pairs = pairs.filter(pl.col("n_precursors") == 3)
    if n3_pairs.height > 0:
        smallest_picked = n3_pairs["log2fc"].min()
        assert smallest_picked == pytest.approx(1.0, abs=1e-6)


def test_stratification_keeps_high_n_bin(simple_prot_stats):
    # Without stratification, the small n=7 bin (only one possible pair Q1↔Q2)
    # could be drowned by n=3 candidates. With stratification + spillover off,
    # the n=7 bin gets its quota first.
    pairs = build_pairs(simple_prot_stats, swap_fraction=0.5,
                        min_log2fc=0.5, seed=42)
    n_values = sorted(pairs["n_precursors"].unique().to_list())
    assert 7 in n_values, \
        "n_precursors=7 bin should have at least one pair, not be drowned"


def test_spillover_when_high_n_bin_has_no_candidates(spill_prot_stats):
    # spill_prot_stats: n=7 has 2 proteins whose only candidate pair has
    # log2fc 0.3 < min_log2fc=1.0; quota of 1 spills into n=3.
    # swap_fraction = 0.3 ⇒ n=3 self_quota = round(0.3*6) = 2.
    # n=3 has six proteins arranged so the matcher can produce 3 disjoint
    # pairs at log2fc = 1.0 (P2↔P1, P4↔P3, P6↔P5). Without spillover the
    # matcher would stop at 2 (the self-quota). With spillover the quota
    # becomes 3, so n=3 should produce 3 pairs.
    pairs = build_pairs(spill_prot_stats, swap_fraction=0.3,
                        min_log2fc=1.0, seed=42)
    assert (pairs["n_precursors"] != 7).all()
    n3_pairs = pairs.filter(pl.col("n_precursors") == 3)
    assert n3_pairs.height == 3, (
        f"expected spillover to push n=3 pair count to 3 (self_quota=2 + "
        f"spill=1), got {n3_pairs.height}"
    )


def test_empty_input_returns_empty():
    empty = pl.DataFrame(schema={
        "PG.ProteinGroups": pl.String, "prot_mean_log2": pl.Float64,
        "n_precursors": pl.Int32, "n_peptides": pl.Int32,
    })
    pairs = build_pairs(empty, swap_fraction=0.1, min_log2fc=0.5, seed=42)
    assert pairs.height == 0
