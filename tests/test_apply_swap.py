"""Tests for the apply_swap path: G1 untouched, G2 picks up partner values."""

import polars as pl
import pytest

from swap_spectronaut_report import (
    apply_swap,
    assign_groups,
    build_pairs,
    build_rank_pairs,
    compute_protein_stats,
)


def _run_pipeline(tiny_report: pl.DataFrame, min_log2fc: float = 0.5,
                    swap_fraction: float = 1.0, seed: int = 1):
    prot_stats, prec_mean, pep_mean = compute_protein_stats(
        tiny_report, blank_cond="blank", min_precursors=2,
    )
    pairs = build_pairs(prot_stats, swap_fraction=swap_fraction,
                          min_log2fc=min_log2fc, seed=seed)
    if pairs.height == 0:
        pytest.skip("pipeline produced no pairs on this fixture")
    prec_match, dropped_prec = build_rank_pairs(pairs, prec_mean, "EG.PrecursorId", "prec_mean")
    pep_match,  dropped_pep  = build_rank_pairs(pairs, pep_mean,  "PEP.StrippedSequence", "pep_mean")
    groups = assign_groups(tiny_report, blank_cond="blank", seed=seed)
    swapped = apply_swap(tiny_report, groups, pairs, prec_match, pep_match,
                          dropped_prec, dropped_pep)
    return tiny_report, swapped, groups, pairs, prec_match


def test_g1_rows_untouched_for_paired_proteins(tiny_report):
    pre, post, groups, pairs, prec_match = _run_pipeline(tiny_report)

    g1_runs = groups.filter(pl.col("group") == "G1")["R.FileName"].to_list()
    paired_pgs = set(pairs["PG.ProteinGroups_hi"].to_list() +
                       pairs["PG.ProteinGroups_lo"].to_list())

    pre_g1 = pre.filter(
        pl.col("R.FileName").is_in(g1_runs) &
        pl.col("PG.ProteinGroups").is_in(list(paired_pgs))
    ).select(["R.FileName", "EG.PrecursorId", "FG.Quantity"])
    post_g1 = post.filter(
        pl.col("R.FileName").is_in(g1_runs) &
        pl.col("PG.ProteinGroups").is_in(list(paired_pgs))
    ).select(["R.FileName", "EG.PrecursorId", "FG.Quantity"])
    # Compare as sorted distinct rows — fragment-level duplication of
    # FG.Quantity makes raw equality flaky.
    pre_keys  = pre_g1.unique().sort(["R.FileName", "EG.PrecursorId"])
    post_keys = post_g1.unique().sort(["R.FileName", "EG.PrecursorId"])
    assert pre_keys.equals(post_keys)


def test_g2_partner_equality_for_paired_proteins(tiny_report):
    pre, post, groups, pairs, prec_match = _run_pipeline(tiny_report)

    g2_runs = groups.filter(pl.col("group") == "G2")["R.FileName"].to_list()
    if not g2_runs:
        pytest.skip("no G2 runs in this seed")

    # For each rank-aligned (hi_prec, lo_prec) pair, post[hi].FG.Quantity
    # in G2 should equal pre[lo].FG.Quantity in the same run (and vice versa).
    pre_fg = (pre.filter(pl.col("R.FileName").is_in(g2_runs))
                  .select(["R.FileName", "EG.PrecursorId", "FG.Quantity"])
                  .unique())
    post_fg = (post.filter(pl.col("R.FileName").is_in(g2_runs))
                    .select(["R.FileName", "EG.PrecursorId", "FG.Quantity"])
                    .unique())

    mismatches = 0
    checked = 0
    for row in prec_match.iter_rows(named=True):
        hi_prec = row["EG.PrecursorId_hi"]
        lo_prec = row["EG.PrecursorId_lo"]
        for run in g2_runs:
            # Look up the partner-before value
            partner_pre = pre_fg.filter(
                (pl.col("R.FileName") == run) &
                (pl.col("EG.PrecursorId") == lo_prec)
            )["FG.Quantity"]
            self_post = post_fg.filter(
                (pl.col("R.FileName") == run) &
                (pl.col("EG.PrecursorId") == hi_prec)
            )["FG.Quantity"]
            if partner_pre.len() == 0 or self_post.len() == 0:
                continue
            p, s = partner_pre.item(), self_post.item()
            if p is None or s is None:
                continue
            checked += 1
            if abs(p - s) > 1e-9:
                mismatches += 1

    assert checked > 0, "no comparable rows; fixture too sparse"
    assert mismatches == 0, f"{mismatches}/{checked} partner-equality mismatches in G2"


def test_blank_rows_not_grouped(tiny_report):
    _, _, groups, _, _ = _run_pipeline(tiny_report)
    # The assign_groups helper only considers non-blank runs.
    blank_rows = groups.filter(pl.col("R.FileName") == "Rb")
    assert blank_rows.height == 0
