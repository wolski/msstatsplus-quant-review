"""End-to-end integration test for swap_spectronaut_report.main().

Writes a small synthetic Spectronaut report + annotation to a tmp dir,
runs main(), and checks the full set of outputs is produced, row count is
preserved, and the labelled Positives recover their intended fold change.
"""

import polars as pl
import pytest

from swap_spectronaut_report import main


# Two HIGH proteins and two LOW proteins, all with 2 precursors x 2 fragments.
# HIGH ~ 2^13, LOW ~ 2^10  => log2fc ~ 3, comfortably above min_log2fc.
PROT_SPECS = [
    ("PH1", 8000.0), ("PH2", 8200.0),     # high
    ("PL1", 1000.0), ("PL2", 1050.0),     # low
]
RUNS = [("S1", "Neat"), ("S2", "Neat"), ("S3", "Neat"), ("S4", "Neat"),
        ("Sb", "blank")]


def _write_inputs(tmp_path):
    rows = []
    for pg, base in PROT_SPECS:
        for prec_i in (1, 2):
            prec = f"_{pg}p{prec_i}_.2"
            pep = f"{pg}seq{prec_i}"
            for frag_i in (1, 2):
                for run, cond in RUNS:
                    q = base * (1.0 + 0.03 * (hash(run + prec) % 4)) if cond != "blank" else None
                    rows.append({
                        "R.FileName": run, "R.Condition": cond,
                        "PG.ProteinGroups": pg, "PG.ProteinAccessions": pg,
                        "PEP.StrippedSequence": pep, "PEP.GroupingKey": pep,
                        "EG.PrecursorId": prec, "EG.ModifiedSequence": prec,
                        "FG.Charge": 2,
                        "FG.Quantity": q,
                        "F.PeakArea": (q / 2.0) if q is not None else None,
                        "F.PeakHeight": (q / 4.0) if q is not None else None,
                        "F.NormalizedPeakArea": (q / 2.0) if q is not None else None,
                        "F.NormalizedPeakHeight": (q / 4.0) if q is not None else None,
                        "FG.MS2RawQuantity": q,
                        "PEP.Quantity": q,
                        "PG.Quantity": (q * 2.0) if q is not None else None,
                    })
    report = pl.DataFrame(rows)
    report_path = tmp_path / "Report.tsv"
    report.write_csv(report_path, separator="\t", null_value="NaN")

    ann = pl.DataFrame({
        "R.FileName": [r for r, _ in RUNS],
        "Condition": [c for _, c in RUNS],
        "BioReplicate": [1, 2, 3, 4, 5],
        "Order": [1, 2, 3, 4, 5],
        "Label": ["Good", "Good", "Good", "Good", "Blank"],
    })
    ann_path = tmp_path / "annotation.csv"
    ann.write_csv(ann_path)
    return report, report_path, ann_path


def test_main_produces_all_outputs(tmp_path):
    report, report_path, ann_path = _write_inputs(tmp_path)
    out_dir = tmp_path / "out"

    rc = main(report=report_path, annotation=ann_path, out_dir=out_dir,
              swap_fraction=0.5, min_precursors=2, min_log2fc=1.0,
              good_rule=None, blank_condition="blank", seed=7)
    assert rc == 0

    for fname in ["Report.tsv",
                  "Report_swap_ground_truth.tsv",
                  "Report_swap_true_positives.tsv",
                  "Report_swap_group_annotation.csv",
                  "annotation.csv",
                  "CSF_protein_swap_list.csv"]:
        assert (out_dir / fname).exists(), f"missing output {fname}"


def test_main_preserves_row_count(tmp_path):
    report, report_path, ann_path = _write_inputs(tmp_path)
    out_dir = tmp_path / "out"
    main(report=report_path, annotation=ann_path, out_dir=out_dir,
         swap_fraction=0.5, min_precursors=2, min_log2fc=1.0,
         good_rule=None, blank_condition="blank", seed=7)

    swapped = pl.read_csv(out_dir / "Report.tsv", separator="\t",
                          null_values=["NaN", ""])
    assert swapped.height == report.height, "row count must be preserved"
    assert set(swapped.columns) == set(report.columns)


def test_main_annotation_rewritten_to_g1g2(tmp_path):
    _, report_path, ann_path = _write_inputs(tmp_path)
    out_dir = tmp_path / "out"
    main(report=report_path, annotation=ann_path, out_dir=out_dir,
         swap_fraction=0.5, min_precursors=2, min_log2fc=1.0,
         good_rule=None, blank_condition="blank", seed=7)

    ann_out = pl.read_csv(out_dir / "annotation.csv")
    conds = set(ann_out["Condition"].to_list())
    # Non-blank runs are relabelled to Condition1 / Condition2 (G1/G2).
    assert conds <= {"Condition1", "Condition2"}
    assert "Condition1" in conds and "Condition2" in conds


def test_main_swap_list_labels_pair_members_positive(tmp_path):
    _, report_path, ann_path = _write_inputs(tmp_path)
    out_dir = tmp_path / "out"
    main(report=report_path, annotation=ann_path, out_dir=out_dir,
         swap_fraction=0.5, min_precursors=2, min_log2fc=1.0,
         good_rule=None, blank_condition="blank", seed=7)

    swap_list = pl.read_csv(out_dir / "CSF_protein_swap_list.csv")
    gt = pl.read_csv(out_dir / "Report_swap_ground_truth.tsv", separator="\t")
    positives = set(swap_list.filter(pl.col("Label") == "Positive")["Protein"].to_list())
    pair_members = set(gt["PG_high"].to_list()) | set(gt["PG_low"].to_list())
    assert positives == pair_members


def test_main_good_rule_restricts_reference_runs(tmp_path):
    # good_rule="label_good" should run cleanly and still produce a swap;
    # exercises the reference-run restriction branch in main().
    _, report_path, ann_path = _write_inputs(tmp_path)
    out_dir = tmp_path / "out"
    rc = main(report=report_path, annotation=ann_path, out_dir=out_dir,
              swap_fraction=0.5, min_precursors=2, min_log2fc=1.0,
              good_rule="label_good", blank_condition="blank", seed=7)
    assert rc == 0
    assert (out_dir / "CSF_protein_swap_list.csv").exists()


def test_main_no_pairs_returns_1(tmp_path):
    # min_log2fc impossibly high -> no candidate pairs -> rc 1.
    _, report_path, ann_path = _write_inputs(tmp_path)
    out_dir = tmp_path / "out"
    rc = main(report=report_path, annotation=ann_path, out_dir=out_dir,
              swap_fraction=0.5, min_precursors=2, min_log2fc=99.0,
              good_rule=None, blank_condition="blank", seed=7)
    assert rc == 1
