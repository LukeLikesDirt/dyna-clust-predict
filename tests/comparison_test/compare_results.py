#!/usr/bin/env python3
"""Compare .predicted JSON files from dyna-clust-predict (R) and dnabarcoder (Python).

Usage:
    python3 compare_results.py <r_results_dir> <py_results_dir>

Compares F-measure profiles at every threshold and reports:
  - Whether F-measure curves match (they should, if the clustering logic is equivalent)
  - Whether selected cutoffs differ (expected due to tiebreaking: R picks middle, Python picks lowest)
"""

import json
import os
import sys
from pathlib import Path


def load_predicted(filepath):
    with open(filepath) as f:
        return json.load(f)


def normalise_fm_keys(fm):
    """Normalise threshold keys to 4-decimal format for consistent comparison."""
    return {f"{float(k):.4f}": v for k, v in fm.items()}


def compare_fmeasures(r_fm, py_fm):
    """Compare two fmeasure dictionaries {threshold_str: fmeasure_value}."""
    r_fm = normalise_fm_keys(r_fm)
    py_fm = normalise_fm_keys(py_fm)
    all_thresholds = sorted(set(list(r_fm.keys()) + list(py_fm.keys())), key=float)
    mismatches = []
    for t in all_thresholds:
        r_val = r_fm.get(t)
        py_val = py_fm.get(t)
        if r_val is None:
            mismatches.append((t, "MISSING in R", py_val))
        elif py_val is None:
            mismatches.append((t, r_val, "MISSING in Python"))
        elif abs(float(r_val) - float(py_val)) > 1e-4:
            mismatches.append((t, r_val, py_val))
    return mismatches, len(all_thresholds)


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <r_results_dir> <py_results_dir>")
        sys.exit(1)

    r_dir = Path(sys.argv[1])
    py_dir = Path(sys.argv[2])

    ranks = ["genus", "family", "order", "class", "phylum"]
    all_match = True

    print("=" * 78)
    print("COMPARISON: dyna-clust-predict (R) vs dnabarcoder (Python)")
    print("Both using the SAME vsearch-derived similarity matrix")
    print("=" * 78)
    print()

    for rank in ranks:
        r_file = r_dir / f"test_{rank}.predicted"
        py_file = py_dir / f"test_{rank}.predicted"

        print(f"── Rank: {rank} ──")

        if not r_file.exists():
            print(f"  R result not found: {r_file}")
            continue
        if not py_file.exists():
            print(f"  Python result not found: {py_file}")
            continue

        r_pred = load_predicted(r_file)
        py_pred = load_predicted(py_file)

        # Both should have the rank as a top-level key
        if rank not in r_pred:
            print(f"  R: rank '{rank}' not in prediction file")
            continue
        if rank not in py_pred:
            print(f"  Python: rank '{rank}' not in prediction file")
            continue

        r_datasets = r_pred[rank]
        py_datasets = py_pred[rank]

        # For global prediction, dataset should be "All"
        for dataset_name in sorted(set(list(r_datasets.keys()) + list(py_datasets.keys()))):
            if dataset_name not in r_datasets:
                print(f"  Dataset '{dataset_name}': MISSING in R results")
                continue
            if dataset_name not in py_datasets:
                print(f"  Dataset '{dataset_name}': MISSING in Python results")
                continue

            r_ds = r_datasets[dataset_name]
            py_ds = py_datasets[dataset_name]

            # Compare metadata
            r_cutoff = r_ds.get("cut-off", "N/A")
            py_cutoff = py_ds.get("cut-off", "N/A")
            r_conf = r_ds.get("confidence", "N/A")
            py_conf = py_ds.get("confidence", "N/A")
            r_seqno = r_ds.get("sequence number", "N/A")
            py_seqno = py_ds.get("sequence number", "N/A")
            r_groupno = r_ds.get("group number", "N/A")
            py_groupno = py_ds.get("group number", "N/A")

            print(f"  Dataset: {dataset_name}")
            print(f"    Sequences:  R={r_seqno}  Python={py_seqno}  {'MATCH' if r_seqno == py_seqno else 'DIFFER !'}")
            print(f"    Groups:     R={r_groupno}  Python={py_groupno}  {'MATCH' if r_groupno == py_groupno else 'DIFFER !'}")
            print(f"    Cutoff:     R={r_cutoff}  Python={py_cutoff}  {'MATCH' if r_cutoff == py_cutoff else 'DIFFER (may be tiebreaking)'}")
            print(f"    Confidence: R={r_conf}  Python={py_conf}  {'MATCH' if r_conf == py_conf else 'DIFFER !'}")

            if r_seqno != py_seqno or r_groupno != py_groupno:
                all_match = False
                print(f"    >>> WARNING: Different sequence/group counts indicate filtering differences!")

            # Compare F-measure curves
            r_fm = r_ds.get("fmeasures", {})
            py_fm = py_ds.get("fmeasures", {})

            if not r_fm and not py_fm:
                print(f"    No fmeasure data available for comparison")
                continue

            mismatches, total = compare_fmeasures(r_fm, py_fm)

            if len(mismatches) == 0:
                print(f"    F-measure curve: PERFECT MATCH ({total} thresholds)")
            else:
                all_match = False
                print(f"    F-measure curve: {len(mismatches)} MISMATCHES out of {total} thresholds:")
                # Show first 10 mismatches
                for t, r_val, py_val in mismatches[:10]:
                    print(f"      threshold={t}: R={r_val}  Python={py_val}")
                if len(mismatches) > 10:
                    print(f"      ... and {len(mismatches) - 10} more")
        print()

    print("=" * 78)
    if all_match:
        print("CONCLUSION: F-measure curves MATCH across all ranks.")
        print("Any cutoff differences are from tiebreaking strategy only")
        print("(R picks middle of tied range, Python picks lowest).")
        print("Your unexpected HPC results are likely from vsearch vs blastn,")
        print("NOT from other code changes.")
    else:
        print("CONCLUSION: Some F-measure curves DO NOT MATCH.")
        print("This means code changes beyond vsearch have altered the prediction logic.")
        print("Review the mismatched thresholds above to diagnose the issue.")
    print("=" * 78)


if __name__ == "__main__":
    main()
