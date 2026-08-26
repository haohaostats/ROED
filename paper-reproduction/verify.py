"""Integrity and numerical checks for the paper-reproduction archive."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def verify_hashes() -> None:
    for row in read_rows(ROOT / "SHA256SUMS.csv"):
        path = ROOT / row["path"]
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != row["sha256"]:
            raise RuntimeError(f"SHA-256 mismatch: {row['path']}")


def numeric(value: str):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def compare_csv(reference: Path, generated: Path) -> None:
    expected, actual = read_rows(reference), read_rows(generated)
    if len(expected) != len(actual):
        raise RuntimeError(f"Row-count mismatch: {reference.relative_to(ROOT)}")
    for index, (left, right) in enumerate(zip(expected, actual), 2):
        if left.keys() != right.keys():
            raise RuntimeError(f"Column mismatch: {reference.relative_to(ROOT)}")
        for column in left:
            # Wall-clock values are hardware-specific and are not numerical claims.
            if "time" in column or "elapsed" in column:
                continue
            x, y = numeric(left[column]), numeric(right[column])
            if x is not None and y is not None:
                if math.isnan(x) and math.isnan(y):
                    continue
                ref_name = reference.as_posix().lower()
                tolerance = 5e-4 if "merit" in ref_name else (5e-6 if "case-application" in ref_name else 1e-9)
                if not math.isclose(x, y, rel_tol=tolerance, abs_tol=tolerance):
                    raise RuntimeError(f"Numeric mismatch in {reference.name}:{index}:{column}")
            elif left[column] != right[column]:
                raise RuntimeError(f"Value mismatch in {reference.name}:{index}:{column}")


def verify_reference_structure() -> None:
    expected_counts = {
        "study1/study1_all_designs.csv": 120,
        "study2/study2_power_sample_metrics.csv": 120,
        "study3/study3_ablation_by_configuration.csv": 24,
        "study3/study3_ablation_aggregated.csv": 6,
        "study4/roed_s_comparison_by_configuration.csv": 24,
        "study4/roed_s_comparison_summary.csv": 6,
        "case-application/fixed_n_selection.csv": 2,
    }
    root = ROOT / "reference-results"
    for relative, count in expected_counts.items():
        rows = read_rows(root / relative)
        if len(rows) != count:
            raise RuntimeError(f"Expected {count} rows in {relative}; found {len(rows)}")
    selection = read_rows(root / "case-application" / "fixed_n_selection.csv")
    if [row["selected"] for row in selection] != ["TRUE", "FALSE"]:
        raise RuntimeError("Case-study selection decisions do not match the manuscript")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-only", action="store_true")
    args = parser.parse_args()
    verify_hashes()
    verify_reference_structure()
    if not args.reference_only:
        for reference in (ROOT / "reference-results").rglob("*.csv"):
            generated = ROOT / "results" / reference.relative_to(ROOT / "reference-results")
            if not generated.exists():
                raise RuntimeError(f"Missing generated CSV: {generated.relative_to(ROOT)}")
            compare_csv(reference, generated)
    print("Verification passed.")


if __name__ == "__main__":
    main()
