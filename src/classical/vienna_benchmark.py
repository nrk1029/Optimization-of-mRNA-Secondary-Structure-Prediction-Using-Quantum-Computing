"""Generate reproducible ViennaRNA MFE reference benchmarks."""

from __future__ import annotations

import argparse
import csv
import time
from pathlib import Path

from src.rna_utils import normalize_sequence, validate_structure


def fold_sequence(sequence: str) -> tuple[str, float, float]:
    sequence = normalize_sequence(sequence)
    try:
        import RNA
    except ImportError as exc:
        raise RuntimeError(
            "ViennaRNA is not installed. Run: python -m pip install -r requirements.txt"
        ) from exc
    started = time.perf_counter()
    structure, mfe = RNA.fold(sequence)
    runtime_ms = (time.perf_counter() - started) * 1000
    validate_structure(sequence, structure)
    return structure, float(mfe), runtime_ms


def generate_benchmarks(input_path: Path, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with input_path.open(newline="", encoding="utf-8") as source:
        rows = list(csv.DictReader(source))
    required = {"id", "sequence"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("Input CSV must contain id and sequence columns")
    fields = ["id", "sequence", "structure", "mfe_kcal_mol", "length", "runtime_ms"]
    with output_path.open("w", newline="", encoding="utf-8") as target:
        writer = csv.DictWriter(target, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            sequence = normalize_sequence(row["sequence"])
            structure, mfe, runtime_ms = fold_sequence(sequence)
            writer.writerow(
                {
                    "id": row["id"],
                    "sequence": sequence,
                    "structure": structure,
                    "mfe_kcal_mol": f"{mfe:.2f}",
                    "length": len(sequence),
                    "runtime_ms": f"{runtime_ms:.3f}",
                }
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    generate_benchmarks(args.input, args.output)


if __name__ == "__main__":
    main()
