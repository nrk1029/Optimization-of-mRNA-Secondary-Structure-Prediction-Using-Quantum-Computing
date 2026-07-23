"""Compare direct base-pair and contiguous-stem QUBO encodings."""

from __future__ import annotations

import csv
import time
from pathlib import Path

import RNA

OFFICIAL_SEQUENCE = "GGAGCAAAACUUGUCGAUUGAGAACAAAAUACAGAAUUUGCUUG"
SEQUENCES = [
    ("toy_9", "GGGAAACCC"),
    *[(f"official_prefix_{n}", OFFICIAL_SEQUENCE[:n]) for n in (12, 16, 20, 28, 36, 44)],
]
VALID_PAIRS = {"AU", "UA", "GC", "CG", "GU", "UG"}
MIN_LOOP_SIZE = 3
MIN_STEM_LENGTH = 2
ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "results" / "encoding_comparison.csv"


def direct_pairs(sequence: str) -> list[tuple[int, int]]:
    return [
        (i, j)
        for i in range(len(sequence))
        for j in range(i + 1, len(sequence))
        if j - i - 1 >= MIN_LOOP_SIZE and sequence[i] + sequence[j] in VALID_PAIRS
    ]


def stems(sequence: str, pairs: list[tuple[int, int]]) -> list[tuple[tuple[int, int], ...]]:
    result = []
    for i, j in pairs:
        maximum = 1
        while i + maximum < j - maximum and sequence[i + maximum] + sequence[j - maximum] in VALID_PAIRS:
            maximum += 1
        for length in range(MIN_STEM_LENGTH, maximum + 1):
            result.append(tuple((i + offset, j - offset) for offset in range(length)))
    return result


def dot_bracket_pairs(structure: str) -> set[tuple[int, int]]:
    stack = []
    result = set()
    for index, symbol in enumerate(structure):
        if symbol == "(":
            stack.append(index)
        elif symbol == ")":
            result.add((stack.pop(), index))
    return result


def pair_conflict(first: tuple[int, int], second: tuple[int, int]) -> bool:
    i, j = first
    k, l = second
    return bool({i, j} & {k, l}) or i < k < j < l or k < i < l < j


def object_conflict(first: tuple[tuple[int, int], ...], second: tuple[tuple[int, int], ...]) -> bool:
    return any(pair_conflict(a, b) for a in first for b in second)


def conflict_count(objects: list[tuple[tuple[int, int], ...]]) -> int:
    count = 0
    for index, first in enumerate(objects):
        for second in objects[index + 1 :]:
            if object_conflict(first, second):
                count += 1
    return count


def main() -> None:
    rows = []
    for sequence_id, sequence in SEQUENCES:
        started = time.perf_counter()
        pairs = direct_pairs(sequence)
        pair_objects = [((i, j),) for i, j in pairs]
        pair_conflicts = conflict_count(pair_objects)
        pair_runtime = (time.perf_counter() - started) * 1000

        started = time.perf_counter()
        stem_objects = stems(sequence, pairs)
        stem_conflicts = conflict_count(stem_objects)
        stem_runtime = (time.perf_counter() - started) * 1000

        mfe_structure, _ = RNA.fold(sequence)
        mfe_pairs = dot_bracket_pairs(mfe_structure)
        representable_by_stems = set().union(*(set(stem) for stem in stem_objects)) if stem_objects else set()
        coverage = len(mfe_pairs & representable_by_stems) / len(mfe_pairs) if mfe_pairs else 1.0
        reduction = 1.0 - len(stem_objects) / len(pairs) if pairs else 0.0

        rows.append({
            "sequence_id": sequence_id,
            "sequence_length": len(sequence),
            "direct_pair_variables": len(pairs),
            "direct_pair_conflicts": pair_conflicts,
            "direct_formulation_ms": round(pair_runtime, 3),
            "stem_variables": len(stem_objects),
            "stem_conflicts": stem_conflicts,
            "stem_formulation_ms": round(stem_runtime, 3),
            "variable_reduction_percent": round(100 * reduction, 2),
            "mfe_pair_coverage_percent": round(100 * coverage, 2),
        })
        print(
            f"n={len(sequence):2d}: direct={len(pairs):3d}, stems={len(stem_objects):3d}, "
            f"reduction={100 * reduction:5.1f}%, MFE coverage={100 * coverage:5.1f}%"
        )

    OUTPUT.parent.mkdir(exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved: {OUTPUT}")


if __name__ == "__main__":
    main()
