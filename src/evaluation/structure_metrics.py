"""Base-pair metrics used to compare candidate and reference structures."""

from __future__ import annotations

from dataclasses import asdict, dataclass

from src.rna_utils import parse_dot_bracket


@dataclass(frozen=True)
class StructureMetrics:
    true_positives: int
    false_positives: int
    false_negatives: int
    precision: float
    recall: float
    f1: float
    exact_match: bool

    def as_dict(self) -> dict[str, int | float | bool]:
        return asdict(self)


def compare_structures(reference: str, candidate: str) -> StructureMetrics:
    if len(reference.strip()) != len(candidate.strip()):
        raise ValueError("Reference and candidate structure lengths differ")
    expected = parse_dot_bracket(reference)
    predicted = parse_dot_bracket(candidate)
    tp = len(expected & predicted)
    fp = len(predicted - expected)
    fn = len(expected - predicted)
    precision = tp / (tp + fp) if tp + fp else float(not expected)
    recall = tp / (tp + fn) if tp + fn else float(not predicted)
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return StructureMetrics(tp, fp, fn, precision, recall, f1, expected == predicted)
