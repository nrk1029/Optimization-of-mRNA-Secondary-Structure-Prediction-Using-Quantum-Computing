"""Dependency-free RNA and dot-bracket utilities."""

from __future__ import annotations

from collections.abc import Iterable

VALID_NUCLEOTIDES = frozenset("AUCG")
VALID_PAIRS = frozenset({"AU", "UA", "GC", "CG", "GU", "UG"})


def normalize_sequence(sequence: str) -> str:
    """Normalize an RNA sequence and reject unsupported characters."""
    normalized = "".join(sequence.upper().split()).replace("T", "U")
    if not normalized:
        raise ValueError("RNA sequence cannot be empty")
    invalid = sorted(set(normalized) - VALID_NUCLEOTIDES)
    if invalid:
        raise ValueError(f"Invalid RNA nucleotide(s): {', '.join(invalid)}")
    return normalized


def parse_dot_bracket(structure: str) -> set[tuple[int, int]]:
    """Return zero-based base pairs from pseudoknot-free dot-bracket notation."""
    stack: list[int] = []
    pairs: set[tuple[int, int]] = set()
    for index, symbol in enumerate(structure.strip()):
        if symbol == "(":
            stack.append(index)
        elif symbol == ")":
            if not stack:
                raise ValueError(f"Unmatched ')' at position {index}")
            pairs.add((stack.pop(), index))
        elif symbol != ".":
            raise ValueError(f"Unsupported dot-bracket symbol {symbol!r}")
    if stack:
        raise ValueError(f"Unmatched '(' at position(s): {stack}")
    return pairs


def validate_structure(sequence: str, structure: str) -> set[tuple[int, int]]:
    """Validate length and canonical/wobble pairing, then return base pairs."""
    sequence = normalize_sequence(sequence)
    structure = structure.strip()
    if len(sequence) != len(structure):
        raise ValueError("Sequence and structure lengths differ")
    pairs = parse_dot_bracket(structure)
    invalid = [(i, j) for i, j in pairs if sequence[i] + sequence[j] not in VALID_PAIRS]
    if invalid:
        raise ValueError(f"Invalid base pair(s): {invalid}")
    return pairs


def pairs_to_dot_bracket(length: int, pairs: Iterable[tuple[int, int]]) -> str:
    """Convert noncrossing, non-overlapping zero-based pairs to dot-bracket."""
    symbols = ["."] * length
    normalized = sorted(pairs)
    used: set[int] = set()
    for i, j in normalized:
        if not 0 <= i < j < length:
            raise ValueError(f"Pair {(i, j)} is outside sequence bounds")
        if i in used or j in used:
            raise ValueError("A nucleotide cannot occur in multiple pairs")
        if any(i < k < j < l or k < i < l < j for k, l in normalized):
            raise ValueError("Crossing pairs require pseudoknot notation")
        used.update((i, j))
        symbols[i], symbols[j] = "(", ")"
    return "".join(symbols)
