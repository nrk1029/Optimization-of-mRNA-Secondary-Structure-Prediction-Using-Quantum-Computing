"""Compare QAOA depths and seeds on the validated 9-qubit RNA QUBO."""

from __future__ import annotations

import csv
import time
from pathlib import Path

import numpy as np
import RNA
from qiskit.primitives import StatevectorSampler
from qiskit_optimization import QuadraticProgram
from qiskit_optimization.algorithms import MinimumEigenOptimizer
from qiskit_optimization.minimum_eigensolvers import NumPyMinimumEigensolver, QAOA
from qiskit_optimization.optimizers import COBYLA
from qiskit_optimization.utils import algorithm_globals

SEQUENCE = "GGGAAACCC"
PAIR_ENERGIES = {"GC": -3.0, "CG": -3.0, "AU": -2.0, "UA": -2.0, "GU": -1.0, "UG": -1.0}
MIN_LOOP_SIZE = 3
PENALTY = 8.0
DEPTHS = (1, 2, 3)
SEEDS = (7, 21, 42)
MAX_ITERATIONS = 250
ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "results" / "qaoa_depth_seed_results.csv"


def candidates() -> list[dict[str, int | float | str]]:
    result = []
    for i in range(len(SEQUENCE)):
        for j in range(i + 1, len(SEQUENCE)):
            pair_type = SEQUENCE[i] + SEQUENCE[j]
            if j - i - 1 >= MIN_LOOP_SIZE and pair_type in PAIR_ENERGIES:
                result.append({"name": f"x_{i + 1}_{j + 1}", "i": i + 1, "j": j + 1, "weight": PAIR_ENERGIES[pair_type]})
    return result


def conflict(a: dict, b: dict) -> bool:
    i, j, k, l = int(a["i"]), int(a["j"]), int(b["i"]), int(b["j"])
    shared = bool({i, j} & {k, l})
    crossing = i < k < j < l or k < i < l < j
    return shared or crossing


def build_problem(items: list[dict]) -> QuadraticProgram:
    problem = QuadraticProgram("rna_qaoa_sweep")
    linear = {}
    quadratic = {}
    for item in items:
        problem.binary_var(str(item["name"]))
        linear[str(item["name"])] = float(item["weight"])
    for index, first in enumerate(items):
        for second in items[index + 1 :]:
            if conflict(first, second):
                quadratic[(str(first["name"]), str(second["name"]))] = PENALTY
    problem.minimize(linear=linear, quadratic=quadratic)
    return problem


def selected_names(result, items: list[dict]) -> list[str]:
    return [str(item["name"]) for bit, item in zip(result.x, items) if bit > 0.5]


def structure_from(selected: list[str], items: list[dict]) -> str:
    lookup = {str(item["name"]): item for item in items}
    structure = ["."] * len(SEQUENCE)
    for name in selected:
        item = lookup[name]
        structure[int(item["i"]) - 1] = "("
        structure[int(item["j"]) - 1] = ")"
    return "".join(structure)


def main() -> None:
    items = candidates()
    problem = build_problem(items)
    exact = MinimumEigenOptimizer(NumPyMinimumEigensolver()).solve(problem)
    exact_selected = selected_names(exact, items)
    exact_objective = float(exact.fval)
    mfe_structure, mfe_energy = RNA.fold(SEQUENCE)
    fold_compound = RNA.fold_compound(SEQUENCE)
    rows = []

    print(f"Exact objective: {exact_objective:.1f}; exact variables: {exact_selected}")
    for depth in DEPTHS:
        for seed in SEEDS:
            algorithm_globals.random_seed = seed
            rng = np.random.default_rng(seed)
            initial_point = rng.uniform(0.0, np.pi, size=2 * depth)
            solver = QAOA(
                sampler=StatevectorSampler(seed=seed),
                optimizer=COBYLA(maxiter=MAX_ITERATIONS),
                reps=depth,
                initial_point=initial_point,
            )
            started = time.perf_counter()
            result = MinimumEigenOptimizer(solver).solve(problem)
            runtime_seconds = time.perf_counter() - started
            selected = selected_names(result, items)
            structure = structure_from(selected, items)
            actual_energy = float(fold_compound.eval_structure(structure))
            matched = set(selected) == set(exact_selected)
            row = {
                "depth_p": depth,
                "seed": seed,
                "qubits": len(items),
                "objective": round(float(result.fval), 3),
                "objective_gap": round(float(result.fval) - exact_objective, 3),
                "structure": structure,
                "vienna_energy_kcal_mol": round(actual_energy, 3),
                "vienna_mfe_kcal_mol": round(float(mfe_energy), 3),
                "vienna_energy_gap_kcal_mol": round(actual_energy - float(mfe_energy), 3),
                "matched_exact_qubo": matched,
                "runtime_seconds": round(runtime_seconds, 3),
            }
            rows.append(row)
            print(f"p={depth} seed={seed}: objective={result.fval:.1f}, match={matched}, runtime={runtime_seconds:.2f}s")

    OUTPUT.parent.mkdir(exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    print(f"Reference structure: {mfe_structure}; MFE: {mfe_energy:.2f} kcal/mol")
    print(f"Saved: {OUTPUT}")


if __name__ == "__main__":
    main()
