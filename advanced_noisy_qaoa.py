"""Evaluate reduced RNA QAOA under ideal and depolarizing noise models."""

from __future__ import annotations

import csv
import time
from pathlib import Path

import numpy as np
import RNA
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error
from qiskit_aer.primitives import SamplerV2
from qiskit_optimization import QuadraticProgram
from qiskit_optimization.algorithms import MinimumEigenOptimizer
from qiskit_optimization.minimum_eigensolvers import NumPyMinimumEigensolver, QAOA
from qiskit_optimization.optimizers import COBYLA
from qiskit_optimization.utils import algorithm_globals

SEQUENCE = "GGGAAACCC"
PAIR_ENERGIES = {"GC": -3.0, "CG": -3.0, "AU": -2.0, "UA": -2.0, "GU": -1.0, "UG": -1.0}
MIN_LOOP_SIZE = 3
PENALTY = 8.0
SEEDS = (7, 21, 42)
SHOTS = 2048
MAX_ITERATIONS = 150
ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "results" / "noisy_qaoa_results.csv"
CONDITIONS = (
    ("ideal", 0.0, 0.0),
    ("low_noise", 0.001, 0.01),
    ("high_noise", 0.005, 0.03),
)


def candidate_pairs() -> list[dict]:
    result = []
    for i in range(len(SEQUENCE)):
        for j in range(i + 1, len(SEQUENCE)):
            pair_type = SEQUENCE[i] + SEQUENCE[j]
            if j - i - 1 >= MIN_LOOP_SIZE and pair_type in PAIR_ENERGIES:
                result.append({"name": f"x_{i + 1}_{j + 1}", "i": i + 1, "j": j + 1, "weight": PAIR_ENERGIES[pair_type]})
    return result


def incompatible(a: dict, b: dict) -> bool:
    i, j, k, l = a["i"], a["j"], b["i"], b["j"]
    return bool({i, j} & {k, l}) or i < k < j < l or k < i < l < j


def build_problem(items: list[dict]) -> QuadraticProgram:
    problem = QuadraticProgram("noisy_rna_qaoa")
    linear, quadratic = {}, {}
    for item in items:
        problem.binary_var(item["name"])
        linear[item["name"]] = item["weight"]
    for index, first in enumerate(items):
        for second in items[index + 1 :]:
            if incompatible(first, second):
                quadratic[(first["name"], second["name"])] = PENALTY
    problem.minimize(linear=linear, quadratic=quadratic)
    return problem


def make_noise_model(one_qubit_error: float, two_qubit_error: float) -> NoiseModel | None:
    if one_qubit_error == 0 and two_qubit_error == 0:
        return None
    model = NoiseModel()
    model.add_all_qubit_quantum_error(depolarizing_error(one_qubit_error, 1), ["sx", "x"])
    model.add_all_qubit_quantum_error(depolarizing_error(two_qubit_error, 2), ["cx"])
    return model


def selected_names(result, items: list[dict]) -> list[str]:
    return [item["name"] for bit, item in zip(result.x, items) if bit > 0.5]


def structure_from(selected: list[str], items: list[dict]) -> str:
    lookup = {item["name"]: item for item in items}
    structure = ["."] * len(SEQUENCE)
    for name in selected:
        item = lookup[name]
        structure[item["i"] - 1] = "("
        structure[item["j"] - 1] = ")"
    return "".join(structure)


def main() -> None:
    items = candidate_pairs()
    problem = build_problem(items)
    exact = MinimumEigenOptimizer(NumPyMinimumEigensolver()).solve(problem)
    exact_selected = selected_names(exact, items)
    exact_objective = float(exact.fval)
    mfe_structure, mfe_energy = RNA.fold(SEQUENCE)
    fold_compound = RNA.fold_compound(SEQUENCE)
    rows = []

    for condition, error_1q, error_2q in CONDITIONS:
        noise_model = make_noise_model(error_1q, error_2q)
        backend = AerSimulator(noise_model=noise_model)
        pass_manager = generate_preset_pass_manager(optimization_level=1, backend=backend)
        for seed in SEEDS:
            algorithm_globals.random_seed = seed
            initial_point = np.random.default_rng(seed).uniform(0.0, np.pi, size=2)
            backend_options = {"method": "automatic"}
            if noise_model is not None:
                backend_options["noise_model"] = noise_model
            sampler = SamplerV2(
                default_shots=SHOTS,
                seed=seed,
                options={"backend_options": backend_options},
            )
            qaoa = QAOA(
                sampler=sampler,
                optimizer=COBYLA(maxiter=MAX_ITERATIONS),
                reps=1,
                initial_point=initial_point,
                pass_manager=pass_manager,
            )
            started = time.perf_counter()
            result = MinimumEigenOptimizer(qaoa).solve(problem)
            runtime = time.perf_counter() - started
            selected = selected_names(result, items)
            structure = structure_from(selected, items)
            actual_energy = float(fold_compound.eval_structure(structure))
            matched = set(selected) == set(exact_selected)
            rows.append({
                "condition": condition,
                "seed": seed,
                "shots": SHOTS,
                "one_qubit_error": error_1q,
                "two_qubit_error": error_2q,
                "objective": round(float(result.fval), 3),
                "objective_gap": round(float(result.fval) - exact_objective, 3),
                "structure": structure,
                "vienna_energy_kcal_mol": round(actual_energy, 3),
                "vienna_energy_gap_kcal_mol": round(actual_energy - float(mfe_energy), 3),
                "matched_exact_qubo": matched,
                "runtime_seconds": round(runtime, 3),
            })
            print(
                f"{condition:10s} seed={seed}: objective={result.fval:5.1f}, "
                f"match={matched}, Vienna gap={actual_energy - mfe_energy:4.1f}, runtime={runtime:.2f}s"
            )

    OUTPUT.parent.mkdir(exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"Reference: {mfe_structure}, {mfe_energy:.2f} kcal/mol")
    print(f"Saved: {OUTPUT}")


if __name__ == "__main__":
    main()
