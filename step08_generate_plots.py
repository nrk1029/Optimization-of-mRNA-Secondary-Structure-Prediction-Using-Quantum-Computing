"""Generate final scaling plots and a consolidated experiment summary."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
SCALING_CSV = ROOT / "scaling_results.csv"
RESULTS_DIR = ROOT / "results"


def load_scaling_rows() -> list[dict[str, str]]:
    with SCALING_CSV.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def save_experiment_summary() -> Path:
    output = RESULTS_DIR / "experiment_summary.csv"
    rows = [
        {
            "experiment": "official_document_candidate",
            "sequence": "GGAGCAAAACUUGUCGAUUGAGAACAAAAUACAGAAUUUGCUUG",
            "method": "provided_candidate",
            "structure": ".................(((....))).................",
            "objective": "",
            "vienna_energy_kcal_mol": "1.20",
            "mfe_energy_kcal_mol": "-7.90",
            "energy_gap_kcal_mol": "9.10",
            "matched_reference": "False",
        },
        {
            "experiment": "reduced_9nt",
            "sequence": "GGGAAACCC",
            "method": "exact_qubo",
            "structure": "(((...)))",
            "objective": "-9.0",
            "vienna_energy_kcal_mol": "-1.20",
            "mfe_energy_kcal_mol": "-1.20",
            "energy_gap_kcal_mol": "0.00",
            "matched_reference": "True",
        },
        {
            "experiment": "reduced_9nt",
            "sequence": "GGGAAACCC",
            "method": "qaoa_p1_seed42",
            "structure": ".((...).)",
            "objective": "-6.0",
            "vienna_energy_kcal_mol": "5.90",
            "mfe_energy_kcal_mol": "-1.20",
            "energy_gap_kcal_mol": "7.10",
            "matched_reference": "False",
        },
    ]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    return output


def plot_resource_scaling(rows: list[dict[str, str]]) -> Path:
    lengths = [int(row["sequence_length"]) for row in rows]
    qubits = [int(row["estimated_qubits"]) for row in rows]
    quadratic = [int(row["quadratic_terms"]) for row in rows]
    depth = [int(row["serial_depth_upper_bound"]) for row in rows]

    figure, axis = plt.subplots(figsize=(9, 5.5))
    axis.plot(lengths, qubits, marker="o", label="Logical qubits")
    axis.plot(lengths, quadratic, marker="s", label="Quadratic terms")
    axis.plot(lengths, depth, marker="^", label="Serial depth upper bound")
    axis.set_yscale("log")
    axis.set_xlabel("RNA sequence length (nt)")
    axis.set_ylabel("Resource count (log scale)")
    axis.set_title("Direct pair-variable QUBO resource scaling")
    axis.grid(True, which="both", alpha=0.25)
    axis.legend()
    figure.tight_layout()
    output = RESULTS_DIR / "resource_scaling.png"
    figure.savefig(output, dpi=220)
    plt.close(figure)
    return output


def plot_runtime(rows: list[dict[str, str]]) -> Path:
    lengths = [int(row["sequence_length"]) for row in rows]
    formulation = [float(row["formulation_runtime_ms"]) for row in rows]
    vienna = [float(row["vienna_runtime_ms"]) for row in rows]

    figure, axis = plt.subplots(figsize=(9, 5.5))
    axis.plot(lengths, formulation, marker="o", label="QUBO formulation")
    axis.plot(lengths, vienna, marker="s", label="ViennaRNA folding")
    axis.set_xlabel("RNA sequence length (nt)")
    axis.set_ylabel("Runtime (ms)")
    axis.set_title("Classical preprocessing runtime")
    axis.grid(True, alpha=0.25)
    axis.legend()
    figure.tight_layout()
    output = RESULTS_DIR / "runtime_scaling.png"
    figure.savefig(output, dpi=220)
    plt.close(figure)
    return output


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    rows = load_scaling_rows()
    outputs = [
        save_experiment_summary(),
        plot_resource_scaling(rows),
        plot_runtime(rows),
    ]
    print("Step 8 completed")
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
