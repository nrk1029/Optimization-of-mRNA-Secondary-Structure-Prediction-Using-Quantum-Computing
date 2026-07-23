# Quantum-Assisted RNA Secondary Structure Prediction

WISER x Moderna Quantum Challenge 2026 submission project.

Prepared by **NELLORE RAVI KUMAR**.

## Project summary

This project formulates pseudoknot-free RNA secondary-structure prediction as
a QUBO, benchmarks candidate structures with ViennaRNA, validates a reduced
9-qubit instance with exact optimization and QAOA, and analyzes scaling to the
44-nucleotide challenge sequence.

Key findings:

- ViennaRNA MFE for the challenge sequence: `-7.90 kcal/mol`
- Supplied candidate energy gap: `9.10 kcal/mol`
- Reduced exact QUBO reproduces `(((...)))`
- QAOA depth/seed sweep: 7 exact matches in 9 runs
- Direct 44-nt encoding: 313 qubits and 23,187 quadratic terms
- Stem encoding reduces variables from 313 to 200
- Tested noisy 9-qubit runs: 9 exact matches in 9 runs
- Full simplified QUBO is valid but thermodynamically weak: F1 `0.50`, energy
  gap `21.7 kcal/mol`

## Setup

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## Reproduce core outputs

```powershell
python -m unittest discover -s tests -v
python -m src.classical.vienna_benchmark --input data/sequences.csv --output results/vienna_benchmarks.csv
python step08_generate_plots.py
python step10_qaoa_sweep.py
python advanced_encoding_comparison.py
python advanced_noisy_qaoa.py
```

QAOA and noisy-simulator experiments may take several minutes.

## Main artifacts

- `REPORT.md`: full technical report
- `docs/Moderna_Stepwise_Documentation.docx`: step-by-step project documentation with code
- `docs/Moderna_Quantum_RNA_Presentation.pptx`: viva-ready project presentation
- `results/experiment_summary.csv`: consolidated benchmark summary
- `results/full_sequence_classical_qubo.csv`: 44-nt classical QUBO result
- `results/qaoa_depth_seed_results.csv`: QAOA robustness sweep
- `results/encoding_comparison.csv`: direct-vs-stem encoding comparison
- `results/noisy_qaoa_results.csv`: ideal and noisy simulator comparison
- `results/resource_scaling.png`: qubit/interaction/depth scaling
- `results/runtime_scaling.png`: classical preprocessing runtime

## Model scope

The direct QUBO uses simplified pair rewards and structural conflict penalties.
It is not a replacement for ViennaRNA's nearest-neighbor thermodynamic model.
Pseudoknots are excluded in the initial model. Only public and synthetic RNA
sequences are used.
