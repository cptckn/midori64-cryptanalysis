# Reproducibility protocol

## What is independently checked

The included evidence has three different evidentiary roles. Do not conflate them.

1. **Feasible primal assignment:** `results/active_solution.sol` contains 1,280 binary values. `python audit_solution.py` regenerates the *same* four-round active-S-box model from the checked-in 19 H-representation inequalities and verifies all 2,049 constraints. This check does **not** by itself establish global optimality.
2. **Computational optimality:** `results/active_solver.log` records Gurobi 13.0.3 returning `OPTIMAL`, objective 16, best bound 16, and gap zero after 794.455 seconds. `results/active_optimality.json` summarizes this completed run. This lower-bound result relies on the solver and correct MILP encoding.
3. **Concrete four-round witness:** `results/probability_trail.json` contains 4 round-by-round difference states. `results/standard_key_witness.json` contains a corresponding 128-bit master key, plaintext/ciphertext pair, whitening key and 3 round keys. Test `test_standard_key_witness_replays_full_trail` independently reruns the reduced four-round core, checking all internal differences.

## Fast validation without Gurobi

Run in the repository root:

```bash
python build_models.py
python -u -m unittest discover -s tests -v
python audit_solution.py
python independent_checks.py
python verify_trail.py --report results/probability_trail.json --samples 0
```

The expected outputs are `14 tests ... OK`, 97 allowed/927 disallowed ten-bit S-box patterns, 2,049 constraints, 1,280 binary variables, objective 16, and 0 constraint violations.

## Independent MILP optimization (requires Gurobi)

Python 3.10+ and Gurobi 13 are recommended. Install solver bindings separately via `python -m pip install -r requirements-gurobi.txt`, and configure an appropriate license. The original 2,049-row LP requires a non-size-limited license; `python compact_models.py` eliminates 192 redundant equality rows, leaving an equivalent 1,857-row, 1,088-variable LP compatible with the size limits of the pip license. Re-run using `python solve.py --objective active --compact --time-limit 7200`. The archived solver `OPTIMAL` record corresponds to the uncompressed model, not a new compact-model optimization.

```bash
python build_models.py
python overnight_midori_active.py --hours 6
```

The last command configures `MIPFocus=3`, `MIPGap=0`, `Threads=8`, `NodefileStart=0.5GB`. It uses a previously saved feasible solution as a MIP start if one is present; if not, the script falls back to the archived primal solution. Differences in hardware, solver version, parallelism and MIP starts may lead to substantial run-time differences. A time limit is not a guarantee of reaching `OPTIMAL`.

The probability-only search is a *separate* MILP solve and need not finish for the reported characteristic bound: once the active-S-box minimum is certified to be 16, every active S-box transition contributes at least 2 to negative log2 cost. The concrete trail attains cost 32, hence maximizes the local-DDT product among modeled single characteristics. This is not a result about exact fixed-key, aggregate differential probability.

## Four-round experiment setup

- S-box: Midori-64 4-bit Sb0.
- State encoding: sixteen hex nibbles, column-major.
- Four SubCell stages and three ShuffleCell/MixColumn inter-round stages.
- Three standard K0/K1-derived intermediate round keys with first three constants alpha0, alpha1, alpha2, together with pre/post whitening.
- **This is deliberately a four-round reduced cipher**, not a complete 16-round encryption implementation.
- Both plaintexts use the *same* master key; intermediate differences follow the recorded trail exactly.

## Provenance details

The files `results/active_optimality.json`, `results/active_solution.sol` and `results/active_solver.log` preserve a completed Gurobi experiment dated 2026-10-20. Machine-specific filesystem paths were normalized in publicly distributed text. The probability trail JSON preserves the earlier, interrupted probability search. The archived standard-key witness was preserved from a prior feasibility
result and independently checked by the replay test. For a fresh solver-generated
JSON artifact, execute `find_midori_standard_key_pair.py` against the checked-in
trail; the model may return another equally valid master-key/plaintext pair. Additional independent checks and exploratory fixed-key sampling were added subsequently; the latter are distinct from the original optimization records.

The original SageMath derivation was performed separately. `generate_hrep_sage.py` contains the optional derivation from the convex hull and set-cover MILP; depending on ties in the set-cover optimum it may return a different equivalent inequality subset. The checked-in `selected_inequalities.txt` is the **fixed** inequality set corresponding to the preserved computational results.
