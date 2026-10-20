# Reproducible Four-Round Midori-64 Differential Characteristic Analysis

**Research software and technical report** — DDT-aware MILP modeling, solver and word-level optimality checks, concrete software verification, and exploratory fixed-key measurements.

> **Scope.** This is a technical implementation and reproduction of known four-round Midori-64 bounds, not a newly published attack or a claim about full 16-round Midori-64 security. Results concern *single differential characteristics* in the specified four-round model.

## Results

| Experiment | Result | Evidence |
| --- | --- | --- |
| Minimum active S-boxes, four SubCell rounds | **16** | Gurobi `OPTIMAL`, objective and best bound 16, zero gap |
| Maximum DDT-product probability of a single characteristic | **2^-32** | Analytical consequence of the certified 16-S-box bound plus explicit attaining trail |
| Concrete path realization | **PASS** | Master key and plaintext pair for reduced four-round Midori64 with normal key schedule and whitening |
| DDT representation | **97 accepted of 1024 10-bit vectors** | Exhaustive check of 19 inequalities per S-box |

The archived active-S-box certification used **2,049 constraints**, **1,280 binary variables**, Gurobi 13.0.3, and took **794.455 seconds** on a Windows AMD Ryzen AI 7 350 machine (8 threads). The probability-only optimization was interrupted, with best primal objective 32 and dual bound 18. **We do not claim that interrupted run proved optimality.** The characteristic optimum instead follows from the independently completed active-S-box optimization and the Midori S-box's 1/4 maximum active transition probability.

## What is new in this repository (versus the published result)

The **16 active S-box** threshold is a known four-round Midori-64 benchmark
with a wide-trail explanation (Daemen--Rijmen, *The Design of Rijndael*, 2002,
Four-Round Propagation Theorem). Midori's MixColumn has branch number 4,
and its optimal ShuffleCell arrangement gives the four-round bound of 16.
This value is not a discovery here.
Besides the archived Gurobi run, `independent_checks.py` recomputes a **word-level
lower bound** by enumerating the exact activity transitions of a four-nibble
MixColumn and applying dynamic programming over 16-bit activity patterns.
It returns `[1, 4, 7, 16, 23, 30, 35, 38]` for 1–8 rounds. Since this
word-level model permits at least every bit-level differential trail, its
four-round lower bound of 16, together with the explicit 16-active trail,
provides a separate combinatorial cross-check of optimality.

The report distinguishes the DDT-product `2^-32` from the **fixed-key
occurrence rate** of that same path. A conditional Monte Carlo run for one
compatible master key estimates the latter near `2^-22`; details, sample counts,
and limitations are in [fixed-key exploration](docs/FIXED_KEY_EXPLORATION.md).
This phenomenon is not claimed as novel and should be viewed in the context
of Canteaut–Fruchon (ASIACRYPT 2025).

## Model

A state contains 16 four-bit cells in column-major order. The model contains four SubCell layers and three intervening ShuffleCell/MixColumn layers, with a nonzero initial difference. Both objective functions use the *same* DDT-aware constraints and only change the objective:

- **Activity:** seven classical inequalities per S-box (448 total).
- **Exact probability:** the Sasaki–Todo reduced H-representation of input/output differences and probability bits (19 × 64 = 1,216 constraints).
- **ShuffleCell:** 192 binary variable equalities over three transitions.
- **MixColumn:** 192 exact XOR parity equalities with binary helper variables.
- **Nonzero input:** one inequality.

The encoded Midori S-box admits active transition counts 2 or 4 out of 16, giving costs 3 or 2 respectively; inactive transitions cost 0. Thus every valid four-round characteristic with N active S-boxes has DDT-product probability at most `(1/4)^N`. With the model-certified N ≥ 16 and an explicit feasible trail with N = 16 and cost 32, the maximum single-trail DDT-product is `2^-32` **for this model**.

## Quick start

Python 3.10+ recommended. Install the small non-solver dependency; all validation/model-generation tests work **without Gurobi**:

```bash
python -m pip install -r requirements.txt
python build_models.py
python compact_models.py
python -u -m unittest discover -s tests -v
python independent_checks.py
python audit_solution.py
python verify_trail.py --report results/probability_trail.json --samples 0
```

The optional `generate_hrep_sage.py` script reconstructs the 10-dimensional H-representation and set-cover reduction in SageMath; the committed 19 inequalities are the fixed input to the reproduced experiments.

The build command writes `generated/midori_4round_active.lp` and `generated/midori_4round_probability.lp`. The optional compact generator also writes `_compact.lp` versions with the 192 redundant ShuffleCell intermediates substituted away. This leaves **1,857 rows and 1,088 binary variables**, and the recorded active-S-box assignment was checked against all remaining constraints. The audit independently checks all 1,280 saved binary assignments against all 2,049 regenerated MILP constraints and recomputes the 16-active-S-box objective.

To **rerun optimizations**, install Gurobi separately and configure a suitable license:

```bash
python -m pip install -r requirements-gurobi.txt
```

Gurobi's size-limited pip license allows at most 2,000 rows and 2,000 variables;
the original 2,049-row model exceeds that limit, while the mathematically
equivalent 1,857-row compact model is within it. Both are included for
reproducibility. The following commands use the original/full model:

```bash
python solve.py --objective active --time-limit 7200
python solve.py --objective probability --time-limit 7200
python overnight_midori_active.py --hours 6
```

To try the license-friendly compact version instead:

```bash
python solve.py --objective active --compact --time-limit 7200
```

The saved `OPTIMAL` log belongs to the **original** model; equivalent
constraints and a projected feasible solution have been checked for the
compact model, but a new optimization of it is not claimed here.

The overnight script sets `MIPFocus=3`, `Threads=8`, and uses optional `generated/probability.sol` as a MIP start; when absent it falls back to the archived feasible `results/active_solution.sol`. These are experimental settings; run time is not deterministic across machines or licenses. Results may require several hours, and early termination must not be confused with optimality.

To replay the recorded probability trail using the independent reference core,
without a probabilistically uninformative rare-event Monte Carlo loop:

```bash
python verify_trail.py --report results/probability_trail.json --samples 0
```

The **independent_checks.py** script also validates full 16-round Midori64
against the two official test vectors and reproduces the four-round witness.

To search for a concrete realization of this particular trail using the reduced four-round Midori64 key schedule:

```bash
python find_midori_standard_key_pair.py --report results/probability_trail.json --search-master-key --time-limit 120
```

That search **can** find a master key/plaintext pair, but is not guaranteed to select the same one on every run. A fixed witness is recorded in `results/standard_key_witness.json` and replayed in the automated unit tests.

## Four-round characteristic

| Round | S-box input Δ (hex) | S-box output Δ (hex) | Active | DDT-product |
| --- | --- | --- | ---: | ---: |
| 1 | `2000000000200002` | `4000000000400004` | 3 | `2^-6` |
| 2 | `0040000000000000` | `0020000000000000` | 1 | `2^-2` |
| 3 | `0000000000002202` | `0000000000004404` | 3 | `2^-6` |
| 4 | `4440000044044044` | `2220000022022022` | 9 | `2^-18` |
| **Total** | | | **16** | **`2^-32`** |

The saved `results/standard_key_witness.json` contains an explicit 128-bit master key and plaintext/ciphertext pair for the chosen **four-round truncated** Midori construction (three linear layers, standard K0/K1-derived round keys, and whitening). The corresponding test replays all 4 rounds and checks every internal difference. **One realized pair does not experimentally establish a fixed-key probability of `2^-32`; that is the conventional DDT-product trail weight.**

## Repository guide

| File/folder | Purpose |
| --- | --- |
| `docs/midori64_technical_report.pdf`, `.tex` | Standalone technical report |
| `build_models.py`, `compact_models.py`, `selected_inequalities.txt`, `generate_hrep_sage.py` | Four-round MILP generator, equivalent compact LP, saved DDT inequalities, optional SageMath derivation |
| `solve.py`, `overnight_midori_active.py` | Solver entry points and reproducible Gurobi options |
| `midori_reference.py`, `verify_trail.py`, `independent_checks.py` | Independent software-layer, published-vector, and word-level lower-bound checks |
| `find_midori_standard_key_pair.py` | Key/plaintext feasibility search under the 4-round key schedule |
| `find_midori_full_path_pair.py` | Optional simplified zero-intermediate-key four-round path witness search; **not** the standard key-schedule experiment |
| `audit_solution.py` | Constraint audit of recorded optimal primal assignment, no Gurobi required |
| `results/` | Saved optimality log, solutions, DDT trail, standard-key witness and fixed-key pilot sampling |
| `experiments/fixed_key_trail.c` | Optional conditional Monte Carlo tool for fixed-key path occurrence |
| `docs/FIXED_KEY_EXPLORATION.md` | Mathematical condition, measured results and caveats for fixed-key phenomenon |
| `tests/` | Exhaustive DDT tests, randomized trail tests, saved-result replay |
| `generated/` | Locally regenerated models and solver outputs (not versioned) |

## Reproducibility and limitations

- The solver's `OPTIMAL` log is a *computational* optimality certification, subject to MILP tolerances and correct model specification; the saved `.sol` file by itself certifies only feasibility, not a lower bound.
- The independent witness search optimizes a zero objective: its `OPTIMAL` status means a feasible pair was found, not that a global trail bound was proven.
- The model proves a result for **four SubCell layers** with **three** intervening linear layers. It is not a 16-round Midori security evaluation.
- Statistical full-differential probabilities, summed across distinct paths, may differ from a single DDT-product characteristic weight.
- The archived solver log and result JSON are preserved copies of the completed optimization, with machine-specific paths normalized. The archived witness is independently verified by the included tests; re-running its feasibility search may produce a different witness.

## Literature and attribution

The model is adapted from the published approach in M. B. İlter and A. A. Selçuk, *MILP-Aided Cryptanalysis of the FUTURE Block Cipher*, **SecITC 2022 / LNCS 13809 (2023), pp. 153–167** ([DOI](https://doi.org/10.1007/978-3-031-32636-3_9)). We adapt the general S-box, permutation and XOR formulation to **Midori**; we do not reuse FUTURE's S-box or state transformation.

- S. Banik et al., *Midori: A Block Cipher for Low Energy*, **ASIACRYPT 2015** ([preprint](https://eprint.iacr.org/2015/1142)). The four-round active-S-box count reproduces the designers' benchmark.
- Y. Sasaki and Y. Todo, *New Algorithm for Modeling S-box in MILP Based Differential and Division Trail Search*, **SecITC 2017**, pp. 150–165 ([DOI](https://doi.org/10.1007/978-3-319-69284-5_11)).

Additional foundational and recent context:

- J. Daemen and V. Rijmen, *The Design of Rijndael*, Springer, 2002, Ch. 9, Four-Round Propagation Theorem ([Springer](https://doi.org/10.1007/978-3-662-04722-4)).

- N. Mouha, Q. Wang, D. Gu and B. Preneel, *Differential and Linear Cryptanalysis Using Mixed-Integer Linear Programming*, Inscrypt 2011 / LNCS 7537 (2012), pp. 57–76 ([DOI](https://doi.org/10.1007/978-3-642-34704-7_5)).
- S. Sun et al., *Automatic Security Evaluation and (Related-key) Differential Characteristic Search*, ASIACRYPT 2014 ([preprint](https://eprint.iacr.org/2013/676)).
- A. Canteaut and M. Fruchon, *Understanding Unexpected Fixed-Key Differential Behaviours*, ASIACRYPT 2025 ([ePrint](https://eprint.iacr.org/2025/1553)).

**Status:** Standalone reproducible four-round technical project and exploratory
fixed-key experiment; not peer-reviewed and not claimed to introduce a new
cryptanalytic bound.
