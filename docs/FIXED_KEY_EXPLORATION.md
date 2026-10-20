# Fixed-key behaviour of the four-round trail (exploratory)

## What is established

The characteristic `2000000000200002 -> ... -> 2220000022022022`
contains 16 active S-boxes and has a **local DDT-product weight** of `2^-32`.
That number is *not* a probability measured under one particular master key.
`independent_checks.py` confirms the trail and independently replays the
standard-key witness against a full-Midori implementation validated by two
published test vectors.

For Midori's Sb0, the left-input conforming sets are

- `2 -> 4`: `{C,D,E,F}`, with outputs `{0,2,4,6}`;
- `4 -> 2`: `{0,2,4,6}`, with outputs `{C,D,E,F}`.

Both sets are two-dimensional affine subspaces of $\mathbb{F}_2^4$ (viewing nibbles as four-bit vectors). In particular, for the recorded
trail, after the first round the active input at position 2 of round 2 is

`x1[2] = S(x0[0]) xor S(x0[10]) xor S(x0[15]) xor RK0[2]`.

The first three S-box outputs belong to `E = {0,2,4,6}`, a linear subgroup
under nibble XOR. Consequently, the round-2 transition `4 -> 2` is possible
for a first-round conforming triple **if and only if** `RK0[2]` belongs to E.
With a uniformly chosen 128-bit master key, RK0[2] is uniform, so this single
round-2 condition rules out three quarters of master keys. This is a *necessary*
condition, not a complete description of four-round compatible keys.

The master key in `results/standard_key_witness.json` is compatible;
Gurobi found a full path for this key, and the independent reference replay
verifies every intermediate difference. A separate fixed-zero-key Gurobi
feasibility model returned `INFEASIBLE` for this **specific path**.

## Conditional Monte Carlo experiment

Compile using GCC or Clang:

```sh
gcc -O3 -o fixed_key_trail experiments/fixed_key_trail.c -lm
./fixed_key_trail key 00000101000100016177166671066660 100000000 2026
./fixed_key_trail key 00000101000100016177166671066660 100000000 31415
```

The program samples a uniform internal state after pre-whitening, with the
three active first-round nibbles forced into one of the conforming triples
for rounds 1 and 2. It corrects for the exact probability `2^-6` of that
conditioning. For each seed, 100 million **conditional** trials yielded:

| Seed | Full-trail hits | Estimated fixed-key path probability |
|---|---:|---:|
| `2026` | 1,541 | approximately `2^-21.99` |
| `31415` | 1,543 | approximately `2^-21.98` |

The two runs are consistent with a fixed-key probability near `2^-22`,
substantially above the DDT-product `2^-32`. The values are **estimates**;
we have not proved an exact probability or the fraction of compatible keys.
The raw outputs and estimation method are recorded in
`results/fixed_key_sampling.json`. Simulations over many keys require far more
conditional samples per key than small pilot studies: *no observed hits does
not establish that a key is impossible*.

## Scientific context and next questions

This behaviour should be discussed in the light of:

- A. Canteaut and M. Fruchon, *Understanding Unexpected Fixed-Key Differential
  Behaviours: How to Avoid Major Weaknesses in Lightweight Designs*,
  ASIACRYPT 2025, https://eprint.iacr.org/2025/1553.
- T. Beyne and V. Rijmen, *Differential Cryptanalysis in the Fixed-Key Model*,
  https://eprint.iacr.org/2022/837.

The underlying fixed-key phenomena are known in the literature; this repository
reproduces and illustrates them for a specified four-round trail. A promising
extension is to derive exact linear conditions on the key and then count
compatible initial states algebraically, instead of estimating probabilities.
There is **no novelty claim** about fixed-key clustering or the original
four-round activity bound.
