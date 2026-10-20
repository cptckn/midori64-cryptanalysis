"""Build exact four-round Midori-64 MILP models (active S-boxes / probability).

Modeling follows Ilter--Selcuk, SecITC 2022 (LNCS 13809, 2023), Sections 4.1--4.4:
  - S-box activity: seven standard inequalities;
  - exact differential transitions and probability: binary H-representation
    reduced using Sasaki--Todo set cover (19 provided inequalities);
  - ShuffleCell: binary equality permutation;
  - MixColumn: XOR parity equalities with two binary dummy variables.

The three transitions between four S-box layers are modeled; a terminal
unconstrained linear layer is intentionally omitted. See README.

Requires Python standard library. Gurobi is only required for solving LPs.
"""
from __future__ import annotations

import argparse
from collections import Counter
from itertools import product
from pathlib import Path
import re

SBOX = (0xC, 0xA, 0xD, 0x3, 0xE, 0xB, 0xF, 0x7,
        0x8, 0x9, 0x1, 0x5, 0x0, 0x2, 0x4, 0x6)
SHUFFLE = (0, 10, 5, 15, 14, 4, 11, 1, 9, 3, 12, 6, 7, 13, 2, 8)
VARIABLE_ORDER = ("x3", "x2", "x1", "x0", "y3", "y2", "y1", "y0", "P0", "P1")
ROUNDS = 4
CELLS = 16
BITS = 4


def difference_distribution_table():
    table = [[0] * 16 for _ in range(16)]
    for dx in range(16):
        for x in range(16):
            table[dx][SBOX[x] ^ SBOX[x ^ dx]] += 1
    return table


def fourbits(n):
    return tuple((n >> bit) & 1 for bit in (3, 2, 1, 0))


def allowed_vectors():
    ddt = difference_distribution_table()
    probability_encoding = {16: (0, 0), 4: (0, 1), 2: (1, 1)}
    return {fourbits(dx) + fourbits(dy) + probability_encoding[ddt[dx][dy]]
            for dx in range(16) for dy in range(16) if ddt[dx][dy]}


def read_hrep(path: Path):
    """Read the verified 19 H-representation inequalities as (coeffs, rhs)."""
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith("Inequality "):
            continue
        _, equation = line.split(":", 1)
        lhs, rhs = equation.split(">=", 1)
        rhs = int(rhs.strip())
        coeffs = [0] * len(VARIABLE_ORDER)
        tokens = list(re.finditer(r"([+-]?)\s*(\d*)\s*(x[0-3]|y[0-3]|P[01])", lhs))
        assert tokens, line
        for t in tokens:
            sign, magnitude, variable = t.groups()
            idx = VARIABLE_ORDER.index(variable)
            assert coeffs[idx] == 0, f"Repeated variable {variable}"
            coeffs[idx] = (-1 if sign == "-" else 1) * int(magnitude or 1)
        # Refuse silently lost or unparsed tokens.
        remainder = re.sub(r"([+-]?)\s*(\d*)\s*(x[0-3]|y[0-3]|P[01])", "", lhs)
        assert not remainder.strip(), f"Unparsed text: {remainder!r}"
        rows.append((tuple(coeffs), rhs))
    assert len(rows) == 19, f"Expected 19 inequalities, got {len(rows)}"
    return rows


def verify_hrep(rows):
    """Exhaustively validate against all 2**10 binary patterns."""
    allowed = allowed_vectors()
    assert len(allowed) == 97
    accepted = set()
    for pt in product((0, 1), repeat=10):
        if all(sum(a * b for a, b in zip(weights, pt)) >= rhs
               for weights, rhs in rows):
            accepted.add(pt)
    assert accepted == allowed, (
        f"Invalid: extra={len(accepted-allowed)}, missing={len(allowed-accepted)}")


def verify_xor3():
    """Exhaustively validate the paper-style XOR encoding."""
    for a, b, c, out in product((0, 1), repeat=4):
        has_dummies = any(a+b+c+out -4*d1+2*d2 == 0
                          for d1, d2 in product((0, 1), repeat=2))
        assert has_dummies == (out == (a ^ b ^ c))


def lp_expression(terms):
    """Convert mapping variable -> integer coefficient to LP syntax."""
    chunks = []
    for name, coef in terms.items():
        if not coef:
            continue
        mag = "" if abs(coef) == 1 else f"{abs(coef)} "
        chunks.append(f"{'+' if coef > 0 else '-'} {mag}{name}")
    assert chunks, "Empty LP expression"
    return (" ".join(chunks))[2:] if chunks[0].startswith("+ ") else " ".join(chunks)


def x(r, c, b): return f"x_{r}_{c}_{b}"
def y(r, c, b): return f"y_{r}_{c}_{b}"
def z(r, c, b): return f"z_{r}_{c}_{b}"
def a(r, c): return f"A_{r}_{c}"
def pp(i, r, c): return f"P{i}_{r}_{c}"
def d(r, col, row, bit, i): return f"mc_r{r}_col{col}_row{row}_b{bit}_d{i}"


def model_text(objective, hrep):
    assert objective in {"active", "probability"}
    output = ["\\ Model built from Ilter--Selcuk 2022: exact S-box H-rep and binary XOR", "Minimize"]
    obj_vars = ([a(r, c) for r in range(ROUNDS) for c in range(CELLS)]
                if objective == "active" else
                [pp(0, r, c) for r in range(ROUNDS) for c in range(CELLS)] +
                [pp(1, r, c) for r in range(ROUNDS) for c in range(CELLS)])
    if objective == "active":
        obj = {name: 1 for name in obj_vars}
    else:
        obj = {name: (1 if name.startswith("P0") else 2) for name in obj_vars}
    output.append(" obj: " + lp_expression(obj))
    output.append("Subject To")
    counts = Counter()

    def add(label, terms, op, rhs, category):
        output.append(f" {label}: {lp_expression(terms)} {op} {rhs}")
        counts[category] += 1

    for r in range(ROUNDS):
        for c in range(CELLS):
            # Seven activity/bijectivity inequalities (Ilter--Selcuk, Sec. 4.1).
            for b in range(BITS):
                add(f"act_r{r}_c{c}_b{b}", {a(r,c):1,x(r,c,b):-1}, ">=", 0, "activity")
            add(f"act_r{r}_c{c}_lower", {**{x(r,c,b):1 for b in range(BITS)}, a(r,c):-1}, ">=", 0, "activity")
            add(f"act_r{r}_c{c}_xy", {**{x(r,c,b):4 for b in range(BITS)},
                    **{y(r,c,b):-1 for b in range(BITS)}}, ">=", 0, "activity")
            add(f"act_r{r}_c{c}_yx", {**{y(r,c,b):4 for b in range(BITS)},
                    **{x(r,c,b):-1 for b in range(BITS)}}, ">=", 0, "activity")
            # Exact DDT and probability bits, reused under either objective.
            vars10 = [x(r,c,b) for b in (3,2,1,0)] + \
                     [y(r,c,b) for b in (3,2,1,0)] + [pp(0,r,c),pp(1,r,c)]
            for j, (weights, rhs) in enumerate(hrep):
                add(f"hrep_r{r}_c{c}_i{j}", dict(zip(vars10, weights)), ">=", rhs, "hrep")

    for r in range(ROUNDS-1):
        # ShuffleCell permutation equalities (Sec. 4.2).
        for out_cell, in_cell in enumerate(SHUFFLE):
            for b in range(BITS):
                add(f"sh_r{r}_c{out_cell}_b{b}",
                    {z(r,out_cell,b):1, y(r,in_cell,b):-1}, "=", 0, "shuffle")
        # XOR/parity encoding (Sec. 4.4).
        for col in range(4):
            cells = tuple(4*col+i for i in range(4))
            for row in range(4):
                for b in range(BITS):
                    terms = {z(r, cells[j], b):1 for j in range(4) if j != row}
                    terms[x(r+1,cells[row],b)] = 1
                    terms[d(r,col,row,b,1)] = -4
                    terms[d(r,col,row,b,2)] = 2
                    add(f"mix_r{r}_col{col}_row{row}_b{b}", terms, "=", 0, "mix")

    add("nonzero_input", {x(0,c,b):1 for c in range(CELLS) for b in range(BITS)},
        ">=", 1, "nonzero")

    binary = set()
    for r in range(ROUNDS):
        for c in range(CELLS):
            binary.add(a(r,c))
            for i in range(2):binary.add(pp(i,r,c))
            for b in range(BITS):
                binary.add(x(r,c,b));binary.add(y(r,c,b))
    for r in range(ROUNDS-1):
        for c in range(CELLS):
            for b in range(BITS):binary.add(z(r,c,b))
        for col in range(4):
            for row in range(4):
                for b in range(BITS):
                    for i in (1,2):binary.add(d(r,col,row,b,i))
    output.append("Binary")
    output.extend(" " + name for name in sorted(binary))
    output.append("End")
    assert sum(counts.values()) == 2049, counts
    assert len(binary) == 1280, len(binary)
    return "\n".join(output) + "\n", counts, len(binary)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir",type=Path,default=Path(__file__).resolve().parent/"generated")
    args=parser.parse_args()
    root=Path(__file__).resolve().parent
    hrep=read_hrep(root/"selected_inequalities.txt")
    verify_hrep(hrep)
    verify_xor3()
    assert sorted(SHUFFLE)==list(range(16))
    args.output_dir.mkdir(parents=True,exist_ok=True)
    print("Independent checks passed: 97 allowed/927 forbidden DDT vectors; 3-input XOR exact.")
    for mode in ("active", "probability"):
        text, counts, nvars=model_text(mode,hrep)
        filepath=args.output_dir/f"midori_4round_{mode}.lp"
        filepath.write_text(text,encoding="utf-8")
        print(f"{filepath}: {sum(counts.values())} constraints, {nvars} binary variables; breakdown {dict(counts)}")

if __name__ == "__main__":
    main()
