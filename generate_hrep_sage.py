"""Optional SageMath derivation of the reduced Midori64 DDT H-representation.

Requires SageMath. May select a different but equivalent 19-inequality subset
when optimal set cover has several solutions. The checked-in inequalities
are the authoritative reproducible model input.
"""

from itertools import product
from sage.all import QQ, ZZ, lcm, Polyhedron
from sage.numerical.mip import MixedIntegerLinearProgram


s_box = {
    0x0: 0xc,
    0x1: 0xa,
    0x2: 0xd,
    0x3: 0x3,
    0x4: 0xe,
    0x5: 0xb,
    0x6: 0xf,
    0x7: 0x7,
    0x8: 0x8,
    0x9: 0x9,
    0xa: 0x1,
    0xb: 0x5,
    0xc: 0x0,
    0xd: 0x2,
    0xe: 0x4,
    0xf: 0x6,
}



def bits4(n):

    return [
        (n >> 3) & 1,
        (n >> 2) & 1,
        (n >> 1) & 1,
        n & 1,
    ]


def lp_term(coeff, var):
    coeff = int(coeff)

    if coeff == 0:
        return None
    elif coeff == 1:
        return f"+ {var}"
    elif coeff == -1:
        return f"- {var}"
    elif coeff > 0:
        return f"+ {coeff} {var}"
    else:
        return f"- {-coeff} {var}"


def integerize_inequality(A, b):

    values = list(A) + [b]
    denominators = [QQ(v).denominator() for v in values]
    scale = lcm(denominators)

    A_int = [ZZ(scale * QQ(a)) for a in A]
    b_int = ZZ(scale * QQ(b))

    return A_int, b_int


ddt = [[0 for _ in range(16)] for _ in range(16)]

for dx in range(16):
    for x in range(16):
        x2 = x ^ dx
        dy = s_box[x] ^ s_box[x2]
        ddt[dx][dy] += 1


nonzero_count = 0

for dx in range(16):
    for dy in range(16):
        if ddt[dx][dy] != 0:
            nonzero_count += 1

print("DDT nonzero entries:", nonzero_count)


allowed_vectors = []

for dx in range(16):
    for dy in range(16):
        count = ddt[dx][dy]

        if count == 0:
            continue

        if count == 16:
            probability_bits = [0, 0]      # probability = 1, cost = 0
        elif count == 4:
            probability_bits = [0, 1]      # probability = 2^-2, cost = 2
        elif count == 2:
            probability_bits = [1, 1]      # probability = 2^-3, cost = 3
        else:
            raise ValueError(
                f"Unexpected DDT count {count} for dx={dx:x}, dy={dy:x}"
            )

        vector = bits4(dx) + bits4(dy) + probability_bits
        allowed_vectors.append(vector)


allowed_set = {tuple(v) for v in allowed_vectors}

print("Allowed DDT-probability vectors:", len(allowed_vectors))


invalid_binary_points = []

for point in product([0, 1], repeat=10):
    if tuple(point) not in allowed_set:
        invalid_binary_points.append(list(point))

print("Invalid binary points:", len(invalid_binary_points))


P = Polyhedron(vertices=allowed_vectors)
hrep = P.Hrepresentation()

inequalities = []

for h in hrep:
    b = h.b()
    A = list(h.A())

    A_int, b_int = integerize_inequality(A, b)
    inequalities.append((A_int, b_int))

print("Initial H-representation inequalities:", len(inequalities))


covering_inequalities = []

for point in invalid_binary_points:
    covering = []

    for j, (A, b) in enumerate(inequalities):
        value = b + sum(A[i] * point[i] for i in range(10))

        if value < 0:
            covering.append(j)

    if len(covering) == 0:
        raise ValueError(f"No inequality blocks invalid point {point}")

    covering_inequalities.append(covering)



mip = MixedIntegerLinearProgram(maximization=False)

s = mip.new_variable(binary=True)

mip.set_objective(sum(s[j] for j in range(len(inequalities))))

for covering in covering_inequalities:
    mip.add_constraint(sum(s[j] for j in covering) >= 1)

opt_value = mip.solve()

selected_indices = []

for j in range(len(inequalities)):
    if mip.get_values(s[j]) > 0.5:
        selected_indices.append(j)

selected_inequalities = [inequalities[j] for j in selected_indices]

print("Selected inequalities after reduction:", len(selected_inequalities))

variable_names = [
    "x3", "x2", "x1", "x0",
    "y3", "y2", "y1", "y0",
    "P0", "P1"
]

with open("selected_probability_inequalities.txt", "w") as file:
    file.write("Selected inequalities after Sasaki-Todo style reduction\n")
    file.write("Variable order: E = (x3,x2,x1,x0,y3,y2,y1,y0,P0,P1)\n\n")

    for idx, (A, b) in enumerate(selected_inequalities):
        terms = []

        for coeff, var in zip(A, variable_names):
            term = lp_term(coeff, var)

            if term is not None:
                terms.append(term)

        left_side = " ".join(terms)

        if left_side.startswith("+ "):
            left_side = left_side[2:]

        file.write(f"Inequality {idx + 1}: {left_side} >= {-b}\n")


for point in allowed_vectors:
    for A, b in selected_inequalities:
        value = b + sum(A[i] * point[i] for i in range(10))

        if value < 0:
            raise ValueError(f"Selected inequality cuts allowed point {point}")


for point in invalid_binary_points:
    blocked = False

    for A, b in selected_inequalities:
        value = b + sum(A[i] * point[i] for i in range(10))

        if value < 0:
            blocked = True
            break

    if not blocked:
        raise ValueError(f"Invalid point is not blocked: {point}")


print("Verification passed.")


