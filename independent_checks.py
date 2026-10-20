"""Independent cross-checks for the four-round Midori-64 results (no MILP solver needed).

1. A full 16-round Midori-64, written from the specification and validated
   against both official test vectors (Banik et al. 2015, Appendix B).
2. DDT statistics of Sb0.
3. The four-round characteristic: DDT-valid transitions + exact linear layers.
4. The master-key witness, replayed with the validated implementation
   (reduced to four SubCell layers, with whitening and real round constants).
5. A word-level lower bound on active S-boxes: exact column transition sets of
   MixColumn + dynamic programming over 16-bit activity patterns. This is a
   relaxation of the bit-level model, so (lower bound 16) + (explicit 16-trail)
   proves min = 16 without trusting a MILP solver's tolerances.

Requires numpy. Runs in about a second.
"""
import collections
import itertools

import numpy as np

SB0 = [0xC, 0xA, 0xD, 0x3, 0xE, 0xB, 0xF, 0x7, 0x8, 0x9, 0x1, 0x5, 0x0, 0x2, 0x4, 0x6]
PERM = [0, 10, 5, 15, 14, 4, 11, 1, 9, 3, 12, 6, 7, 13, 2, 8]   # ShuffleCell: out[i] = in[PERM[i]]
# Round constants beta_0..beta_14, 16 bits each, column-major cell order (added to each cell's LSB).
BETA = [
    [0,0,0,1,0,1,0,1,1,0,1,1,0,0,1,1], [0,1,1,1,1,0,0,0,1,1,0,0,0,0,0,0],
    [1,0,1,0,0,1,0,0,0,0,1,1,0,1,0,1], [0,1,1,0,0,0,1,0,0,0,0,1,0,0,1,1],
    [0,0,0,1,0,0,0,0,0,1,0,0,1,1,1,1], [1,1,0,1,0,0,0,1,0,1,1,1,0,0,0,0],
    [0,0,0,0,0,0,1,0,0,1,1,0,0,1,1,0], [0,0,0,0,1,0,1,1,1,1,0,0,1,1,0,0],
    [1,0,0,1,0,1,0,0,1,0,0,0,0,0,0,1], [0,1,0,0,0,0,0,0,1,0,1,1,1,0,0,0],
    [0,1,1,1,0,0,0,1,1,0,0,1,0,1,1,1], [0,0,1,0,0,0,1,0,1,0,0,0,1,1,1,0],
    [0,1,0,1,0,0,0,1,0,0,1,1,0,0,0,0], [1,1,1,1,1,0,0,0,1,1,0,0,1,0,1,0],
    [1,1,0,1,1,1,1,1,1,0,0,1,0,0,0,0],
]


def st(h): return [int(c, 16) for c in h]
def hx(s): return ''.join('%X' % v for v in s)
def xor(a, b): return [x ^ y for x, y in zip(a, b)]
def sub(s): return [SB0[v] for v in s]
def shuf(s): return [s[PERM[i]] for i in range(16)]


def mix(s):
    out = [0] * 16
    for c in range(4):
        col = s[4*c:4*c+4]
        t = col[0] ^ col[1] ^ col[2] ^ col[3]
        for r in range(4):
            out[4*c+r] = t ^ col[r]      # row r = XOR of the other three cells
    return out


def keys(master_hex):
    k0, k1 = st(master_hex[:16]), st(master_hex[16:])
    return xor(k0, k1), [xor((k0, k1)[i % 2], BETA[i]) for i in range(15)]


def encrypt(p_hex, k_hex, sbox_layers=16, trace=None):
    """sbox_layers=16 is full Midori-64; sbox_layers=4 is the report's reduced cipher."""
    wk, rk = keys(k_hex)
    s = xor(st(p_hex), wk)
    for i in range(sbox_layers - 1):
        if trace is not None: trace.append(s)
        s = sub(s)
        if trace is not None: trace.append(s)
        s = xor(mix(shuf(s)), rk[i])
    if trace is not None: trace.append(s)
    s = sub(s)
    if trace is not None: trace.append(s)
    return hx(xor(s, wk))


def main():
    # 1. Official test vectors
    for p, k, c in [('0000000000000000', '00000000000000000000000000000000', '3C9CCEDA2BBD449A'),
                    ('42C20FD3B586879E', '687DED3B3C85B3F35B1009863E2A8CBF', '66BCDC6270D901CD')]:
        assert encrypt(p, k) == c, (p, k)
    print('[1] full Midori-64 matches both official test vectors')

    # 2. DDT
    ddt = [[0] * 16 for _ in range(16)]
    for a in range(16):
        for x in range(16):
            ddt[a][SB0[x] ^ SB0[x ^ a]] += 1
    hist = collections.Counter(v for row in ddt for v in row)
    print('[2] DDT histogram', dict(sorted(hist.items())),
          '| max active entry', max(ddt[a][b] for a in range(1, 16) for b in range(16)))

    # 3. Trail
    trail = [('2000000000200002', '4000000000400004'), ('0040000000000000', '0020000000000000'),
             ('0000000000002202', '0000000000004404'), ('4440000044044044', '2220000022022022')]
    active = cost = 0
    for r, (i, o) in enumerate(trail):
        for a, b in zip(st(i), st(o)):
            assert ddt[a][b] > 0
            if a:
                active += 1
                cost += {4: 2, 2: 3}[ddt[a][b]]
        if r < 3:
            assert mix(shuf(st(o))) == st(trail[r + 1][0])
    print(f'[3] trail valid: {active} active S-boxes, DDT-product 2^-{cost}')

    # 4. Witness
    key = '00000101000100016177166671066660'
    t0, t1 = [], []
    c0 = encrypt('ADBBDBABBDCBAAAD', key, 4, t0)
    c1 = encrypt('8DBBDBABBDEBAAAF', key, 4, t1)
    for r in range(4):
        assert hx(xor(t0[2*r], t1[2*r])) == trail[r][0]
        assert hx(xor(t0[2*r+1], t1[2*r+1])) == trail[r][1]
    assert (c0, c1) == ('ADBDDBADBDABACAD', '8F9DDBAD9FA98C8F')
    print('[4] witness pair follows every internal difference; ciphertexts match')

    # 5. Word-level DP lower bound
    trans = collections.defaultdict(set)          # column in-pattern -> possible out-patterns
    for v in itertools.product(range(16), repeat=4):
        t = v[0] ^ v[1] ^ v[2] ^ v[3]
        out = [t ^ x for x in v]
        trans[sum((x != 0) << k for k, x in enumerate(v))].add(sum((x != 0) << k for k, x in enumerate(out)))
    inf = 10**6
    idx = np.arange(1 << 16)
    weight = np.array([bin(p).count('1') for p in range(1 << 16)])
    shuffled = np.zeros(1 << 16, dtype=np.int64)
    for i in range(16):
        shuffled |= ((idx >> PERM[i]) & 1) << i
    f = np.where(idx > 0, weight, inf)            # f[p] = min #active so far, p = current SubCell pattern
    best = [int(f.min())]
    for _ in range(7):
        g = np.full(1 << 16, inf)
        g[shuffled] = f
        A = g.reshape(16, 16, 16, 16)
        for ax in range(4):                       # apply column transitions one column at a time
            Bn = np.full_like(A, inf)
            zero = [slice(None)] * 4; zero[ax] = 0
            Bn[tuple(zero)] = np.take(A, 0, axis=ax)
            for s in range(1, 16):
                src = np.take(A, s, axis=ax)
                for q in trans[s]:
                    sl = [slice(None)] * 4; sl[ax] = q
                    Bn[tuple(sl)] = np.minimum(Bn[tuple(sl)], src)
            A = Bn
        f = A.reshape(-1) + weight
        f[0] = inf
        best.append(int(f.min()))
    print('[5] word-level minimum active S-boxes for 1..8 rounds:', best,
          '(designers, Table 8: 4..8 rounds = 16, 23, 30, 35, 38)')


if __name__ == '__main__':
    main()
