"""Search a *specific saved MILP trail* under Midori64's real 4-round key schedule.

Usage (from folder containing midori_reference.py and generated/probability_result.json):
    py -u find_midori_standard_key_pair.py --master-key 00000000000000000000000000000000 --time-limit 120
    py -u find_midori_standard_key_pair.py --search-master-key --time-limit 120

Midori64 uses 128-bit K=K0||K1, WK=K0 XOR K1, RKi=K_(i mod 2) XOR alpha_i.
The four-round reduced cipher uses 3 full rounds + final SubCell, then WK.
The three alpha matrices (i=0..2) are encoded here in *column-major* nibble order,
following the published Midori64 state S = [[s0,s4,s8,s12],...].

A feasible pair is verified by independently executing both plaintexts through
four SubCell rounds, the three genuine round-key additions, and two whitenings.
INFEASIBLE means this *specific* differential trail cannot be realized for the
selected key mode; it does not disprove its MILP differential-transition feasibility.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from midori_reference import (
    SBOX_REF, SHUFFLE_REF, parse_hex_state, hex_state, xor_states,
    round_trace_ref, validate_solver_trail,
)

# As published in Midori's specifications (row-major matrix notation); the
# state serialization used here is column-major, indexed s0,s1,...,s15.
ALPHA_ROWS = (
    ((0,0,1,0), (0,1,0,0), (0,0,1,1), (1,1,1,1)),
    ((0,1,1,0), (1,0,1,0), (1,0,0,0), (1,0,0,0)),
    ((1,0,0,0), (0,1,0,1), (1,0,1,0), (0,0,1,1)),
)
ALPHAS = tuple(tuple(rows[r][c] for c in range(4) for r in range(4))
               for rows in ALPHA_ROWS)
assert tuple(hex_state(a) for a in ALPHAS) == (
    '0001010110110011', '0111100011000000', '1010010000110101'
)


def split_master_key(hex32):
    hex32 = hex32.upper()
    if len(hex32) != 32:
        raise ValueError('The Midori64 master key must be 32 hexadecimal digits (128 bits)')
    return parse_hex_state(hex32[:16]), parse_hex_state(hex32[16:])


def derive_keys(k0, k1):
    wk = xor_states(k0, k1)
    rks = tuple(xor_states((k0,k1)[r%2], ALPHAS[r]) for r in range(3))
    return wk, rks


def run():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--report', type=Path, default=Path('generated/probability_result.json'))
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument('--master-key', default='0'*32, metavar='32_HEX_DIGITS')
    mode.add_argument('--search-master-key', action='store_true')
    ap.add_argument('--time-limit', type=float, default=120)
    ap.add_argument('--out', type=Path, default=Path('generated/standard_key_pair.json'))
    args = ap.parse_args()
    if args.time_limit <= 0:
        ap.error('Time limit must be positive')
    trail = json.loads(args.report.read_text(encoding='utf-8'))
    validate_solver_trail(trail)
    fixed_halves = None if args.search_master_key else split_master_key(args.master_key)

    try:
        import gurobipy as gp
        from gurobipy import GRB
    except ImportError as e:
        raise SystemExit('Install gurobipy: py -m pip install gurobipy') from e

    m = gp.Model('midori64_4round_official_keys_exact_path_feasibility')
    m.Params.TimeLimit = args.time_limit
    m.Params.MIPFocus = 1
    selected, allowed = {}, {}
    for r, rec in enumerate(trail['per_round_trail']):
        for c in range(16):
            dx = int(rec['input_differences'][c], 16)
            dy = int(rec['sbox_output_differences'][c], 16)
            vals = [x for x in range(16) if SBOX_REF[x] ^ SBOX_REF[x^dx] == dy]
            if not vals:
                raise ValueError(f'Forbidden differential at round={r}, cell={c}')
            allowed[r,c] = vals
            for x in vals:
                selected[r,c,x] = m.addVar(vtype=GRB.BINARY, name=f'v_{r}_{c}_{x:X}')
            m.addConstr(gp.quicksum(selected[r,c,x] for x in vals)==1,
                        name=f'onehot_{r}_{c}')

    # K0/K1 bit variables are reused across alternating rounds (0 and 2 use K0).
    kb = {}
    if args.search_master_key:
        for h in range(2):
            for c in range(16):
                for b in range(4):
                    kb[h,c,b] = m.addVar(vtype=GRB.BINARY, name=f'K{h}_{c}_{b}')

    def selected_bit(r,c,b,sbox):
        return gp.quicksum(((SBOX_REF[x] if sbox else x) >> b & 1)*selected[r,c,x]
                           for x in allowed[r,c])

    for r in range(3):
        for c in range(16):
            column_start = c//4*4
            contributors = [SHUFFLE_REF[j] for j in range(column_start,column_start+4) if j!=c]
            for b in range(4):
                prev_bits = gp.quicksum(selected_bit(r,j,b,True) for j in contributors)
                next_bit = selected_bit(r+1,c,b,False)
                if fixed_halves is None:
                    # Output = XOR(three source nibbles, K_r, alpha_r)
                    alpha_b = (ALPHAS[r][c] >> b) & 1
                    key_b = kb[r%2,c,b]
                    q = m.addVar(vtype=GRB.INTEGER, lb=0, ub=3,
                                 name=f'parity_{r}_{c}_{b}')
                    m.addConstr(prev_bits + next_bit + key_b + alpha_b == 2*q,
                                name=f'xor_{r}_{c}_{b}')
                else:
                    rkey = xor_states(fixed_halves[r%2], ALPHAS[r])
                    key_b = (rkey[c] >> b) & 1
                    q = m.addVar(vtype=GRB.INTEGER, lb=0, ub=3,
                                 name=f'parity_{r}_{c}_{b}')
                    m.addConstr(prev_bits + next_bit + key_b == 2*q,
                                name=f'xor_{r}_{c}_{b}')

    m.setObjective(0, GRB.MINIMIZE)
    m.optimize()
    names = {GRB.OPTIMAL: 'OPTIMAL', GRB.INFEASIBLE: 'INFEASIBLE',
             GRB.TIME_LIMIT: 'TIME_LIMIT', GRB.INTERRUPTED: 'INTERRUPTED'}
    result = {'gurobi_status': names.get(m.Status, str(m.Status)),
              'sol_count': m.SolCount,
              'key_mode': 'SEARCH_STANDARD_MASTER_KEY' if args.search_master_key else 'FIXED_STANDARD_MASTER_KEY',
              'trail_probability_from_ddt': '2^-32',
              'global_probability_optimality_proven': False}
    if m.SolCount:
        if args.search_master_key:
            halves = tuple(tuple(sum((int(kb[h,c,b].X > 0.5) << b) for b in range(4))
                                 for c in range(16)) for h in range(2))
        else:
            halves = fixed_halves
        wk, rks = derive_keys(*halves)
        internal_states = []
        for r in range(4):
            internal_states.append(tuple(max(allowed[r,c],key=lambda x:selected[r,c,x].X)
                                         for c in range(16)))
        internal_left = internal_states[0]
        delta = parse_hex_state(trail['per_round_trail'][0]['input_differences'])
        internal_right = xor_states(internal_left, delta)
        # Whitening before round 0; the actual plaintexts are WK-XOR of the
        # chosen first SubCell input states, not those internal states themselves.
        p_left, p_right = xor_states(internal_left, wk), xor_states(internal_right, wk)
        trace_l = round_trace_ref(xor_states(p_left, wk), rks)
        trace_r = round_trace_ref(xor_states(p_right, wk), rks)
        for r, expected in enumerate(trail['per_round_trail']):
            observed_in = hex_state(xor_states(trace_l[r]['input'], trace_r[r]['input']))
            observed_out = hex_state(xor_states(trace_l[r]['sbox_output'], trace_r[r]['sbox_output']))
            if observed_in != expected['input_differences'].upper() or observed_out != expected['sbox_output_differences'].upper():
                raise RuntimeError(f'Independent four-round trace mismatch at round {r}')
            if trace_l[r]['input'] != internal_states[r]:
                raise RuntimeError(f'Independent intermediate state mismatch at round {r}')
        c_left = xor_states(trace_l[-1]['sbox_output'],wk)
        c_right = xor_states(trace_r[-1]['sbox_output'],wk)
        result.update({
            'master_key': hex_state(halves[0])+hex_state(halves[1]),
            'whitening_key': hex_state(wk),
            'round_keys': [hex_state(k) for k in rks],
            'plaintext_left':hex_state(p_left), 'plaintext_right':hex_state(p_right),
            'ciphertext_left':hex_state(c_left), 'ciphertext_right':hex_state(c_right),
            'full_trail_independent_check':'PASS',
            'per_round_inputs': [hex_state(x['input']) for x in trace_l],
            'notes':'Complete 4-round reduced Midori64 with standard K0/K1, 3 alpha round constants, initial/final WK.'
        })
    else:
        result['notes']= ('INFEASIBLE excludes a pair following THIS EXACT trail for this key mode; '
                          'TIME_LIMIT/INTERRUPTED is inconclusive. Other differential trails may work.')
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))
    print('Saved report:',args.out)

if __name__=='__main__':
    run()
