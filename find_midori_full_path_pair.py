"""Find a full four-round plaintext-pair witness for a saved Midori MILP trail.

Run in repository root (requires gurobipy and midori_reference.py):
    py -u find_midori_full_path_pair.py --time-limit 600

By default, this uses the existing independent zero-key reduced-round core:
three SubCell/ShuffleCell/MixColumn transitions and final SubCell.
This is NOT a reproduction of the complete Midori-64 round-constant/key schedule.

A feasible solution is an actual 64-bit pair satisfying the entire specified
four-round differential characteristic in the selected reduced-round core.
An unknown/time-limit outcome does NOT disprove the trail.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from midori_reference import (
    SBOX_REF, SHUFFLE_REF, ZERO, parse_hex_state, hex_state,
    xor_states, round_trace_ref, validate_solver_trail,
)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--report', type=Path,
                   default=Path(__file__).resolve().parent / 'generated' / 'probability_result.json')
    p.add_argument('--time-limit', type=float, default=600.0)
    p.add_argument('--out', type=Path,
                   default=Path(__file__).resolve().parent / 'generated' / 'full_path_pair.json')
    p.add_argument('--round-keys', nargs=3, metavar='HEX16', default=['0'*16]*3,
                   help='Three explicit 16-nibble XOR additions between S-box rounds (default: zero)')
    args = p.parse_args()

    if args.time_limit <= 0:
        p.error('--time-limit must be positive')
    trail = json.loads(args.report.read_text(encoding='utf-8'))
    validate_solver_trail(trail)
    keys = tuple(parse_hex_state(s) for s in args.round_keys)

    try:
        import gurobipy as gp
        from gurobipy import GRB
    except ImportError as e:
        raise SystemExit('Install gurobipy first: py -m pip install gurobipy') from e

    m = gp.Model('midori_four_round_full_pair_feasibility')
    m.Params.TimeLimit = args.time_limit
    m.Params.MIPFocus = 1

    # v[r,c,x] selects the *left-hand* input nibble to S-box c in round r.
    # The right-hand input is necessarily x XOR the prescribed input difference.
    # Retain only selections with exactly the prescribed S-box output difference.
    selected = {}
    allowed = {}
    for r, data in enumerate(trail['per_round_trail']):
        for c in range(16):
            dx = int(data['input_differences'][c], 16)
            dy = int(data['sbox_output_differences'][c], 16)
            vals = [x for x in range(16)
                    if (SBOX_REF[x] ^ SBOX_REF[x ^ dx]) == dy]
            assert vals, f'Empty DDT support at round {r}, cell {c}'
            allowed[r,c] = vals
            for x in vals:
                selected[r,c,x] = m.addVar(vtype=GRB.BINARY, name=f'v_{r}_{c}_{x}')
            m.addConstr(gp.quicksum(selected[r,c,x] for x in vals) == 1,
                        name=f'onehot_{r}_{c}')

    def selected_bit(r, c, b, after_sbox):
        return gp.quicksum(((SBOX_REF[x] if after_sbox else x) >> b & 1)
                           * selected[r,c,x] for x in allowed[r,c])

    # Inter-round equality: x_{r+1} = MixColumn(ShuffleCell(S(x_r))) XOR key_r.
    # Each MixColumn output bit is XOR of the other 3 column entries.
    # Enforce parity exactly with an integer helper variable (0,1,2).
    for r in range(3):
        for c in range(16):
            col_start = (c // 4) * 4
            contributing_sbox_cells = [SHUFFLE_REF[j] for j in range(col_start, col_start+4)
                                        if j != c]
            for b in range(4):
                kb = (keys[r][c] >> b) & 1
                helper = m.addVar(vtype=GRB.INTEGER, lb=0, ub=2, name=f'parity_{r}_{c}_{b}')
                expr = gp.quicksum(selected_bit(r, j, b, True) for j in contributing_sbox_cells)
                expr += selected_bit(r+1, c, b, False) + kb
                m.addConstr(expr == 2*helper, name=f'xor_{r}_{c}_{b}')

    m.setObjective(0, GRB.MINIMIZE)  # Pure feasibility search, not probability optimization.
    m.optimize()

    if m.SolCount == 0:
        print('No complete pair found; solver status:', m.Status)
        print('If TIME_LIMIT, this does NOT mean the characteristic is impossible.')
        return

    states = []
    for r in range(4):
        nibbles = []
        for c in range(16):
            xs = allowed[r,c]
            x = max(xs, key=lambda q: selected[r,c,q].X)
            if selected[r,c,x].X < 0.5:
                raise RuntimeError('Gurobi returned no integral selected S-box input')
            nibbles.append(x)
        states.append(tuple(nibbles))

    left = states[0]
    delta = parse_hex_state(trail['per_round_trail'][0]['input_differences'])
    right = xor_states(left, delta)
    trace_left, trace_right = round_trace_ref(left, keys), round_trace_ref(right, keys)
    for r, wanted in enumerate(trail['per_round_trail']):
        actual_in = hex_state(xor_states(trace_left[r]['input'], trace_right[r]['input']))
        actual_out = hex_state(xor_states(trace_left[r]['sbox_output'], trace_right[r]['sbox_output']))
        if (actual_in.lower() != wanted['input_differences'].lower() or
                actual_out.lower() != wanted['sbox_output_differences'].lower() or
                trace_left[r]['input'] != states[r]):
            raise RuntimeError(f'Unexpected witness mismatch in round {r}')

    evidence = {
        'status': 'FEASIBLE_COMPLETE_FOUR_ROUND_PAIR',
        'model': 'four S-box rounds, three ShuffleCell/MixColumn layers with specified XOR round keys',
        'round_keys': [hex_state(k) for k in keys],
        'plaintext_left': hex_state(left),
        'plaintext_right': hex_state(right),
        'ciphertext_left': hex_state(trace_left[-1]['sbox_output']),
        'ciphertext_right': hex_state(trace_right[-1]['sbox_output']),
        'verified_trail': [
            {
                'round': r,
                'left_sbox_input': hex_state(trace_left[r]['input']),
                'right_sbox_input': hex_state(trace_right[r]['input']),
                'left_sbox_output': hex_state(trace_left[r]['sbox_output']),
                'right_sbox_output': hex_state(trace_right[r]['sbox_output']),
                'input_difference': hex_state(xor_states(trace_left[r]['input'],trace_right[r]['input'])),
                'sbox_output_difference': hex_state(xor_states(trace_left[r]['sbox_output'],trace_right[r]['sbox_output'])),
            } for r in range(4)
        ],
        'notes': ('This is an explicit witness in the chosen reduced-round core. '
                  'It does not prove a global maximum differential probability, '
                  'nor does it validate an official Midori-64 round-key schedule.'),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(evidence, indent=2), encoding='utf-8')
    print(json.dumps(evidence, indent=2))
    print('Saved witness:', args.out)


if __name__ == '__main__':
    main()
