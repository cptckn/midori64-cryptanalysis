"""Validate solver trails using an independent four-round Midori core.

Use: python verify_trail.py --report results/probability_trail.json --samples 0
(Does not require Gurobi.)

Crucial limitation: end-to-end Monte Carlo will normally find zero hits for a
very-low-probability full differential characteristic. This neither falsifies nor
validates the full-trail probability. The verified properties are the exact DDT
counts, concrete S-box witnesses, inter-round algebra, and operation of a
separate reduced-round implementation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import random
from midori_reference import (round_trace_ref, xor_states, hex_state, parse_hex_state,
                              mixcolumn_ref, shufflecell_ref, validate_solver_trail,
                              encrypt_four_round_core)


def monte_carlo_check(report, samples, seed):
    """Check the reference encryptor's difference traces on random input pairs.

    A hit is recorded only when ALL four S-box-layer input and output
    differences equal the specific MILP trail. Failure to observe any hit is
    expected if the path probability is tiny.
    """
    if samples == 0:
        return {'samples': 0, 'skipped': True,
                'note': 'Random-pair Monte Carlo skipped by request. The independently validated DDT, inter-round propagation, and archived standard-key witness are checked separately.'}
    rng = random.Random(seed)
    wanted = [(item['input_differences'].upper(), item['sbox_output_differences'].upper())
              for item in report['per_round_trail']]
    dx0 = tuple(int(c, 16) for c in wanted[0][0])
    hits = 0
    consistency_checks = 0
    first_pair = None
    for _ in range(samples):
        left = tuple(rng.randrange(16) for _ in range(16))
        right = xor_states(left, dx0)
        left_trace = round_trace_ref(left)
        right_trace = round_trace_ref(right)
        assert encrypt_four_round_core(left) == left_trace[-1]['sbox_output']
        assert encrypt_four_round_core(right) == right_trace[-1]['sbox_output']
        observed = [(hex_state(xor_states(la['input'], ra['input'])),
                     hex_state(xor_states(la['sbox_output'], ra['sbox_output'])))
                    for la, ra in zip(left_trace, right_trace)]
        assert observed[0][0] == wanted[0][0]
        # Check reference cipher satisfies per-round bitwise difference propagation.
        for i in range(1, 4):
            assert parse_hex_state(observed[i][0]) == mixcolumn_ref(
                shufflecell_ref(parse_hex_state(observed[i-1][1])))
        consistency_checks += 1
        if observed == wanted:
            hits += 1
            if first_pair is None:
                first_pair = {'plaintext_left': hex_state(left),
                              'plaintext_right': hex_state(right)}
    return {'samples': samples, 'seed': seed, 'all_sampled_core_traces_consistent': consistency_checks == samples,
            'full_4round_path_hits': hits, 'example_full_path_pair': first_pair,
            'note': ('No hits at a rare probability do NOT disprove the trail. '
                     'The zero-round-key core is an explicitly chosen reduced-round variant, '
                     'not a claimed reproduction of all Midori-64 round constants.')}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--objective', choices=('active', 'probability'), default='probability')
    p.add_argument('--report', type=Path, help='Override generated/{objective}_result.json')
    p.add_argument('--samples', type=int, default=1000)
    p.add_argument('--seed', type=int, default=2026)
    args = p.parse_args()
    if args.samples < 0:
        p.error('--samples must be nonnegative')
    folder = Path(__file__).resolve().parent / 'generated'
    folder.mkdir(parents=True, exist_ok=True)
    report_file = args.report or (folder / f'{args.objective}_result.json')
    if not report_file.exists():
        raise SystemExit(f'No solver report found: {report_file}. Run py solve.py --objective {args.objective} first.')
    report = json.loads(report_file.read_text(encoding='utf-8'))
    if report.get('solution_count', 0) < 1:
        raise SystemExit('Gurobi has not reported a feasible trail in the given JSON file.')
    result = validate_solver_trail(report)
    result['reduced_round_software_test'] = monte_carlo_check(report, args.samples, args.seed)
    output = folder / f'{args.objective}_verification.json'
    output.write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2))
    print('Verification report:', output)

if __name__ == '__main__':
    main()
