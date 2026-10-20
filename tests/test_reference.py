"""Independent reference implementation, trail-checker, and failure tests."""
import sys
from pathlib import Path
import random
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import build_models as b
import midori_reference as ref


def trail_from_witness(witness, objective):
    inp_state = lambda name, r: tuple(
        sum(witness[f'{name}_{r}_{c}_{bit}'] << bit for bit in range(4)) for c in range(16))
    ddt = b.difference_distribution_table()
    details = []
    for r in range(4):
        xi, yi = inp_state('x', r), inp_state('y', r)
        details.append({'round': r, 'input_differences': ref.hex_state(xi),
                        'sbox_output_differences': ref.hex_state(yi),
                        'active_sboxes': sum(bool(x) for x in xi)})
    active = sum(row['active_sboxes'] for row in details)
    cost = sum({16: 0, 4: 2, 2: 3}[ddt[x][y]] for row in details
               for x, y in zip(ref.parse_hex_state(row['input_differences']),
                               ref.parse_hex_state(row['sbox_output_differences'])))
    return {'objective': objective, 'solution_count': 1,
            'incumbent': active if objective == 'active' else cost,
            'active_sboxes': active,
            'negative_log2_differential_probability': cost,
            'per_round_trail': details}


class TestIndependentReference(unittest.TestCase):
    def test_constants_and_ddt_against_builder(self):
        self.assertEqual(ref.SBOX_REF, b.SBOX)
        self.assertEqual(ref.SHUFFLE_REF, b.SHUFFLE)
        self.assertEqual(ref.difference_distribution_table_ref(),
                         b.difference_distribution_table())
        for row in ref.difference_distribution_table_ref():
            self.assertEqual(sum(row), 16)

    def test_shufflecell_and_mixcolumn_known_index_and_basis(self):
        inp = tuple(range(16))
        self.assertEqual(ref.shufflecell_ref(inp), ref.SHUFFLE_REF)
        for c in range(4):
            for row in range(4):
                state = [0] * 16
                state[4*c+row] = 0xB
                out = ref.mixcolumn_ref(state)
                self.assertEqual(tuple(out[4*c:4*c+4]),
                                 tuple(0 if i == row else 0xB for i in range(4)))

    def test_each_local_transition_has_concrete_nibble_pair_witness(self):
        ddt = ref.difference_distribution_table_ref()
        for dx in range(16):
            for dy in range(16):
                p = ref.nibble_pair_for_transition(dx, dy)
                self.assertEqual(p is not None, ddt[dx][dy] > 0)
                if p:
                    x, y = p
                    self.assertEqual(x ^ y, dx)
                    self.assertEqual(ref.SBOX_REF[x] ^ ref.SBOX_REF[y], dy)

    def test_all_constructed_trails_independently_validate(self):
        from test_models import random_witness
        for seed in range(30):
            witness = random_witness(seed)
            for objective in ('active', 'probability'):
                report = trail_from_witness(witness, objective)
                checked = ref.validate_solver_trail(report)
                self.assertEqual(checked['independent_reference_validation'], 'PASS')
                self.assertEqual(checked['active_sboxes'], report['active_sboxes'])

    def test_tampered_trail_is_rejected(self):
        from test_models import random_witness
        report = trail_from_witness(random_witness(42), 'active')
        report['per_round_trail'][2]['input_differences'] = 'F' * 16
        with self.assertRaises(ValueError):
            ref.validate_solver_trail(report)
        report = trail_from_witness(random_witness(42), 'active')
        report['per_round_trail'][0]['sbox_output_differences'] = 'G' * 16
        with self.assertRaises(ValueError):
            ref.validate_solver_trail(report)

    def test_real_reduced_round_encryption_pair_trace(self):
        rng = random.Random(112)
        for _ in range(100):
            a = tuple(rng.randrange(16) for _ in range(16))
            b_ = tuple(rng.randrange(16) for _ in range(16))
            a_trace, b_trace = ref.round_trace_ref(a), ref.round_trace_ref(b_)
            self.assertEqual(ref.encrypt_four_round_core(a), a_trace[-1]['sbox_output'])
            for i in range(3):
                diff = ref.xor_states(a_trace[i]['sbox_output'], b_trace[i]['sbox_output'])
                next_expected = ref.mixcolumn_ref(ref.shufflecell_ref(diff))
                self.assertEqual(ref.xor_states(a_trace[i+1]['input'], b_trace[i+1]['input']),
                                 next_expected)


if __name__ == '__main__':
    unittest.main(verbosity=2)
