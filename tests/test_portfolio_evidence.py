"""Recorded results: check against independently implemented Midori operations."""
import sys
import json
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from audit_solution import audit_solution
from midori_reference import (parse_hex_state, hex_state, xor_states,
                              round_trace_ref, validate_solver_trail)
from find_midori_standard_key_pair import split_master_key, derive_keys

class PortfolioEvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.trail=json.loads((ROOT/'results'/'probability_trail.json').read_text())
        cls.witness=json.loads((ROOT/'results'/'standard_key_witness.json').read_text())
    def test_recorded_active_solution_satisfies_every_lp_row(self):
        self.assertEqual(audit_solution(ROOT/'results'/'active_solution.sol'),16)
    def test_probability_trail_is_ddt_and_diffusion_consistent(self):
        report=validate_solver_trail(self.trail)
        self.assertEqual(report['independent_reference_validation'],'PASS')
        self.assertEqual(report['negative_log2_probability'],32)
    def test_standard_key_witness_replays_full_trail(self):
        w=self.witness
        k0,k1=split_master_key(w['master_key']); wk,rks=derive_keys(k0,k1)
        self.assertEqual(hex_state(wk),w['whitening_key'])
        self.assertEqual([hex_state(x) for x in rks],w['round_keys'])
        p0=parse_hex_state(w['plaintext_left']);p1=parse_hex_state(w['plaintext_right'])
        t0=round_trace_ref(xor_states(p0,wk),rks)
        t1=round_trace_ref(xor_states(p1,wk),rks)
        for index,(left,right,expected) in enumerate(zip(t0,t1,self.trail['per_round_trail'])):
            self.assertEqual(hex_state(xor_states(left['input'],right['input'])),expected['input_differences'],index)
            self.assertEqual(hex_state(xor_states(left['sbox_output'],right['sbox_output'])),expected['sbox_output_differences'],index)
            self.assertEqual(hex_state(left['input']),w['per_round_inputs'][index])
        self.assertEqual(hex_state(xor_states(t0[-1]['sbox_output'],wk)),w['ciphertext_left'])
        self.assertEqual(hex_state(xor_states(t1[-1]['sbox_output'],wk)),w['ciphertext_right'])

if __name__=='__main__':unittest.main()
