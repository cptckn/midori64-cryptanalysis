"""Independent, small Midori-64 four-round CORE implementation for verification.

This module intentionally does not import build_models.py. Constants and state
ordering have been cross-checked against Banik et al., 'Midori: A Block Cipher
for Low Energy' (2015), Sections 3.1--3.3. It implements the four-SubCell
truncated Midori core: three (SubCell, ShuffleCell, MixColumn, KeyAdd)
transitions and a final SubCell. It is NOT a complete 16-round cipher and is
not intended to be used for encryption or real-world security.

State is a tuple of 16 nibbles with cells ordered column-major; each column
consists of consecutive indices (0..3, 4..7, 8..11, 12..15).
"""
from __future__ import annotations

SBOX_REF = (0xC, 0xA, 0xD, 0x3, 0xE, 0xB, 0xF, 0x7,
            0x8, 0x9, 0x1, 0x5, 0x0, 0x2, 0x4, 0x6)
SHUFFLE_REF = (0, 10, 5, 15, 14, 4, 11, 1,
               9, 3, 12, 6, 7, 13, 2, 8)
ZERO = (0,) * 16


def validate_state(state):
    v = tuple(state)
    if len(v) != 16 or any(not isinstance(x, int) or not 0 <= x < 16 for x in v):
        raise ValueError('A state must contain exactly 16 nibbles (integers 0..15)')
    return v


def parse_hex_state(s):
    if not isinstance(s, str) or len(s) != 16:
        raise ValueError('Expected 16 hex digits encoding 16 Midori-64 nibbles')
    try:
        return validate_state(tuple(int(ch, 16) for ch in s))
    except ValueError as exc:
        raise ValueError('Invalid 16-nibble hex state') from exc


def hex_state(state):
    return ''.join(f'{n:X}' for n in validate_state(state))


def xor_states(left, right):
    left, right = validate_state(left), validate_state(right)
    return tuple(a ^ b for a, b in zip(left, right))


def subcell_ref(state):
    return tuple(SBOX_REF[n] for n in validate_state(state))


def shufflecell_ref(state):
    state = validate_state(state)
    return tuple(state[j] for j in SHUFFLE_REF)


def mixcolumn_ref(state):
    state = validate_state(state)
    output = []
    for c in (0, 4, 8, 12):
        a, b, d, e = state[c:c+4]
        output.extend((b ^ d ^ e, a ^ d ^ e, a ^ b ^ e, a ^ b ^ d))
    return tuple(output)


def round_trace_ref(state, round_keys=(ZERO, ZERO, ZERO)):
    """Return state before each SubCell and after each SubCell for 4 rounds.

    For the differential model, identical round keys are used for both texts.
    Midori's real round constants are key-dependent affine translations: they
    cancel in the immediate difference but may influence multi-round empirical
    hit rates. The caller may supply any three 16-nibble round keys.
    """
    if len(round_keys) != 3:
        raise ValueError('Exactly three intermediate 16-nibble round keys are required')
    curr = validate_state(state)
    trace = []
    for r in range(4):
        out = subcell_ref(curr)
        trace.append({'round': r, 'input': curr, 'sbox_output': out})
        if r < 3:
            curr = xor_states(mixcolumn_ref(shufflecell_ref(out)),
                              validate_state(round_keys[r]))
    return trace


def encrypt_four_round_core(state, round_keys=(ZERO, ZERO, ZERO)):
    return round_trace_ref(state, round_keys)[-1]['sbox_output']


def difference_distribution_table_ref():
    ddt = [[0] * 16 for _ in range(16)]
    for alpha in range(16):
        for x in range(16):
            beta = SBOX_REF[x] ^ SBOX_REF[x ^ alpha]
            ddt[alpha][beta] += 1
    return ddt


def nibble_pair_for_transition(dx, dy):
    """Return concrete (x,x^dx) with output difference dy, or None."""
    for x in range(16):
        if (SBOX_REF[x] ^ SBOX_REF[x ^ dx]) == dy:
            return x, x ^ dx
    return None


def validate_solver_trail(report):
    """Check every nibble transition and inter-round diffusion independently.

    Does not prove existence of a plaintext pair following the ENTIRE trail:
    local S-box DDT feasibility need not imply full multi-round feasibility.
    """
    t = report.get('per_round_trail')
    if not isinstance(t, list) or len(t) != 4:
        raise ValueError('Solver report must contain four rounds in per_round_trail')
    ddt = difference_distribution_table_ref()
    sum_active = sum_cost = 0
    earlier_out = None
    verification = []
    for r, record in enumerate(t):
        if record.get('round') != r:
            raise ValueError(f'Incorrect round index at entry {r}')
        inp = parse_hex_state(record['input_differences'])
        out = parse_hex_state(record['sbox_output_differences'])
        active = 0
        cost = 0
        local_witnesses = []
        for c, (dx, dy) in enumerate(zip(inp, out)):
            count = ddt[dx][dy]
            if count not in (2, 4, 16):
                raise ValueError(f'Round {r}, cell {c}: forbidden S-box differential {dx:X}->{dy:X}')
            pair = nibble_pair_for_transition(dx, dy)
            if pair is None:
                raise ValueError(f'Round {r}, cell {c}: no plaintext-pair witness')
            local_witnesses.append(pair)
            active += int(dx != 0)
            cost += {2: 3, 4: 2, 16: 0}[count]
        if active != record.get('active_sboxes'):
            raise ValueError(f'Round {r}: active S-box count does not match report')
        if earlier_out is not None:
            expected = mixcolumn_ref(shufflecell_ref(earlier_out))
            if inp != expected:
                raise ValueError(f'Round {r}: input difference fails ShuffleCell/MixColumn')
        sum_active += active
        sum_cost += cost
        verification.append({'round': r, 'active_sboxes': active,
                             'minus_log2_round_probability': cost,
                             'nibble_transitions_with_local_pair_witnesses': len(local_witnesses)})
        earlier_out = out
    if report.get('active_sboxes') != sum_active:
        raise ValueError('Global active S-box count does not match independently recomputed result')
    if report.get('negative_log2_differential_probability') != sum_cost:
        raise ValueError('Global probability cost does not match independently recomputed result')
    obj = report.get('objective')
    if obj in ('active', 'probability') and 'incumbent' in report:
        expected = sum_active if obj == 'active' else sum_cost
        if abs(report['incumbent'] - expected) > 1e-5:
            raise ValueError('Gurobi incumbent differs from independent objective')
    return {'independent_reference_validation': 'PASS',
            'end_to_end_plaintext_pair_existence': None,
            'end_to_end_explanation': 'Local DDT validation alone does not construct a full four-round plaintext pair.',
            'active_sboxes': sum_active, 'negative_log2_probability': sum_cost,
            'approximate_trail_probability': f'2^-{sum_cost}',
            'rounds': verification}
