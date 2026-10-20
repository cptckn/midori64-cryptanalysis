"""Generate an equivalent, size-limited-license-friendly LP by eliminating ShuffleCell intermediates.

Original: 2,049 rows, 1,280 binary variables. Reduced: 1,857 rows,
1,088 binary variables. The equality z[r,out,bit] = y[r,SHUFFLE[out],bit]
is substituted throughout, so every feasible compact assignment uniquely
extends to a feasible original assignment with the same objective.

Use ``python compact_models.py`` after ``python build_models.py``.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import re
from build_models import SHUFFLE

_Z = re.compile(r"\bz_(\d+)_(\d+)_(\d+)\b")

def eliminate_shuffle(original: str) -> str:
    """Eliminate precisely the 192 z=y equality rows and their 192 z variables."""
    lines = original.splitlines()
    assert lines.count('Subject To') == 1 and lines.count('Binary') == 1
    section = None
    out = []
    removed_rows = removed_variables = 0
    for line in lines:
        if line == 'Subject To':
            section = 'rows'
        if line == 'Binary':
            section = 'binary'
        if section == 'rows' and re.match(r'^\s+sh_r\d+_c\d+_b\d+:', line):
            removed_rows += 1
            continue
        if section == 'binary' and re.fullmatch(r'\s+z_\d+_\d+_\d+', line):
            removed_variables += 1
            continue
        if section == 'rows':
            line = _Z.sub(lambda m: f'y_{m.group(1)}_{SHUFFLE[int(m.group(2))]}_{m.group(3)}', line)
        out.append(line)
    assert removed_rows == 192 and removed_variables == 192, (removed_rows, removed_variables)
    compact = '\n'.join(out) + '\n'
    assert not _Z.search(compact)
    assert sum(1 for line in out if re.match(r'^\s+[a-z]+_.*:', line)) == 1857
    assert len([line for line in compact.split('\nBinary\n')[1].split('\nEnd')[0].splitlines() if line.strip()]) == 1088
    return compact

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--output-dir', type=Path, default=Path(__file__).resolve().parent/'generated')
    args=p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for objective in ('active','probability'):
        inp=args.output_dir/f'midori_4round_{objective}.lp'
        if not inp.exists(): p.error(f'Missing {inp}; run python build_models.py first')
        dst=args.output_dir/f'midori_4round_{objective}_compact.lp'
        dst.write_text(eliminate_shuffle(inp.read_text(encoding='utf-8')), encoding='utf-8')
        print(f'{dst}: 1857 constraints, 1088 binary variables (ShuffleCell equality elimination)')

if __name__ == '__main__': main()
