"""Independently check a saved Gurobi .sol against the *regenerated* LP.

Only verifies the primal solution and objective; the branch-and-bound optimality
certificate is the separately preserved solver log and OPTIMAL status.
No Gurobi installation required for this audit.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import re
from build_models import model_text, read_hrep

ROOT = Path(__file__).resolve().parent

def parse_solution(path: Path) -> tuple[dict[str, int], int]:
    lines = path.read_text(encoding='utf-8').splitlines()
    match = next((re.search(r'Objective value\s*=\s*([0-9.]+)', x) for x in lines if 'Objective value' in x), None)
    if match is None or match.group(1) is None:
        raise ValueError('Missing objective header')
    output = {}
    for line in lines:
        if line.startswith('#') or not line.strip(): continue
        fields = line.split()
        if len(fields) != 2: raise ValueError('Invalid solution line: ' + line)
        name, raw_value = fields
        v = float(raw_value)
        if abs(v - round(v)) > 1e-6 or round(v) not in (0,1):
            raise ValueError(f'Nonbinary variable {name}={v}')
        if name in output: raise ValueError(f'Duplicate variable {name}')
        output[name] = round(v)
    return output, int(round(float(match.group(1))))

def parse_rows(lp: str):
    in_rows = False
    for raw in lp.splitlines():
        line=raw.strip()
        if line == 'Subject To': in_rows = True; continue
        if line == 'Binary': break
        if not in_rows or not line: continue
        label, statement = line.split(': ', 1)
        match=re.fullmatch(r'(.*?)\s*(>=|=)\s*(-?\d+)',statement)
        if not match: raise ValueError(f'Cannot parse LP row: {label}')
        expr, op, rhs = match.groups()
        terms=re.findall(r'([+-]?)\s*(\d*)\s*([A-Za-z][\w]*)',expr)
        yield label, [(name,(-1 if sign=='-' else 1)*int(num or 1)) for sign,num,name in terms], op, int(rhs)

def audit_solution(solution_path: Path):
    values, expected_objective = parse_solution(solution_path)
    lp, counts, n_variables = model_text('active',read_hrep(ROOT/'selected_inequalities.txt'))
    binary = lp.split('\nBinary\n',1)[1].split('\nEnd',1)[0].splitlines()
    names = {line.strip() for line in binary if line.strip()}
    if set(values) != names:
        raise AssertionError(f'Missing variables: {len(names-set(values))}; extra variables: {len(set(values)-names)}')
    violations=[]
    rows=list(parse_rows(lp))
    for label,terms,op,rhs in rows:
        lhs=sum(coef*values[name] for name,coef in terms)
        if not ((lhs>=rhs) if op=='>=' else (lhs==rhs)):
            violations.append((label,lhs,op,rhs))
    active_count=sum(values[f'A_{r}_{c}'] for r in range(4) for c in range(16))
    if active_count != expected_objective: violations.append(('objective',active_count,'=',expected_objective))
    print(f'Checked {len(rows)} LP rows and {len(values)} binary variables; objective={active_count}; violations={len(violations)}')
    if violations: raise AssertionError(violations[:5])
    assert len(rows)==2049 and len(values)==1280 and active_count==16
    return active_count

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--solution',type=Path,default=ROOT/'results'/'active_solution.sol')
    args=ap.parse_args()
    audit_solution(args.solution)

if __name__=='__main__':main()
