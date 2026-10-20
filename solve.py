"""Solve and independently validate a generated Midori-64 MILP instance.

Requires a working Gurobi installation/license and gurobipy.
Does not assume previous objective bound 16 is optimal.
"""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
import json
from pathlib import Path

from build_models import difference_distribution_table, SHUFFLE
from midori_reference import validate_solver_trail


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--objective',choices=['active','probability'],default='active')
    parser.add_argument('--compact', action='store_true', help='use equality-eliminated LP (1857 rows, 1088 vars)')
    parser.add_argument('--time-limit',type=float,default=7200)
    parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent/'generated')
    args=parser.parse_args()
    try:
        import gurobipy as gp
    except ModuleNotFoundError:
        raise SystemExit('Install gurobipy in your Python environment to solve LPs.')
    directory=args.output_dir.resolve()
    suffix = '_compact' if args.compact else ''
    lp=directory/f'midori_4round_{args.objective}{suffix}.lp'
    if not lp.exists():
        raise SystemExit(f'Model not built: {lp}. First run: python build_models.py')
    m=gp.read(str(lp))
    m.setParam('TimeLimit',args.time_limit)
    m.setParam('MIPGap',0.0)
    m.setParam('MIPGapAbs',0.0)
    m.setParam('LogFile',str(directory/f'{args.objective}{suffix}.log'))
    m.optimize()
    result={
        'time_utc':datetime.now(timezone.utc).isoformat(),
        'objective':args.objective,
        'gurobi_version':'.'.join(map(str,gp.gurobi.version())),
        'status_code':m.Status,
        'status':{gp.GRB.OPTIMAL:'OPTIMAL', gp.GRB.TIME_LIMIT:'TIME_LIMIT', gp.GRB.INFEASIBLE:'INFEASIBLE', gp.GRB.INTERRUPTED:'INTERRUPTED'}.get(m.Status, str(m.Status)),
        'optimization_runtime_sec':m.Runtime,
        'solution_count':m.SolCount,
        'constraint_count':m.NumConstrs,
        'variable_count':m.NumVars,
        'warning':'OPTIMAL is established only when Gurobi returns OPTIMAL, subject to solver numerical tolerances.'
    }
    if m.SolCount:
        result['incumbent']=m.ObjVal
        result['best_bound']=m.ObjBound
        result['mip_gap']=m.MIPGap
        m.write(str(directory/f'{args.objective}{suffix}.sol'))
        assignment={v.VarName:int(round(v.X)) for v in m.getVars()}
        assert all(abs(v.X-round(v.X))<1e-5 for v in m.getVars())
        nibble=lambda prefix,r,c: sum(assignment[f'{prefix}_{r}_{c}_{b}']<<b for b in range(4))
        ddt=difference_distribution_table()
        trail=[]
        total_active=0
        total_cost=0
        for r in range(4):
            inp=[nibble('x',r,c) for c in range(16)]
            out=[nibble('y',r,c) for c in range(16)]
            for c,(dx,dy) in enumerate(zip(inp,out)):
                assert ddt[dx][dy] in (2,4,16),('Impossible DDT transition',r,c,dx,dy)
                active=int(dx!=0)
                cost={16:0,4:2,2:3}[ddt[dx][dy]]
                assert assignment[f'A_{r}_{c}']==active
                assert (assignment[f'P0_{r}_{c}'],assignment[f'P1_{r}_{c}'])=={16:(0,0),4:(0,1),2:(1,1)}[ddt[dx][dy]]
                total_active+=active;total_cost+=cost
            trail.append({'round':r,'input_differences':''.join(f'{d:X}' for d in inp),
                          'sbox_output_differences':''.join(f'{d:X}' for d in out),
                          'active_sboxes':sum(1 for d in inp if d)})
            if r<3:
                shuffled=[out[SHUFFLE[c]] for c in range(16)]
                if not args.compact:
                    assert all(nibble('z',r,c)==shuffled[c] for c in range(16))
                expected=[0]*16
                for col in range(4):
                    cells=[4*col+j for j in range(4)]
                    for row in range(4):
                        val=0
                        for j in range(4):
                            if j!=row:val ^= shuffled[cells[j]]
                        expected[cells[row]]=val
                assert all(nibble('x',r+1,c)==expected[c] for c in range(16))
        result['active_sboxes']=total_active
        result['negative_log2_differential_probability']=total_cost
        result['per_round_trail']=trail
        assert round(m.ObjVal)==(total_active if args.objective=='active' else total_cost)
        result['independently_verified_differential_transitions']=True
        result['independent_core_check']=validate_solver_trail(result)['independent_reference_validation']
    else:
        result['best_bound']=m.ObjBound if m.Status not in (gp.GRB.INFEASIBLE,gp.GRB.INF_OR_UNBD) else None
    destination=directory/f'{args.objective}{suffix}_result.json'
    destination.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='per_round_trail'},indent=2))
    print('Report:',destination)

if __name__=='__main__': main()
