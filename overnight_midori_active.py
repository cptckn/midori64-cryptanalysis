"""Overnight attempt to certify four-round Midori-64 active-S-box minimum.

Uses the unchanged generated/midori_4round_active.lp model and a solution
from generated/probability.sol as an optional MIP start.  Does not modify any
other reports, LP files, or solutions.

Run on Windows (from this repository root):
    py -u overnight_midori_active.py --hours 8

Gurobi 13 with an academic license is required.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import tempfile


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hours', type=float, default=8.0, help='time limit in hours (default: 8)')
    parser.add_argument('--threads', type=int, default=8, help='Gurobi threads (default: 8)')
    parser.add_argument('--nodefile-gb', type=float, default=0.5,
                        help='begin writing search nodes to disk once node storage exceeds this size in GB')
    args = parser.parse_args()
    if args.hours <= 0 or args.threads < 1 or args.nodefile_gb <= 0:
        parser.error('hours, threads and nodefile-gb must be positive')

    try:
        import gurobipy as gp
    except ImportError:
        parser.error('gurobipy is not installed in the Python environment')

    root = Path(__file__).resolve().parent
    generated = root / 'generated'
    generated.mkdir(exist_ok=True)
    lp = generated / 'midori_4round_active.lp'
    if not lp.exists():
        parser.error(f'Missing {lp}. Run: py build_models.py')
    # Store potentially large Gurobi node files in the OS temporary directory.
    nodestore = Path(tempfile.gettempdir()) / 'midori_gurobi_overnight_nodes'
    nodestore.mkdir(parents=True, exist_ok=True)
    report_path = generated / 'overnight_active_result.json'
    progress_path = generated / 'overnight_active_progress.json'
    solution_path = generated / 'overnight_active.sol'

    model = gp.read(str(lp))
    model.setParam('TimeLimit', float(args.hours * 3600))
    model.setParam('MIPFocus', 3)  # focus on tightening the dual bound
    model.setParam('MIPGap', 0)
    model.setParam('MIPGapAbs', 0)
    model.setParam('Threads', args.threads)
    model.setParam('NodefileStart', args.nodefile_gb)
    model.setParam('NodefileDir', str(nodestore))
    model.setParam('LogFile', str(generated / 'overnight_active.log'))

    mip_start = generated / 'probability.sol'
    if not mip_start.exists() and (root/'results'/'active_solution.sol').exists():
        mip_start = root/'results'/'active_solution.sol'
    if mip_start.exists():
        try:
            model.read(str(mip_start))
            print(f'Loaded previously found feasible trail as a MIP start: {mip_start}', flush=True)
        except gp.GurobiError as exc:
            print(f'Could not load optional MIP start ({exc}); solving without it.', flush=True)
    else:
        print('No probability.sol file found. Gurobi will find its own incumbent.', flush=True)

    last_checkpoint = [-600.0]
    def on_progress(m, where):
        if where != gp.GRB.Callback.MIP:
            return
        runtime = m.cbGet(gp.GRB.Callback.RUNTIME)
        if runtime < last_checkpoint[0] + 600.0:
            return
        last_checkpoint[0] = runtime
        primal = m.cbGet(gp.GRB.Callback.MIP_OBJBST)
        dual = m.cbGet(gp.GRB.Callback.MIP_OBJBND)
        snapshot = {
            'time_utc': datetime.now(timezone.utc).isoformat(),
            'runtime_sec': round(runtime, 1),
            'best_objective': None if primal >= gp.GRB.INFINITY / 2 else float(primal),
            'best_bound': None if abs(dual) >= gp.GRB.INFINITY / 2 else float(dual),
            'nodes_explored': int(m.cbGet(gp.GRB.Callback.MIP_NODCNT)),
            'nodes_remaining': int(m.cbGet(gp.GRB.Callback.MIP_NODLFT)),
            'note': 'Intermediate snapshot only; not a proof of optimality.',
        }
        progress_path.write_text(json.dumps(snapshot, indent=2), encoding='utf-8')
        print(f"[Checkpoint] elapsed={runtime/3600:.2f}h incumbent={snapshot['best_objective']} bound={snapshot['best_bound']}", flush=True)

    print('Starting active-S-box minimization; this can run for hours.', flush=True)
    print(f'Results: {report_path}', flush=True)
    print(f'Search-node scratch directory: {nodestore}', flush=True)
    model.optimize(on_progress)

    statuses = {
        gp.GRB.OPTIMAL: 'OPTIMAL',
        gp.GRB.TIME_LIMIT: 'TIME_LIMIT',
        gp.GRB.INTERRUPTED: 'INTERRUPTED',
        gp.GRB.INFEASIBLE: 'INFEASIBLE',
        gp.GRB.INF_OR_UNBD: 'INF_OR_UNBD',
        gp.GRB.NODE_LIMIT: 'NODE_LIMIT',
        gp.GRB.MEM_LIMIT: 'MEM_LIMIT',
    }
    result = {
        'time_utc': datetime.now(timezone.utc).isoformat(),
        'objective': 'minimize_number_of_active_sboxes',
        'status': statuses.get(model.Status, str(model.Status)),
        'status_code': model.Status,
        'runtime_sec': model.Runtime,
        'solution_count': model.SolCount,
        'constraint_count': model.NumConstrs,
        'variable_count': model.NumVars,
        'setting_mipfocus': 3,
        'setting_threads': args.threads,
        'setting_nodefile_start_gb': args.nodefile_gb,
        'note': 'Only OPTIMAL certifies global optimality within Gurobi numerical tolerances.',
    }
    if model.Status not in (gp.GRB.INFEASIBLE, gp.GRB.INF_OR_UNBD):
        result['best_bound'] = float(model.ObjBound)
    if model.SolCount:
        result['incumbent'] = float(model.ObjVal)
        result['mip_gap'] = float(model.MIPGap)
        model.write(str(solution_path))
        result['solution_path'] = str(solution_path)
    report_path.write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2), flush=True)
    print('Report:', report_path, flush=True)

if __name__ == '__main__':
    main()
