"""Independent structural and constructive checks (no Gurobi required)."""
from pathlib import Path
import sys
import random
import re
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import build_models as b


def vector_to_dict(prefix, round_id, cell, value):
    return {f"{prefix}_{round_id}_{cell}_{bit}":(value>>bit)&1 for bit in range(4)}


def random_witness(seed):
    rng=random.Random(seed)
    ddt=b.difference_distribution_table()
    state=[0]*16
    for _ in range(rng.randint(1,4)):
        state[rng.randrange(16)]=rng.randrange(1,16)
    assignments={}
    for r in range(4):
        output=[0]*16
        for c,dx in enumerate(state):
            candidates=[dy for dy in range(16) if ddt[dx][dy]>0]
            dy=rng.choice(candidates)
            output[c]=dy
            assignments.update(vector_to_dict('x',r,c,dx))
            assignments.update(vector_to_dict('y',r,c,dy))
            assignments[b.a(r,c)] = int(dx!=0)
            count=ddt[dx][dy]
            p0,p1={16:(0,0),4:(0,1),2:(1,1)}[count]
            assignments[b.pp(0,r,c)]=p0
            assignments[b.pp(1,r,c)]=p1
        if r<3:
            shuffled=[output[b.SHUFFLE[c]] for c in range(16)]
            for c,value in enumerate(shuffled):assignments.update(vector_to_dict('z',r,c,value))
            next_state=[0]*16
            for col in range(4):
                cells=[4*col+i for i in range(4)]
                for row in range(4):
                    out=0
                    for j in range(4):
                        if j!=row:out ^= shuffled[cells[j]]
                    next_state[cells[row]]=out
                    for bit in range(4):
                        bits=[(shuffled[cells[j]]>>bit)&1 for j in range(4) if j!=row] + [(out>>bit)&1]
                        num=sum(bits)
                        assert num in (0,2,4)
                        d1,d2={0:(0,0),2:(1,1),4:(1,0)}[num]
                        assignments[b.d(r,col,row,bit,1)]=d1
                        assignments[b.d(r,col,row,bit,2)]=d2
            state=next_state
    return assignments


def lp_rows(lines):
    in_constraints=False
    for line in lines:
        line=line.strip()
        if line=='Subject To':in_constraints=True;continue
        if line=='Binary':break
        if not in_constraints or not line:continue
        label,statement=line.split(': ',1)
        match=re.fullmatch(r'(.*?)\s*(>=|=)\s*(-?\d+)',statement)
        assert match, (label,statement)
        expr,op,rhs=match.groups()
        tokens=re.findall(r'([+-]?)\s*(\d*)\s*([A-Za-z][\w]*)',expr)
        yield label,[(name,(-1 if sign=='-' else 1)*int(num or 1)) for sign,num,name in tokens],op,int(rhs)


class TestMidori(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.hrep=b.read_hrep(ROOT/'selected_inequalities.txt')
        b.verify_hrep(cls.hrep)
        b.verify_xor3()
        cls.models={mode:b.model_text(mode,cls.hrep)[0].splitlines() for mode in ('active','probability')}
    def test_models_have_same_constraints(self):
        for mode in ('active','probability'):
            rows=list(lp_rows(self.models[mode]))
            self.assertEqual(len(rows),2049)
        self.assertEqual(list(lp_rows(self.models['active'])),list(lp_rows(self.models['probability'])))
    def test_random_valid_trails_satisfy_every_constraint(self):
        rows=list(lp_rows(self.models['active']))
        for seed in range(20):
            assignment=random_witness(seed)
            for label,terms,op,rhs in rows:
                try: lhs=sum(coef*assignment[var] for var,coef in terms)
                except KeyError as exc:raise AssertionError(f'Missing variable in {label}: {exc}')
                self.assertTrue((lhs==rhs) if op=='=' else (lhs>=rhs),f'{seed} {label} {lhs} {op} {rhs}')
    def test_all_generated_variables_are_binary(self):
        names=set()
        for line in self.models['active'][self.models['active'].index('Binary')+1:]:
            if line=='End':break
            names.add(line.strip())
        self.assertEqual(len(names),1280)
        for seed in range(4):
            assignment=random_witness(seed)
            self.assertEqual(set(assignment),names)
            self.assertEqual(set(assignment.values()).difference({0,1}),set())

if __name__=='__main__':unittest.main(verbosity=2)
