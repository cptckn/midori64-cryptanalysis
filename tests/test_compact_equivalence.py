"""Check compact LP elimination algebraically and with the archived 16-active trail."""
import sys, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from build_models import model_text, read_hrep, SHUFFLE
from compact_models import eliminate_shuffle
from audit_solution import parse_solution,parse_rows

class CompactEquivalenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original,_counts,_vars=model_text('active',read_hrep(ROOT/'selected_inequalities.txt'))
        cls.compact=eliminate_shuffle(cls.original)
        cls.assignment,_=parse_solution(ROOT/'results'/'active_solution.sol')

    def test_compact_dimensions_and_feasibility(self):
        rows=list(parse_rows(self.compact))
        self.assertEqual(len(rows),1857)
        binary=self.compact.split('\nBinary\n')[1].split('\nEnd')[0].splitlines()
        variables={x.strip() for x in binary if x.strip()}
        self.assertEqual(len(variables),1088)
        projected={key:val for key,val in self.assignment.items() if key in variables}
        self.assertEqual(len(projected),1088)
        for name, terms, op, rhs in rows:
            lhs=sum(coef*projected[var] for var,coef in terms)
            self.assertTrue(lhs>=rhs if op=='>=' else lhs==rhs, name)
        self.assertEqual(sum(projected[f'A_{r}_{c}'] for r in range(4) for c in range(16)),16)

    def test_all_shuffle_equalities_stipulate_unique_extension(self):
        for r in range(3):
            for c in range(16):
                for b in range(4):
                    self.assertEqual(self.assignment[f'z_{r}_{c}_{b}'],self.assignment[f'y_{r}_{SHUFFLE[c]}_{b}'])

if __name__=='__main__':unittest.main()
