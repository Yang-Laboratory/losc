import unittest
from psi4_losc.psi4_losc import scf
import psi4
import numpy as np

psi4.core.be_quiet()

class TestSCFFractional(unittest.TestCase):
    def setUp(self):
        self.mol_H = psi4.geometry("""
            0 2
            H
        symmetry c1
        """, name="H")
        self.mol_H_plus = psi4.geometry("""
            1 1
            H
        symmetry c1
        """, name="H+")
        self.mol_H2 = psi4.geometry("""
            0 1
            H
            H 1 1.2
        symmetry c1
        """, name="H2")
        self.mol_H2_plus = psi4.geometry("""
            1 2
            H
            H 1 1.2
        symmetry c1
        """, name="H2+")

        psi4.set_options({'guess':      'core',
                          'basis':      '3-21g',
                          'scf_type':   'df',
                          'e_convergence': 1e-8,
                          'maxiter': 50})

    def test_EN_curve_H_0_1(self):
        print("\n==> Test the E(N) curve for H atom in [0, 1] based on HF:")
        print("      guess=core, basis=3-21g, scf_type=DF")
        psi4.core.set_active_molecule(self.mol_H)
        psi4.set_options({'d_convergence': 1,
                          'reference':  'uhf'})
        name = 'hf'
        E_0 = scf(name, occ={'alpha': {'homo': 0}}).energy()
        E_1 = scf(name, occ={'alpha': {'homo': 1}}).energy()
        print(f"E_0={E_0}")
        print(f"E_1={E_1}")

        def ref_E_n(n):
            return (E_1 - E_0) * n

        n = 0.99
        E_n = scf(name, occ={'alpha': {'homo': n}}).energy()
        print(f'E_n={E_n}')
        self.assertAlmostEqual(ref_E_n(n), E_n)


if __name__ == '__main__':
    unittest.main()
