"""
Test code for the extended SCF procedure, `psi4_losc.scf.scf()` for fractional
systems.
"""

import unittest
import psi4
import numpy as np
from psi4_losc.scf import scf

psi4.core.be_quiet()


class TestSelfSCFFractional(unittest.TestCase):
    def setUp(self):
        self.mol_H = psi4.geometry("""
            0 2
            H
        symmetry c1
        """, name="H")
        self.mol_H2_plus = psi4.geometry("""
            1 2
            H
            H 1 1.2
        symmetry c1
        """, name="H2+")

        psi4.set_options({'guess':      'core',
                          'basis':      '3-21g',
                          'scf_type':   'direct',
                          'e_convergence': 1e-8,
                          'maxiter': 50})

    def test_EN_curve_H_gs_HF(self):
        """
        Test the ground state of E(N) curve of one-electron system: H atom from HF.
        In this case, HF results are exact.
        """
        print(
            "\n==> Test the ground state E(N) curve of H atom in [0, 1]:")
        print("    method=HF, guess=core, basis=3-21g, scf_type=DF")
        psi4.core.set_active_molecule(self.mol_H)
        optstash = psi4.driver.p4util.OptionsState(
            ['reference'])
        psi4.set_options({'reference':  'uhf'})

        name = 'hf'
        E_0 = scf(name, occ={'alpha': {'homo': 0}}).energy()
        E_1 = scf(name, occ={'alpha': {'homo': 1}}).energy()

        def ref_E_n(n):
            return (E_1 - E_0) * n + E_0
        for n in np.arange(0, 1, 0.1):
            E_n = scf(name, occ={'alpha': {'homo': n}}).energy()
            self.assertAlmostEqual(ref_E_n(n), E_n)

        optstash.restore()

    def test_EN_curve_H_gs_b3lyp(self):
        """
        Test the ground state of E(N) curve of one-electron system: H atom from B3LYP.
        """
        print(
            "\n==> Test the ground state E(N) curve of H atom in [0, 1]:")
        print("    method=b3lyp, guess=core, basis=3-21g, scf_type=DF")
        psi4.core.set_active_molecule(self.mol_H)
        optstash = psi4.driver.p4util.OptionsState(
            ['reference'])
        psi4.set_options({'reference':  'uhf'})

        name = 'b3lyp'
        # ref comes from self generation. compared with qm4d calculations
        # and these results look reasonable.
        E_ref = {
            0: 0,
            0.1: -0.06230425262483324,
            0.2: -0.12438454229516446,
            0.3: -0.18386651447939936,
            0.4: -0.24007115938269238,
            0.5: -0.29266784724316897,
            0.6: -0.34147427067192015,
            0.7: -0.38638636474480936,
            0.8: -0.42734614613486516,
            0.9: -0.46432455648727333,
            1.0: -0.4973113890200121,
        }
        for n, ref in E_ref.items():
            E_n = scf(name, occ={'alpha': {'homo': n}}).energy()
            self.assertAlmostEqual(ref, E_n)

        optstash.restore()

    def test_EN_curve_H2_plus_gs_HF(self):
        """
        Test the ground state E(N) curve of one-electron system: H2+ atom from HF.
        This is a case with multiple atoms. In this case, HF results are exact.
        """
        print("\n==> Test the ground state of E(N) curve of H2+ in [0, 1]:")
        print("      method=HF, guess=core, basis=3-21g, scf_type=DF")
        psi4.core.set_active_molecule(self.mol_H2_plus)
        optstash = psi4.driver.p4util.OptionsState(
            ['reference'])
        psi4.set_options({'reference':  'uhf'})

        name = 'hf'
        psi4.core.set_output_file('t', False)
        E_0 = scf(name, occ={'alpha': {'homo': 0}}).energy()
        E_1 = scf(name, occ={'alpha': {'homo': 1}}).energy()

        def ref_E_n(n):
            return (E_1 - E_0) * n + E_0
        for n in np.arange(0, 1, 0.5):
            E_n = scf(name, occ={'alpha': {'homo': n}}).energy()
            self.assertAlmostEqual(ref_E_n(n), E_n)

        optstash.restore()


if __name__ == '__main__':
    unittest.main()
