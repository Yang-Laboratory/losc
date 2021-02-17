"""
Test code for post-SCF-LOSC calculation, `psi4_losc.post_scf_losc()`.
"""

import unittest
import psi4
import numpy as np
import psi4_losc
from psi4_losc import post_scf_losc

psi4.core.be_quiet()

class TestPostSCFLOSCIntegerAufbau(unittest.TestCase):
    """
    Test `psi4_losc.post_scf_losc()` for integer systems with aufbau occupations.
    """
    def setUp(self):
        self.mol_H2_1A = psi4.geometry("""
            0 1
            H
            H 1 1
        symmetry c1
        """, name="H2: 1A")
        self.mol_H2_10A = psi4.geometry("""
            0 1
            H
            H 1 10
        symmetry c1
        """, name="H2: 10A")
        self.mol_H2_plus_1A = psi4.geometry("""
            1 2
            H
            H 1 1
        symmetry c1
        """, name="H2+: 1A")
        self.mol_H2_plus_10A = psi4.geometry("""
            1 2
            H
            H 1 10
        symmetry c1
        """, name="H2+: 10A")
        self.mol_H2O = psi4.geometry("""
            0 1
            O 0 0 0
            H 1 0 0
            H 0 1 0
        symmetry c1
        """, name="H2O")

        psi4.set_options({'guess':      'core',
                          'basis':      '6-31g',
                          'scf_type':   'df',
                          'e_convergence': 1e-8,
                          'maxiter': 100})

    def run_mol_no_correction(self, mol, precision=7):
        """
        Test cases in which LOSC gives no corrections. Compare the total energy
        between `psi4_losc.post_scf_losc()` and `psi4.energy()`. Three
        representative DFAs, svwn, blyp and b3lyp, are used for test.
        """
        psi4.core.set_active_molecule(mol)
        dfa_info = {'svwn': psi4_losc.SVWN,
                    'blyp': psi4_losc.BLYP,
                    'b3lyp': psi4_losc.B3LYP}
        for dfa in 'svwn blyp b3lyp'.split():
            print(f'    Test mol {mol.name()}:  dfa={dfa}')
            E_ref, dfa_wfn = psi4.energy(dfa, return_wfn=True)
            E_calc, _ = post_scf_losc(dfa_info[dfa], dfa_wfn)
            self.assertAlmostEqual(E_ref, E_calc, places=precision)

    def run_mol_correction(self, mol, E_ref, precision=4):
        """
        Test cases in which LOSC gives corrections. Compare the total energy
        between `psi4_losc.post_scf_losc()` and an estimated reference results.
        The estimated references are taken from qm4d calculations. Only use
        b3lyp functional for test.
        """
        psi4.core.set_active_molecule(mol)
        _, dfa_wfn = psi4.energy('b3lyp', return_wfn=True)
        E_calc, _ = post_scf_losc(psi4_losc.B3LYP, dfa_wfn)
        print(f'    Test mol {mol.name()}:  dfa=b3lyp, E_losc={E_calc}')
        self.assertAlmostEqual(E_ref, E_calc, places=precision)

    def test_uks_open_shell(self):
        """
        Test unrestricted calculations for open shell cases.
        """
        print("\n==> Test UKS Open Shell")
        optstash = psi4.driver.p4util.OptionsState(['SCF', 'REFERENCE'])
        psi4.core.set_local_option('SCF', 'REFERENCE', 'UHF')

        # LO = CO
        self.run_mol_no_correction(self.mol_H2_plus_1A)

        # LO != CO
        # Take from psi4_losc.post_scf_losc() self. Compared to qm4d output and
        # it is close to the 4-th digit.
        E_ref = -0.5008643357099776
        self.run_mol_correction(self.mol_H2_plus_10A, E_ref)

        optstash.restore()

    def test_uks_close_shell(self):
        """
        Test unrestricted calculations for close shell cases.
        """
        print("\n==> Test UKS Close Shell")
        optstash = psi4.driver.p4util.OptionsState(['SCF', 'REFERENCE'])
        psi4.core.set_local_option('SCF', 'REFERENCE', 'UHF')

        # LO = CO
        self.run_mol_no_correction(self.mol_H2_1A)

        # LO != CO
        E_ref = -0.7580683355488667
        self.run_mol_correction(self.mol_H2_10A, E_ref)

        optstash.restore()

    def test_rks_close_shell(self):
        """
        Test restricted calculations for close shell cases.
        """
        print("\n==> Test RKS Close Shell")
        optstash = psi4.driver.p4util.OptionsState(['SCF', 'REFERENCE'])
        psi4.core.set_local_option('SCF', 'REFERENCE', 'RHF')

        # LO = CO
        self.run_mol_no_correction(self.mol_H2_1A)

        # LO != CO
        E_ref = -0.7580683355488667
        self.run_mol_correction(self.mol_H2_10A, E_ref)

        optstash.restore()


if __name__ == '__main__':
    unittest.main()
