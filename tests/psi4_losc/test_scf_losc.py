"""
Test code for SCF-LOSC (frozen-LO) calculations.
"""

import unittest
import psi4_losc
import psi4
from psi4_losc import scf_losc
from psi4_losc.scf import scf_losc as my_scf_losc

psi4.core.be_quiet()


class TestSCFLOSCIntegerAufbau(unittest.TestCase):
    """
    Test `psi4_losc.scf_losc()` for integer systems with aufbau occupations.
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
            H 1 1
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
        between `psi4_losc.scf_losc()` and `psi4.energy()`. Three representative
        DFAs, svwn, blyp and b3lyp, are used for test.
        """
        psi4.core.set_active_molecule(mol)
        dfa_info = {'svwn': psi4_losc.SVWN,
                    'blyp': psi4_losc.BLYP,
                    'b3lyp': psi4_losc.B3LYP}
        for dfa in 'svwn blyp b3lyp'.split():
            print(f'    Test mol {mol.name()}:  dfa={dfa}')
            E_ref, dfa_wfn = psi4.energy(dfa, return_wfn=True)
            E_calc = scf_losc(dfa_info[dfa], dfa_wfn).energy()
            my_E_calc = my_scf_losc(dfa_info[dfa], dfa_wfn).energy()
            self.assertAlmostEqual(E_ref, E_calc, places=precision)
            self.assertAlmostEqual(my_E_calc, E_calc, places=precision)

    def test_uks_open_shell(self):
        """
        Test unrestricted calculations for open shell cases.
        """
        print("\n==> Test UKS Open Shell")
        optstash = psi4.driver.p4util.OptionsState(['SCF', 'REFERENCE'])
        psi4.core.set_local_option('SCF', 'REFERENCE', 'UHF')

        # LO = CO
        self.run_mol_no_correction(self.mol_H2_plus_1A)

        # TODO: LO != CO

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

        # TODO: LO != CO

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

        # TODO: LO != CO

        optstash.restore()


if __name__ == '__main__':
    unittest.main()
