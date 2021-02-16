from psi4_losc.utils import form_occ
from psi4_losc.utils import is_aufbau_system
from psi4_losc.utils import is_integer_system
import psi4
import unittest


class TestFormOcc(unittest.TestCase):
    def setUp(self):
        psi4.core.be_quiet()
        self.mol_H = psi4.geometry("""
            0 2
            H
            symmetry c1
            """, name="H")

        self.mol_H2 = psi4.geometry("""
            0 1
            H
            H 1 1
            symmetry c1
            """, name="H2")

        self.wfn_H = psi4.core.Wavefunction.build(self.mol_H, '3-21G')
        self.wfn_H2 = psi4.core.Wavefunction.build(self.mol_H2, '3-21G')

    def test_key_error_spin_key(self):
        print('\n==> Test spin key error:')
        with self.assertRaises(Exception):
            form_occ(self.wfn_H, occ={'wrong_key': {}})

    def test_key_error_occ_key(self):
        print('\n==> Test occupation key error:')
        print('    check float occ index key.')
        with self.assertRaises(Exception):
            form_occ(self.wfn_H, occ={'alpha': {1.2: 1}})
        print('    check wrong str key.')
        with self.assertRaises(Exception):
            form_occ(self.wfn_H, occ={'alpha': {'wrong_key': 1}})
        print('    check out-of-range index key.')
        with self.assertRaises(Exception):
            form_occ(self.wfn_H, occ={'alpha': {-1: 1}})
        print('    check same appearance of HOMO and HOMO index.')
        with self.assertRaises(Exception):
            form_occ(self.wfn_H, occ={'alpha': {0: 1, 'homo': 1}})
        print('    check same appearance of LUMO and LUMO index.')
        with self.assertRaises(Exception):
            form_occ(self.wfn_H, occ={'alpha': {1: 1, 'lumo': 1}})
        print('    check special case without HOMO.')
        with self.assertRaises(Exception):
            form_occ(self.wfn_H, occ={'beta': {'homo': 1}})

    def test_val_error_occ_val(self):
        print('\n==> Test occupation value error:')
        print('    check occ number is out-of-range.')
        with self.assertRaises(Exception):
            form_occ(self.wfn_H, occ={'alpha': {0: 999}})
        with self.assertRaises(Exception):
            form_occ(self.wfn_H, occ={'alpha': {0: -1}})

    def test_cases(self):
        print('\n==> Test the correctness of occupation value:')
        print('    check case #1: electrons on both alpha and beta.')
        nocc, occ_idx, occ_val = form_occ(
            self.wfn_H2, occ={'alpha': {'homo': 0.5, 2: 0.8},
                              'beta': {'lumo': 0.7}})
        self.assertEqual(nocc, (2, 2))
        self.assertEqual(occ_idx, [(0, 2), (0, 1)])
        self.assertEqual(occ_val, [(0.5, 0.8), (1, 0.7)])

        print('    check case #2: with zero electrons')
        nocc, occ_idx, occ_val = form_occ(
            self.wfn_H2, occ={'alpha': {'homo': 0.0, 2: 0.8},
                              'beta': {'homo': 0.0}})
        self.assertEqual(nocc, (1, 0))
        self.assertEqual(occ_idx, [(2,), ()])
        self.assertEqual(occ_val, [(0.8,), ()])

class TestIsIntergerSystem(unittest.TestCase):
    def setUp(self):
        psi4.core.be_quiet()
        self.mol_H = psi4.geometry("""
            0 2
            H
            symmetry c1
            """, name="H")
        self.wfn_H = psi4.core.Wavefunction.build(self.mol_H, '3-21G')

    def test_cases(self):
        print('\n==> Test the correctness:')
        print('    Normal psi4 wavefunction object.')
        rst = is_integer_system(self.wfn_H)
        self.assertTrue(rst)
        print('    Set HOMO to 1. No effects.')
        rst = is_integer_system(self.wfn_H, occ={'alpha': {'homo': 1.0}})
        self.assertTrue(rst)
        print('    Set HOMO to fraction. N-delta system.')
        rst = is_integer_system(self.wfn_H, occ={'alpha': {'homo': 0.1}})
        self.assertFalse(rst)
        print('    Set LUMO to fraction. N+delta system.')
        rst = is_integer_system(self.wfn_H, occ={'alpha': {'lumo': 0.1}})
        self.assertFalse(rst)
        print('    Special case: just a proton without any electron')
        rst = is_integer_system(self.wfn_H, occ={'alpha': {'homo': 0}})
        self.assertTrue(rst)


class TestIsAufbauSystem(unittest.TestCase):
    def setUp(self):
        psi4.core.be_quiet()
        self.mol_H = psi4.geometry("""
            0 2
            H
            symmetry c1
            """, name="H")

        self.mol_H2 = psi4.geometry("""
            0 1
            H
            H 1 1
            symmetry c1
            """, name="H2")

        self.wfn_H = psi4.core.Wavefunction.build(self.mol_H, '3-21G')
        self.wfn_H2 = psi4.core.Wavefunction.build(self.mol_H2, '3-21G')

    def test_aufbau(self):
        print('\n==> Test the correctness of aufbau cases:')
        print('    Default psi4 wfn.')
        rst = is_aufbau_system(self.wfn_H)
        self.assertTrue(rst)
        print('    Set HOMO to 1.')
        rst = is_aufbau_system(self.wfn_H, occ={'alpha': {'homo': 1.0}})
        self.assertTrue(rst)
        print('    Set HOMO to fraction: aufbau N-delta case.')
        rst = is_aufbau_system(self.wfn_H, occ={'alpha': {'homo': 0.1}})
        self.assertTrue(rst)
        print('    Set HOMO to zero: aufbau N-1 case.')
        rst = is_aufbau_system(self.wfn_H, occ={'alpha': {'homo': 0.1}})
        self.assertTrue(rst)
        print('    Set LUMO to fraction: aufbau N+delta case.')
        rst = is_aufbau_system(self.wfn_H, occ={'alpha': {'lumo': 0.1}})
        self.assertTrue(rst)
        print('    Set LUMO to 1: aufbau N+1 case.')
        rst = is_aufbau_system(self.wfn_H, occ={'alpha': {'lumo': 1}})
        self.assertTrue(rst)
        print('    Special case: a naked proton.')
        rst = is_aufbau_system(self.wfn_H, occ={'alpha': {'homo': 0}})
        self.assertTrue(rst)

    def test_nonaufbau(self):
        print('\n==> Test the correctness of non-aufbau cases:')
        print('    Fractional excited state.')
        rst = is_aufbau_system(
            self.wfn_H2, occ={'alpha': {'homo': 0, 'lumo': 0.1}})
        self.assertFalse(rst)
        print('    Integer excited state #1.')
        rst = is_aufbau_system(self.wfn_H2, occ={'alpha': {2: 1}})
        self.assertFalse(rst)
        print('    Integer excited state #2.')
        rst = is_aufbau_system(self.wfn_H2, occ={'beta': {2: 1}})
        self.assertFalse(rst)


if __name__ == '__main__':
    unittest.main()
