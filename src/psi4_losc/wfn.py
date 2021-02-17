"""
LOSC wavefunction class that are derived from `psi4.core.RHF` and `psi4.core.UHF`.
With the inheritance of wavefunction class from psi4 and override
'virtual' functions of `psi4.core.HF` class, it becomes very easy to integrate
LOSC with psi4.
"""

import psi4
import numpy as np
from py_losc import py_losc


def _common_init(self, nspin, curvature, C_lo):
    """
    Common initialization of LOSC wavefunction.
    """
    if not curvature:
        raise ValueError(
            'Need LOSC curvature matrix to build a LOSC wavefunction.')
    if not C_lo:
        raise ValueError(
            'Need LOSC localized orbitals to build a LOSC wavefunction.')
    if len(curvature) != nspin:
        raise ValueError('Wrong size of spin of LOSC curvature matrix.')
    if len(C_lo) != nspin:
        raise ValueError('Wrong size of spin of LOSC localzied orbitals.')
    self._curvature = curvature
    self._C_lo = C_lo

    self._E_losc = 0


def _add_losc(self, nspin):
    """
    1. Build LOSC effective Hamiltonian and add it to internal wavefunction Fock
    matrix.
    2. Build LOSC correction to total energy
    """
    # call psi4.core.RHF.form_F() to build the DFA Fock matrix.
    S = np.asarray(self.S())
    D = [np.asarray(self.Da()), np.asarray(self.Db())]
    F = [np.asarray(self.Fa()), np.asarray(self.Fb())]
    # ==> add LOSC effective Fock matrix <==
    self._E_losc = 0
    for s in range(nspin):
        # build LOSC local occupation matrix
        local_occ = py_losc.local_occupation(self._C_lo[s], S, D[s])
        # build LOSC effective Fock matrix
        H_losc = py_losc.ao_hamiltonian_correction(
            S, self._C_lo[s], self._curvature[s], local_occ)
        F[s][:] += H_losc
        # form LOSC energy correction
        self._E_losc += py_losc.energy_correction(
            self._curvature[s], local_occ)
    if nspin == 1:
        self._E_losc *= 2


class RLOSC(psi4.core.RHF):
    def __init__(self, wfn, superfunc, curvature, C_lo):
        # Call base class constructor.
        super().__init__(wfn, superfunc)
        # Some common init for LOSC wfn.
        _common_init(self, 1, curvature, C_lo)

    def form_F(self):
        # Build DFA Fock matrix. This will update `wfn.Fa()`` and `wfn.Fb()``.
        super().form_F()
        # Add LOSC Hamiltonian into wfn internal Fock matrix.
        _add_losc(self, 1)

    def compute_E(self):
        # Compute DFA total energy
        E_dfa = super().compute_E()
        # Add LOSC energy. Register LOSC energy into wfn energetics table.
        self.set_energies('LOSC energy', self._E_losc)
        E_tot = E_dfa + self._E_losc
        return E_tot


class ULOSC(psi4.core.UHF):
    def __init__(self, wfn, superfunc, curvature, C_lo):
        # Call base class constructor.
        super().__init__(wfn, superfunc)
        # Some common init for LOSC wfn.
        _common_init(self, 2, curvature, C_lo)

    def form_F(self):
        # Build DFA Fock matrix. This will update `wfn.Fa()`` and `wfn.Fb()``.
        super().form_F()
        # Add LOSC Hamiltonian into wfn internal Fock matrix.
        _add_losc(self, 2)

    def compute_E(self):
        # Compute DFA total energy
        E_dfa = super().compute_E()
        # Add LOSC energy. Register LOSC energy into wfn energetics table.
        self.set_energies('LOSC energy', self._E_losc)
        E_tot = E_dfa + self._E_losc
        return E_tot
