"""Import this module will update function `psi4.proc.scf_wavefunction_factory()`
defined in psi4 to be an extended version. Making such update is to enable
`psi4.energy()` to be compatible with LOSC calculation.
"""

import psi4
import py_losc
import numpy as np

# Function `psi4.proc.scf_wavefunction_factory()` will be extended in this
# module. Here, we save the original one first.
_psi4_scf_wavefunction_factory = psi4.driver.scf_wavefunction_factory


def _to_losc_wfn(wfn, losc_data):
    """
    Update psi4 wfn methods to have LOSC contributions by overriding some
    key functions.
    """
    self = wfn
    original_form_F = self.form_F
    origninal_compute_E = self.compute_E
    nspin = 1 if wfn.same_a_b_orbs() and wfn.same_a_b_dens() else 2
    curvature = losc_data.get('curvature', [])
    C_lo = losc_data.get('C_lo', [])
    E_losc = [0] # we need to use list to let inner function has modify access.

    def form_F():
        # Build DFA Fock matrix. This will update `wfn.Fa()` and `wfn.Fb()`.
        original_form_F()

        # ==> LOSC contributions <==
        # 1. Build LOSC effective Hamiltonian and add it to internal
        # wavefunction Fock matrix.
        # 2. Build LOSC correction to total energy
        S = np.asarray(self.S())
        D = [np.asarray(self.Da()), np.asarray(self.Db())]
        F = [np.asarray(self.Fa()), np.asarray(self.Fb())]
        E_losc[0] = 0
        for s in range(nspin):
            # build LOSC local occupation matrix
            local_occ = py_losc.local_occupation(C_lo[s], S, D[s])
            # build LOSC effective Fock matrix
            H_losc = py_losc.ao_hamiltonian_correction(
                S, C_lo[s], curvature[s], local_occ)
            F[s][:] += H_losc
            # form LOSC energy correction
            E_losc[0] += py_losc.energy_correction(curvature[s], local_occ)
        if nspin == 1:
            E_losc[0] *= 2

    def compute_E():
        # Compute DFA total energy
        E_dfa = origninal_compute_E()
        # Add LOSC energy. Register LOSC energy into wfn energetics table.
        self.set_energies('LOSC energy', E_losc[0])
        E_tot = E_dfa + E_losc[0]
        return E_tot

    # override wfn methods.
    self.form_F = form_F
    self.compute_E = compute_E


def _scf_wavefunction_factory_extended_version(name, ref_wfn, reference, **kwargs):
    """
    Build an SCF wavefunction. This is the extended version for
    `psi4.proc.scf_wavefunction_factory()`.
    """
    # build psi4 SCF wfn.
    reference = reference.upper()
    wfn = _psi4_scf_wavefunction_factory(name, ref_wfn, reference, **kwargs)

    # update to LOSC wavefunction behavior if it is necessary.
    losc_data = kwargs.get('losc_data', {})
    if losc_data:
        # same sanity checks
        curvature = losc_data.get('curvature', [])
        C_lo = losc_data.get('C_lo', [])
        if not curvature:
            raise Exception('need LOSC curvature matrix to build LOSC class.')
        if not C_lo:
            raise Exception('need LOSC localized orbital to build LOSC class.')
        # restricted
        if reference in ["RHF", "RKS"]:
            if len(curvature) != 1 or len(C_lo) != 1:
                raise Exception(
                    'restricted case: size of curvature or C_lo is not 1.')
        # unrestricted
        elif reference in ['UHF', 'UKS']:
            if len(curvature) != 2 or len(C_lo) != 2:
                raise Exception(
                    'unrestricted case: size of curvature or C_lo is not 2.')
        # error reference
        else:
            raise Exception(
                f"reference ({reference}) not supported for LOSC wavefunction.")

        # update to losc wfn.
        _to_losc_wfn(wfn, losc_data)
    return wfn


# Here, we update `psi4.proc.scf_wavefunction_factory()` with the extended
# version.
psi4.driver.scf_wavefunction_factory = _scf_wavefunction_factory_extended_version
# we need to update psi4.proc.scf_wavefunction_factory,
# not psi4.driver.scf_wavefunction_factory.
psi4.proc.scf_wavefunction_factory = _scf_wavefunction_factory_extended_version
