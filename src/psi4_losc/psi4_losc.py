"""
Integrate post-SCF-LOSC and SCF-LOSC (frozen-LO) calculation with psi4 package.

Notes
-----
This module ONLY suuports the calculations of integer systems. To perform
calculations of fractional systems in order to test the delocalization error,
you can use the extended SCF procedure provided in `psi4_losc/scf.py` module.
"""

import psi4
import numpy as np
from py_losc import py_losc
from psi4_losc import utils
from psi4.driver.p4util.exceptions import ValidationError
from pkg_resources import parse_version
from qcelemental import constants

if parse_version(psi4.__version__) >= parse_version('1.3a1'):
    build_superfunctional = psi4.driver.dft.build_superfunctional
else:
    build_superfunctional = psi4.driver.dft.build_superfunctional

B3LYP = py_losc.DFAInfo(0.8, 0.2, 'B3LYP')
SVWN = py_losc.DFAInfo(1.0, 0, 'SVWN')
BLYP = py_losc.DFAInfo(1.0, 0, 'BLYP')
PBE = py_losc.DFAInfo(1.0, 0, 'PBE')
PBE0 = py_losc.DFAInfo(0.75, 0.25, 'PBE0')


def _validate_dfa_wfn(dfa_wfn):
    "Validate the input dfa wavefunction for the LOSC calculation."
    if not isinstance(dfa_wfn, psi4.core.HF):
        raise ValidationError('Not a psi4.core.HF object.')
    if dfa_wfn.molecule().schoenflies_symbol() != 'c1':
        raise ValidationError('LOSC only supports C1 symmetry')
    supfunc = dfa_wfn.functional()
    if supfunc.is_x_lrc():
        raise ValidationError(
            'Current LOSC does not support range-separated exchange functional.')
    if supfunc.is_c_hybrid():
        raise ValidationError(
            'Sorry, LOSC does not support double hybrid functional.')
    if supfunc.is_meta():
        raise ValidationError('Sorry, LOSC does not support meta-GGA.')


def post_scf_losc(dfa_info, dfa_wfn, orbital_energy_unit='eV', verbose=1):
    """
    Perform the post-SCF-LOSC calculation based on a DFA wavefunction.

    This function supports post-SCF LOSC calculations for integer/fractional
    systems with aufbau/non-aufbau occupations.
    If you want to calculate integer system with aufbau occupations, use
    `psi4.energy()` or `psi4_losc.scf()` to get the `dfa_wfn`.
    If you want to calculate integer system with non-aufbau occupation, or
    fractional system with aufbau/non-aufbau occupations, use `psi4_losc.scf()`
    to get the `dfa_wfn`.

    Parameters
    ----------
    dfa_info: py_losc.DFAInfo
        The information of the parent DFA, including the weights of exchanges.
    dfa_wfn: psi4.core.HF
        The converged wavefunction from a parent DFA. At return, it will have
        some new dynamic attributes. See the notes.
    orbital_energy_unit: str, default to 'eV'
        The units of orbital energies used to print in the output.
        Valid choices are ['au', 'eV'].
        'au': atomic unit, hartree.
        'eV': electronvolt.
    verbose: int, default to 1
        print level. 0 means print nothing. 1 means normal print level. A larger
        number means more details.

    Returns
    -------
    energy: float
        The total energy of LOSC-DFA.
    eig: [[float, ...], ...], shape=(nspin, nbasis).
        The orbital energies from LOSC-DFA. nspin=1 if it is RKS and nspin=2
        if it is UKS.

    See Also
    --------
    py_losc.DFAInfo(): constructor of the DFA information class.
    psi4.energy(): return a DFA SCF wavefunction. psi4.energy() only supports
        calculations for integer systems with aufbau occupations.
    psi4_losc.scf(): return a DFA SCF wavefunction. psi4_losc.scf() supports
        calculations for integer/fractional systems with aufbau/non-aufbau
        occupations.

    Notes
    -----
    1. Ideally, we should be able to extract the weights of exchanges from the
    psi4 superfunctional objects which can be accessed from psi4 wfn. However,
    it looks like psi4 does not support this functionality well. So we require
    the user take the responsibility to construct the py_losc.DFAInfo object
    manually by himself. We provide several py_losc.DFAInfo objects for common
    DFAs in this module.

    2. At the return, the attribute of dfa_wfn, `dfa_wfn.losc_data` is updated
    or created if not exists:
        dfa_wfn.losc_data: dict
            The data that can be used for SCF-LOSC (frozen-LO) calculations.

            dfa_wfn.losc_data['occ']: dict
                The occupation information. See `psi4_losc.scf()`.
                If It exists, no changes, otherwise, it will be created with
                an empty dict which represents the aufbau occupation numbers.
            dfa_wfn.losc_data['curvature']: [np.array, ...]
                A list (size=1 for RKS and 2 for UKS) of LOSC curvature matrix.
                It will be created.
            dfa_wfn.losc_data['C_lo']: [np.array, ...]
                A list (size=1 for RKS and 2 for UKS) of LOSC LO coefficient
                matrix. It will be created.
    This is just for the implementation purpose. Users should not rely on this.
    """
    def local_print(level, *args):
        if verbose >= level:
            t = [f'{i}' for i in args]
            psi4.core.print_out(f'{" ".join(t)}\n')

    # sanity-check of input dfa wfn.
    _validate_dfa_wfn(dfa_wfn)
    if orbital_energy_unit not in ['eV', 'au']:
        raise Exception(
            'Invalid input for "orbital_energy_unit". Choices are ["eV", "au"].')

    # print header
    local_print(1, "--------------------------------")
    local_print(1, "          post-SCF-LOSC")
    local_print(1, "          by Yuncai Mei")
    local_print(1, "--------------------------------")

    eig_factor = 1.0 if orbital_energy_unit == 'au' else constants.hartree2ev
    is_rks = dfa_wfn.same_a_b_orbs() and dfa_wfn.same_a_b_dens()
    nspin = 1 if is_rks else 2
    nbf = dfa_wfn.basisset().nbf()

    # map needed matrices to DFA wfn.
    C_co = [np.asarray(dfa_wfn.Ca()), np.asarray(dfa_wfn.Cb())]
    H_ao = [np.asarray(dfa_wfn.Fa()), np.asarray(dfa_wfn.Fb())]
    D_ao = [np.asarray(m) for m in dfa_wfn.mintshelper().ao_dipole()]
    S = np.asarray(dfa_wfn.S())
    D = [np.asarray(dfa_wfn.Da()), np.asarray(dfa_wfn.Db())]

    # ====> !!! Start of LOSC !!! <====
    # ==> LOSC localization <==
    C_lo = [None] * nspin
    for s in range(nspin):
        # create losc localizer object
        localizer = py_losc.LocalizerV2(C_co[s], H_ao[s], D_ao)
        # compute LOs
        C_lo[s] = localizer.lo()

    # ==> LOSC curvature matrix <==
    # build matrices related to density fitting
    df_pmn, df_Vpq_inv = utils.form_df_basis_matrix(dfa_wfn)
    df_pii = [None] * nspin
    for s in range(nspin):
        # build three-center integral <fitbasis|lo, lo>
        df_pii[s] = np.einsum('pmn,mi,ni->pi', df_pmn,
                              C_lo[s], C_lo[s], optimize=True)

    # build weights of grid points
    grid_w = utils.form_grid_w(dfa_wfn)

    # build values of LOs on grid points
    grid_lo = [utils.form_grid_lo(dfa_wfn, C_lo_) for C_lo_ in C_lo]

    # build LOSC curvature matrices
    curvature = [None] * nspin
    for s in range(nspin):
        # build losc curvature matrix
        curvature_helper = py_losc.CurvatureV2(
            dfa_info, df_pii[s], df_Vpq_inv, grid_lo[s], grid_w)
        curvature[s] = curvature_helper.kappa()

    # ==> LOSC local occupation matrix <==
    local_occ = [None] * nspin
    for s in range(nspin):
        # build losc local occupation matrix
        local_occ[s] = py_losc.local_occupation(C_lo[s], S, D[s])

    # ==> LOSC corrections <==
    H_losc = [None] * nspin
    E_losc = [None] * nspin
    losc_eig = [None] * nspin
    for s in range(nspin):
        # build losc effective Hamiltonian matrix
        H_losc[s] = py_losc.ao_hamiltonian_correction(
            S, C_lo[s], curvature[s], local_occ[s])
        # calculate losc energy correction
        E_losc[s] = py_losc.energy_correction(curvature[s], local_occ[s])
        # calculate corrected orbital energy from LOSC.
        losc_eig[s] = np.array(py_losc.orbital_energy_post_scf(
            H_ao[s], H_losc[s], C_co[s])) * eig_factor

    E_losc_tot = 2 * E_losc[0] if nspin == 1 else sum(E_losc)
    # ====> !!! End of LOSC !!! <====

    # ==> add/update LOSC related attributes to dfa_wfn <==
    if not hasattr(dfa_wfn, 'losc_data'):
        dfa_wfn.losc_data = {}
    if 'occ' not in dfa_wfn.losc_data:
        dfa_wfn.losc_data['occ'] = {}
    dfa_wfn.losc_data['curvature'] = curvature
    dfa_wfn.losc_data['C_lo'] = C_lo

    # ==> Print energies to output <==
    E_losc_dfa_tot = dfa_wfn.energy() + E_losc_tot
    local_print(1, '\n=> Total Energy <=')
    local_print(1, 'DFA total energy: {:.16f}'.format(dfa_wfn.energy()))
    local_print(1, 'LOSC correction energy: {:.16f}'.format(E_losc_tot))
    local_print(1, 'LOSC-DFA total energy: {:.16f}'.format(E_losc_dfa_tot))

    # print orbital energies
    dfa_eigs = [np.asarray(dfa_wfn.epsilon_a()) * eig_factor,
                np.asarray(dfa_wfn.epsilon_b()) * eig_factor]
    _, occ_idx, occ_val = utils.form_occ(dfa_wfn, dfa_wfn.losc_data['occ'])
    occ_idx_val = [dict(zip(occ_idx[s], occ_val[s])) for s in range(nspin)]

    local_print(1, '\n=> Orbital Energy <=')
    for s in range(nspin):
        if not is_rks:
            local_print(1, f"{'Alpha' if s == 0 else 'Beta'} orbital energy:")
        local_print(1, "{:<5s}  {:<8s}  {:>14s} {:>14s}"
                    .format("Index", "Occ", f"DFA ({orbital_energy_unit})",
                            f"LOSC ({orbital_energy_unit})"))
        for i in range(nbf):
            local_print(1, "{:<5d}  {:<8.5f}  {:>14.6f} {:>14.6f}"
                        .format(i, occ_idx_val[s].get(i, 0), dfa_eigs[s][i],
                                losc_eig[s][i]))
        local_print(1, "")
    return E_losc_dfa_tot, losc_eig


def scf_losc(dfa_info, dfa_wfn):
    """
    Perform the SCF-LOSC (frozen-LO) calculation based on a DFA wavefunction.

    This function calls `psi4.energy()` to do the SCF calculation internally.
    To make the usage of `psi4.energy()` compatibale with LOSC calculation,
    we define the LOSC wavefunction (see `wfn.RLOSC` and `wfn.ULOSC`) that are
    derived from `psi4.core.RHF` or `psi4.core.UHF`. Simply overwrite several
    key functions derived from `psi4.core.HF`, we can make `psi4.energy()`
    function involve the LOSC contributions correctly in the SCF procedure.
    See `wfn.py`.

    This function does not supports calculations for fractional systems.

    Parameters
    ----------
    dfa_info: py_losc.DFAInfo
        The information of the parent DFA, including the weights of exchanges.
    dfa_wfn: psi4.core.HF like psi4 object
        The converged wavefunction from a parent DFA.

    Returns
    -------
    wfn: psi4_losc.wfn.RLOSC or psi4_losc.wfn.ULOSC
        The LOSC wavefunction.

    Notes
    -----
    1. Since `psi4.energy()` is used internally to drive the SCF procedure,
    ideally, this function supports all the types of SCF calculations that
    are supported by psi4, such as integer/aufbau systems, MOM calculations.
    However, this is not fully tested. But for normal SCF, meaning integer
    and aufbau system, this function works fine.
    """
    # Check if the user tries to customize the occupation number.
    if hasattr(dfa_wfn, 'losc_data'):
        if 'occ' in dfa_wfn.losc_data:
            raise Exception('Customizing occupation number is allowed.')

    # Do post-SCF-LOSC to build curvature and LO.
    post_scf_losc(dfa_info, dfa_wfn)

    # Do SCF-LOSC.
    dfa_name = dfa_wfn.functional().name()
    _, wfn = psi4.energy(
        dfa_name, losc_data=dfa_wfn.losc_data, return_wfn=True)

    # Remove dynamical attributes created for LOSC.
    if hasattr(wfn, 'losc_data'):
        delattr(wfn, 'losc_data')
    return wfn
