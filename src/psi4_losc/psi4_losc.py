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
import psi4_losc.options as options
from py_losc import py_losc
from psi4_losc import utils
from qcelemental import constants

B3LYP = py_losc.DFAInfo(0.8, 0.2, 'B3LYP')
SVWN = py_losc.DFAInfo(1.0, 0, 'SVWN')
BLYP = py_losc.DFAInfo(1.0, 0, 'BLYP')
PBE = py_losc.DFAInfo(1.0, 0, 'PBE')
PBE0 = py_losc.DFAInfo(0.75, 0.25, 'PBE0')


def _validate_dfa_wfn(dfa_wfn):
    "Validate the input dfa wavefunction for the LOSC calculation."
    # check type
    if not isinstance(dfa_wfn, psi4.core.HF):
        raise Exception('Not a psi4.core.HF object.')
    # check symmetry
    if dfa_wfn.molecule().schoenflies_symbol() != 'c1':
        raise Exception('LOSC only supports C1 symmetry')
    # check spin
    is_rks = psi4.core.get_option('SCF', 'REFERENCE') in ['RKS', 'RHF']
    is_rks_wfn = dfa_wfn.same_a_b_orbs() and dfa_wfn.same_a_b_dens()
    if is_rks != is_rks_wfn:
        raise Exception(
            'Reference in passed wfn is different to psi4 reference setting.')
    # check super functional
    supfunc = dfa_wfn.functional()
    if supfunc.is_x_lrc():
        raise Exception(
            'Current LOSC does not support range-separated exchange functional.')
    if supfunc.is_c_hybrid():
        raise Exception(
            'Sorry, LOSC does not support double hybrid functional.')
    if supfunc.is_meta():
        raise Exception('Sorry, LOSC does not support meta-GGA.')


def post_scf_losc(dfa_info, dfa_wfn, orbital_energy_unit='eV', verbose=1,
                  return_losc_data=False,
                  window=None):
    """
    Perform the post-SCF-LOSC calculation based on a DFA wavefunction.

    This function supports post-SCF LOSC calculations for integer/fractional
    systems with aufbau/non-aufbau occupations:
    (1) If you want to calculate integer system with aufbau occupations, use
    `psi4.energy()` to get the input `dfa_wfn`.
    (2) If you want to calculate integer/non-aufbau system or fractional system
    (either aufbau or non-aufbau occupation), use `psi4_losc.scf.scf()`
    to get the input `dfa_wfn`.

    Parameters
    ----------
    dfa_info: py_losc.DFAInfo
        The information of the parent DFA, including the weights of exchanges.
    dfa_wfn: psi4.core.HF
        The converged wavefunction from a parent DFA.
    orbital_energy_unit: str, default to 'eV'
        The units of orbital energies used to print in the output.
        Valid choices are ['au', 'eV'].
        'au': atomic unit, hartree.
        'eV': electronvolt.
    verbose: int, default to 1
        print level. 0 means print nothing. 1 means normal print level. A larger
        number means more details.
    return_losc_data: bool, default to false.
        Return the data of LOSC or not.
    window: [float, float], default to None.
        The orbital energy window in eV to select COs to do localization.
        Select all COs whose energies are in the window. Default to None which
        means dismiss the window setting and use all COs to do localization.

    Returns
    -------
    energy: float
        The total energy of LOSC-DFA.
    eig: [np.array, ...]
        All orbital energies (same number to basis set) from LOSC-DFA. For RKS,
        It includes alpha orbital energies. For UKS, it includes alpha and beta
        orbital energies in order.
    losc_data: dict
        Returned if `return_losc_data` is true. It contains the data of LOSC
        calculations.

    See Also
    --------
    py_losc.DFAInfo(): constructor of the DFA information class.
    psi4.energy(): return a DFA SCF wavefunction. psi4.energy() only supports
        calculations for integer systems with aufbau occupations.
    psi4_losc.scf.scf(): return a DFA SCF wavefunction. psi4_losc.scf.scf()
        supports calculations for integer/fractional systems with
        aufbau/non-aufbau occupations.

    Notes
    -----
    1. Ideally, we should be able to extract the weights of exchanges from the
    psi4 superfunctional objects which can be accessed from psi4 wfn. However,
    it looks like psi4 does not support this functionality well. So we require
    the user take the responsibility to construct the py_losc.DFAInfo object
    manually by himself. We provide several `py_losc.DFAInfo` objects for common
    DFAs in `psi4_losc` module.
    """
    # sanity-check of input dfa wfn.
    _validate_dfa_wfn(dfa_wfn)
    if orbital_energy_unit not in ['eV', 'au']:
        raise Exception(
            'Invalid input for "orbital_energy_unit". Choices are ["eV", "au"].')

    local_print = utils.init_local_print(verbose)
    # print header
    local_print(1, "--------------------------------")
    local_print(1, "          post-SCF-LOSC")
    local_print(1, "          by Yuncai Mei")
    local_print(1, "--------------------------------")

    is_rks = dfa_wfn.same_a_b_orbs() and dfa_wfn.same_a_b_dens()
    nspin = 1 if is_rks else 2

    # Need real number of electrons. Cannot rely on `wfn.nalpha()` and
    # `wfn.nbeta()`, because `post_scf_losc()` can support fractional systems.
    occ = {}
    if hasattr(dfa_wfn, 'losc_data'):
        occ = dfa_wfn.losc_data['occ']
    _, _, occ_val = utils.form_occ(dfa_wfn, occ)
    nelec = [sum(x) for x in occ_val]

    # map needed matrices to DFA wfn.
    mintshelper = psi4.core.MintsHelper(dfa_wfn.basisset())
    C_co = [np.asarray(dfa_wfn.Ca()), np.asarray(dfa_wfn.Cb())]
    H_ao = [np.asarray(dfa_wfn.Fa()), np.asarray(dfa_wfn.Fb())]
    D_ao = [np.asarray(m) for m in mintshelper.ao_dipole()]
    S = np.asarray(dfa_wfn.S())
    D = [np.asarray(dfa_wfn.Da()), np.asarray(dfa_wfn.Db())]

    # ====> !!! Start of LOSC !!! <====
    if verbose >= 1:
        options.show_options()
    # ==> LOSC localization <==

    def select_CO(wfn, spin, window):
        """Return the CO index bounds."""
        eig = [np.asarray(wfn.epsilon_a()), np.asarray(wfn.epsilon_b())][spin]
        nbf = wfn.basisset().nbf()
        # return all COs.
        if not window or nelec[spin] == 0:
            return None
        # select COs
        if len(window) != 2:
            raise Exception('Invalid LOSC window: wrong size of window.')
        if window[0] >= window[1]:
            raise ValueError('Invalid LOS window: left bound >= right bound.')
        idx_start = next(filter(
            lambda x: x[1] * constants.hartree2ev >= window[0], enumerate(eig)), [nbf])[0]
        idx_end = next(filter(
            lambda x: x[1] * constants.hartree2ev >= window[1], enumerate(eig)), [nbf])[0]
        if idx_end - idx_start <= 0:
            raise Exception(
                'LOSC window is too tight. No COs selected to do LOSC localization.')
        return (idx_start, idx_end)

    # Get selected COs index from the window setting.
    local_print(1, '\n==> LOSC Localization Window <==')
    selected_co_idx = [None] * nspin
    for s in range(nspin):
        idx = select_CO(dfa_wfn, s, window)
        if not idx:
            local_print(1, f'localization COs index (spin={s}): ALL COs')
        else:
            local_print(
                1, f'localization COs index (spin={s}): [{idx[0]}, {idx[1]})')
        selected_co_idx[s] = idx
    local_print(1, '')

    C_lo = [None] * nspin
    U = [None] * nspin
    for s in range(nspin):
        # Get selected COs.
        if selected_co_idx[s]:
            idx_start, idx_end = selected_co_idx[s]
            C_co_used = C_co[s][:, list(range(idx_start, idx_end))]
        else:
            C_co_used = C_co[s]
        # create losc localizer object
        if options.localization['version'] == 2:
            localizer = py_losc.LocalizerV2(C_co_used, H_ao[s], D_ao)
            localizer.set_c(options.localization['v2_parameter_c'])
            localizer.set_gamma(options.localization['v2_parameter_gamma'])
        else:
            raise Exception('Detect non-supporting localization version.')
        localizer.set_max_iter(options.localization['max_iter'])
        localizer.set_convergence(options.localization['convergence'])
        localizer.set_random_permutation(
            options.localization['random_permutation'])
        # compute LOs
        C_lo[s], U[s] = localizer.lo_U()

    # ==> LOSC curvature matrix <==
    # build matrices related to density fitting
    df_pii, df_Vpq_inv = utils.form_df_matrix(dfa_wfn, C_lo)

    # build weights of grid points
    grid_w = utils.form_grid_w(dfa_wfn)

    # build values of LOs on grid points
    grid_lo = [utils.form_grid_lo(dfa_wfn, C_lo_) for C_lo_ in C_lo]

    # build LOSC curvature matrices
    curvature = [None] * nspin
    for s in range(nspin):
        # build losc curvature matrix
        if options.curvature['version'] == 2:
            curvature_helper = py_losc.CurvatureV2(dfa_info, df_pii[s],
                                                   df_Vpq_inv, grid_lo[s],
                                                   grid_w)
            curvature_helper.set_tau(options.curvature['v2_parameter_tau'])
            curvature_helper.set_zeta(options.curvature['v2_parameter_zeta'])
        elif options.curvature['version'] == 1:
            curvature_helper = py_losc.CurvatureV1(dfa_info, df_pii[s],
                                                   df_Vpq_inv, grid_lo[s],
                                                   grid_w)
            curvature_helper.set_tau(options.curvature['v1_parameter_tau'])
        else:
            raise Exception('Detect un-supported curvature version.')
        # TODO: add curvature setting options.
        curvature[s] = curvature_helper.kappa()

    # ==> LOSC local occupation matrix <==
    local_occ = [None] * nspin
    for s in range(nspin):
        # build losc local occupation matrix
        local_occ[s] = py_losc.local_occupation(C_lo[s], S, D[s])

    # ==> print matrix into psi4 output file <==
    if verbose >= 2:
        local_print(2, '\n==> Details of Matrices in LOSC <==')
        for s in range(nspin):
            local_print(2, f'\nCO coefficient matrix.T: spin={s}')
            utils.print_full_matrix(C_co[s].T)
        for s in range(nspin):
            local_print(2, f'\nLO coefficient matrix.T: spin={s}')
            utils.print_full_matrix(C_lo[s].T)
        for s in range(nspin):
            local_print(2, f'\nLO U matrix: spin={s}')
            utils.print_full_matrix(U[s])
        for s in range(nspin):
            local_print(2, f'\nCurvature: spin={s}')
            utils.print_sym_matrix(curvature[s])
        for s in range(nspin):
            local_print(2, f'\nLocal Occupation matrix: spin={s}')
            utils.print_sym_matrix(local_occ[s])

    # ==> LOSC corrections <==
    H_losc = [None] * nspin
    E_losc = [None] * nspin
    losc_eigs = [None] * nspin
    eig_factor = 1.0 if orbital_energy_unit == 'au' else constants.hartree2ev
    for s in range(nspin):
        # build losc effective Hamiltonian matrix
        H_losc[s] = py_losc.ao_hamiltonian_correction(
            S, C_lo[s], curvature[s], local_occ[s])
        # calculate losc energy correction
        E_losc[s] = py_losc.energy_correction(curvature[s], local_occ[s])
        # calculate corrected orbital energy from LOSC.
        losc_eigs[s] = np.array(py_losc.orbital_energy_post_scf(
            H_ao[s], H_losc[s], C_co[s])) * eig_factor

    E_losc_tot = 2 * E_losc[0] if is_rks else sum(E_losc)
    E_losc_dfa_tot = dfa_wfn.energy() + E_losc_tot
    # ====> !!! End of LOSC !!! <====

    # ==> pack LOSC data into dict <==
    dfa_eigs = [np.asarray(dfa_wfn.epsilon_a()) * eig_factor,
                np.asarray(dfa_wfn.epsilon_b()) * eig_factor]
    occ = {}
    if hasattr(dfa_wfn, 'losc_data'):
        occ = dfa_wfn.losc_data.get('occ', {})
    losc_data = {}
    losc_data['occ'] = occ
    losc_data['losc_type'] = 'post-SCF-LOSC'
    losc_data['orbital_energy_unit'] = orbital_energy_unit
    losc_data['nspin'] = nspin
    losc_data['curvature'] = curvature
    losc_data['C_lo'] = C_lo
    losc_data['dfa_energy'] = dfa_wfn.energy()
    losc_data['dfa_orbital_energy'] = dfa_eigs
    losc_data['losc_energy'] = E_losc_tot
    losc_data['losc_dfa_energy'] = E_losc_tot + dfa_wfn.energy()
    losc_data['losc_dfa_orbital_energy'] = losc_eigs

    # ==> Print energies to output <==
    utils.print_total_energies(1, losc_data)
    utils.print_orbital_energies(1, dfa_wfn, losc_data)

    if return_losc_data:
        return E_losc_dfa_tot, losc_eigs, losc_data
    else:
        return E_losc_dfa_tot, losc_eigs


def scf_losc(dfa_info, dfa_wfn, orbital_energy_unit='eV', verbose=1):
    """
    Perform the SCF-LOSC (frozen-LO) calculation based on a DFA wavefunction.

    This function use `psi4.energy()` to do the SCF procedure.
    This function only supports calculations for integer systems.

    Parameters
    ----------
    dfa_info: py_losc.DFAInfo
        The information of the parent DFA, including the weights of exchanges.
    dfa_wfn: psi4.core.HF like psi4 object
        The converged wavefunction from a parent DFA.
    orbital_energy_unit: str, default to 'eV'
        The units of orbital energies used to print in the output.
        Valid choices are ['au', 'eV'].
        'au': atomic unit, hartree.
        'eV': electronvolt.
    verbose: int
        The print level to control `psi4_losc.post_scf_losc()` and
        `psi4.energy()`.

    Returns
    -------
    wfn: psi4.core.RHF or psi4.core.HHF with LOSC behaviors.

    See Also
    --------
    psi4_losc.build_scf_wfn.py: update psi4 wfn objects to have LOSC hebaviors.

    Notes
    -----
    1. The `psi4.energy()` function is used internally to drive the SCF
    procedure. So, ideally, this function supports all the types of SCF
    calculations that are supported by psi4, such as integer/aufbau systems,
    MOM calculations. However, these are not fully tested. But for normal SCF,
    meaning integer and aufbau system, this function works fine.

    2. SCF-LOSC (frozen-LO) requires to use the DFA wfn as the initial guess.
    The psi4 guess setting for SCF will be ignored. To use DFA wfn as the
    initial guess, currently it is limited by `psi4.energy()` interface in
    which directly passing the wfn as an argument is rejected. So, the only
    choice is to write the DFA wfn into scratch file, then let psi4 use
    `guess read` to set up the initial guess. So, calling this function will
    generate a wfn scratch file! You need to take care of it at the exit.
    """
    # sanity-check of input dfa wfn.
    _validate_dfa_wfn(dfa_wfn)
    if orbital_energy_unit not in ['eV', 'au']:
        raise Exception(
            'Invalid input for "orbital_energy_unit". Choices are ["eV", "au"].')

    # Check if the user tries to customize the occupation number.
    if hasattr(dfa_wfn, 'losc_data'):
        if 'occ' in dfa_wfn.losc_data:
            raise Exception('Customizing occupation number is allowed.')

    local_print = utils.init_local_print(verbose)

    # Do post-SCF-LOSC to build curvature and LO.
    _, _, losc_data = post_scf_losc(dfa_info, dfa_wfn, verbose=verbose,
                                    return_losc_data=True,
                                    orbital_energy_unit=orbital_energy_unit)

    # Do SCF-LOSC.
    # Trying to use ref_wfn keyword in psi4.energy(ref_wfn=dfa_wfn) will cause
    # an error. This issue may be solved in later version of psi4. Currently,
    # the way to let psi4 use a ref_wfn as initial guess is to save the wfn
    # into a file then use psi4 `guess read` option.
    optstash = psi4.driver.p4util.OptionsState(
        ['SCF', 'GUESS'],
        ['SCF', 'PRINT'])

    # update local settings.
    psi4.core.set_local_option('SCF', 'GUESS', 'READ')
    psi4.core.set_local_option('SCF', 'PRINT', verbose)

    # No idea what does 180 mean. Anyways, this is the just the scratch file
    # naming conversion in psi4. The wfn file is binary in which the wfn
    # data is represented as python dict and saved by calling `numpy.save()`.
    # Calling `wfn.get_scratch_filename()` may not return the full file name
    # with an suffix of `.npy`. So we need to manually check it. Again, this
    # is just the convension in psi4.
    ref_wfn_file = dfa_wfn.get_scratch_filename(180)
    if not ref_wfn_file.endswith('.npy'):
        ref_wfn_file += '.npy'
    dfa_wfn.to_file(ref_wfn_file)
    local_print(
        1, '\nNotice: psi4 guess setting is ignored for the SCF-LOSC calculation.')
    local_print(
        1, 'The wfn from the associated DFA is ALWAYS used as the initial guess.')
    local_print(1, f'File of DFA wfn: {ref_wfn_file}')
    local_print(1, '')

    dfa_name = dfa_wfn.functional().name()
    _, wfn = psi4.energy(dfa_name, losc_data=losc_data, return_wfn=True)

    eig_factor = 1.0 if orbital_energy_unit == 'au' else constants.hartree2ev
    losc_data['losc_type'] = 'SCF-LOSC'
    losc_data['dfa_energy'] = wfn.energy() - wfn.get_energies('LOSC energy')
    losc_data['losc_energy'] = wfn.get_energies('LOSC energy')
    losc_data['losc_dfa_energy'] = wfn.energy()
    losc_data['losc_dfa_orbital_energy'] = [np.asarray(wfn.epsilon_a()) * eig_factor,
                                            np.asarray(wfn.epsilon_b()) * eig_factor]

    # ==> Print energies in LOSC style <==
    utils.print_total_energies(verbose, losc_data)
    utils.print_orbital_energies(verbose, wfn, losc_data)

    # Remove dynamical attributes created for LOSC.
    if hasattr(wfn, 'losc_data'):
        delattr(wfn, 'losc_data')

    optstash.restore()
    return wfn
