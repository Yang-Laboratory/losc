"""
Integrate py_losc module with psi4 package to perform LOSC calculation.
"""

import numpy as np
from py_losc import py_losc
from psi4_losc import jk
from psi4_losc import diis
from psi4_losc import utils

import psi4
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


def _scf(name, guess_wfn=None, losc_ref_wfn=None, dfa_info=None, occ={},
         verbose=1, orbital_energy_unit='eV'):
    """
    Perform the SCF (normal SCF-DFA or SCF-LOSC-DFA) calculation for general
    systems.

    Parameters
    ----------
    name: str
        The name of the DFT functional or HF method. The same style as
        psi4.energy() function.
    guess_wfn: psi4.core.RHF or psi4.core.UHF, default to None.
        The wfn used to set the initial guess of the SCF calculation. This is
        equivalent to copy guess_wfn.C (coefficient matrix) and guess_wfn.D
        (density matrix) to the scf wavefunction object to initialize the SCF
        calculation. Note: (1) setting this variable will ignore the psi4
        global guess setting; (2) if you use this variable, make sure you are
        passing a reasonable guess wfn to the function. This function does not
        check the validity the input.
    losc_ref_wfn: psi4.core.RHF or psi4.core.UHF, default to None.
        The wfn used to do SCF-LOSC (frozen-LO) calculation. This variable is
        the flag to trigger SCF-LOSC (frozen-LO) calculation.
    dfa_info: py_losc.DFAInfo, default to None.
        DFA information of the weights of exchanges.
    occ: dict, default to an empty dict.
        A dictionary that specifies the customized occupation number. This
        variable will be ignored and overwrite as `losc_ref_wfn.losc_data['occ']`
        if the `losc_ref_wfn` is provided.
    verbose: int, default to 1
        print level to the psi4 output file.
        0 means print nothing. 1 means normal print level. A larger number means
        more details.

    Returns
    -------
    wfn: psi4.core.RHF or psi4.core.UHF
        The SCF wavefunction.
        1. Following members in the returned wfn object are updated:
            wfn.S(): AO overlap matrix.
            wfn.H(): Core matrix.
            wfn.Ca(), wfn.Cb(): CO coefficient matrix.
            wfn.Fa(), wfn.Fb(): Fock matrix.
            wfn.Da(), wfn.Db(): density matrix.
            wfn.Va(), wfn.Vb(): DFA (just the DFA, if it is LOSC-DFA) Vxc matrix.
            wfn.epsilon_a(), wfn.epsilon_b(): orbital energies.
            wfn.energy(): total energy.

        2. Other members are not touched. Accessing and using other members
        in the returned wfn may be a undefined behavior.

    Notes
    -----
    1. This function only relies on the GLOBAL settings in psi4, NOT any local
    settings for any psi4 modules!

    2. For the updated matrices/vectors in returned wfn (such as C, D, F matrices),
    although psi4.wfn objects provide interface to interact with the internal
    data (such the interface wfn.Fa() that returns a shared_ptr to internal Fock
    matrix), we do not map wfn internal matrices into python to modify these
    matrices in place in python. The reason is that the lifetime of the internal
    psi4.Matrix object managed by the shared_ptr is hard to track in python.
    Take wfn.Fa matrix as example.
    >>>
    # 1. Map internal Fock matrix into python. Now `F` in python refers to the
    # wfn.Fa matrix.
    F = np.asarray(wfn.Fa())

    # 2. psi4.core code internally changes the Fock matrix by doing something
    # like
    # wfn.Fa = std::make_shared(Matrix(nbf, nbf)); # this is in c++ psi4.core.

    # 3. Now `F` in python no longer refers to `wfn.Fa()`, since wfn.Fa is
    # updated.
    >>>
    Because of this issue, all the matrices those will be used to update wfn
    are allocated and managed in python. I know this doubles the memory cost
    (same matrix, such as Fock matrix, is allocated in both python and psi4.core),
    but this makes the logic clearer and less chances to have bug.
    At the return, the python matrix will be copied into wfn through the
    interface. I know we don't copy, but this is the price we have to pay.
    """
    def update_C(occ_idx, occ_val, Cin, C, Cocc, D):
        """
        Update C, Cocc and D matrices.
        """
        C[:] = Cin
        Cocc[:] = Cin[:, occ_idx]
        D[:] = np.einsum('i,ui,vi->uv', np.asarray(occ_val), Cocc, Cocc,
                         optimize=True)

    def is_aufbau(nocc, occ_idx):
        nspin = len(nocc)
        for s in range(nspin):
            if occ_idx[s] and occ_idx[s][-1] >= nocc[s]:
                return False
        return True

    def is_integer(occ_val):
        for v in occ_val:
            for x in v:
                if int(x) != float(x):
                    return False
        return True

    def local_print(level, *args):
        if verbose >= level:
            t = [f'{i}' for i in args]
            psi4.core.print_out(f'{" ".join(t)}\n')

    # If we do DFT, we the global reference should be either 'rks' or 'uks'.
    # Using 'rhf` or `uhf` for DFT calculation will lead to error in
    # psi4.core.HF constructor.
    if name.upper not in ['HF', 'SCF']:
        reference = psi4.core.get_global_option('REFERENCE')
        if reference == 'RHF':
            psi4.set_options({'reference': 'rks'})
        elif reference == 'UHF':
            psi4.set_options({'reference': 'uks'})

    mol = psi4.core.get_active_molecule()

    # ==> Sanity checks <==
    # TODO:
    # Currenty only support c1 symmetry.
    if mol.schoenflies_symbol() != 'c1':
        raise Exception('Current SCF code only supports C1 symmetry')

    if orbital_energy_unit not in ['eV', 'au']:
        raise Exception(
            'Invalid input for "orbital_energy_unit". Choices are ["eV", "au"].')

    if losc_ref_wfn:
        if not dfa_info:
            raise Exception("SCF-LOSC miss argument: dfa_info.")
        occ = losc_ref_wfn.losc_data['occ']
        C_lo = losc_ref_wfn.losc_data['C_lo']
        curvature = losc_ref_wfn.losc_data['curvature']

    if losc_ref_wfn:
        local_print(1, "---------------------------------------")
        local_print(1, "          SCF-LOSC (frozen-LO)")
        local_print(1, "             by Yuncai Mei")
        local_print(1, "---------------------------------------")
    else:
        local_print(1, "---------------------------------------")
        local_print(1, "                  SCF")
        local_print(1, "             by Yuncai Mei")
        local_print(1, "---------------------------------------")

    # Get global settings
    maxiter = psi4.core.get_global_option('MAXITER')
    E_conv = psi4.core.get_global_option('E_CONVERGENCE')
    D_conv = psi4.core.get_global_option('D_CONVERGENCE')
    reference = psi4.core.get_global_option('REFERENCE')
    basis = psi4.core.get_global_option('BASIS')
    local_print(1, '\n==> SCF settings <==')
    local_print(1, f'Reference: {reference}')
    local_print(1, f'Basis set: {basis}')
    local_print(1, f'MaxIter: {maxiter}')
    local_print(1, f'Energy convergence: {E_conv}')
    local_print(1, f'Density matrix convergence: {D_conv}\n')

    is_rks = True if reference == 'RKS' else False
    nspin = 1 if is_rks else 2

    # Build an SCF wavefunction with the corresponding DFT functional.
    optstash = psi4.driver.p4util.OptionsState(
        ['SCF', 'PRINT'])
    psi4.core.set_local_option('SCF', 'PRINT', 0)
    psi4.core.set_global_option('PRINT', 0)
    wfn = psi4.core.Wavefunction.build(mol, basis)
    wfn = psi4.proc.scf_wavefunction_factory(name, wfn, reference)

    supfunc = wfn.functional()
    nbf = wfn.basisset().nbf()

    # Create the occupation number, including the fractional cases.
    nocc, occ_idx, occ_val = utils.form_occ(wfn, occ)
    is_integer = is_integer(occ_val)
    is_aufbau = is_aufbau(nocc, occ_idx)
    local_print(1, "=> Occupation Number <=")
    local_print(1, f"Is integer system: {is_integer}")
    local_print(1, f"Is aufbau occupation: {is_aufbau}")
    for s in range(nspin):
        if not is_rks:
            local_print(
                1, f'{"Alpha" if s == 0 else "Beta"} Occupation Number:')
        for i in range(nocc[s]):
            local_print(1, 'spin={:<2d} idx={:<5d} occ={:<10f}'
                        .format(s, occ_idx[s][i], occ_val[s][i]))
        local_print(1, "")

    # ==> Set up matrices <==
    # S, H matrix
    if guess_wfn:
        S = np.asarray(guess_wfn.S())
        H = np.asarray(guess_wfn.H())
    else:
        S = np.asarray(wfn.mintshelper().ao_overlap())
        V = np.asarray(wfn.mintshelper().ao_potential())
        T = np.asarray(wfn.mintshelper().ao_kinetic())
        H = V + T
    # S^(-1/2)
    A = psi4.core.Matrix(nbf, nbf).from_array(S)
    A.power(-0.5, 1.e-16)
    A = np.asarray(A)

    # D, F, C, Cocc, eigenvalue matrices/vectors.
    D_psi = [psi4.core.Matrix(nbf, nbf) for s in range(nspin)]
    D = [np.asarray(m) for m in D_psi]
    F = [np.zeros((nbf, nbf)) for s in range(nspin)]
    C = [np.zeros((nbf, nbf)) for s in range(nspin)]
    Cocc_psi = [psi4.core.Matrix(nbf, nocc[s]) for s in range(nspin)]
    Cocc = [np.asarray(m) for m in Cocc_psi]
    eig = [np.zeros(nbf) for s in range(nspin)]
    if supfunc.needs_xc():
        Vpot = wfn.V_potential()
        Vpot.initialize()
        Vxc_psi = [psi4.core.Matrix(nbf, nbf) for s in range(nspin)]
        Vxc = [np.asarray(m) for m in Vxc_psi]

    # Set up DIIS helper object
    diis_helper = [diis.diis() for s in range(nspin)]
    diis_error = [np.zeros((nbf, nbf)) for s in range(nspin)]

    # Set up JK builder
    if is_integer and is_aufbau:
        local_print(1, "Use psi4.JK to calculate JK matrix.")
        jk_helper = jk.JK_psi4_jk(wfn, Cocc_psi)
    else:
        # To make fractional or non-aufbau calculation available, we need to
        # transform J and K matrices based on density matrix (explicitly
        # constructed based on occupation number). psi4 JK class does not
        # provide interface to accept density matrix as input to build J and K
        # matrices. So here we only the choice to transform from AO ERI to MO
        # ERI for J and K with help of psi4.MinsHelper library.
        local_print(1, "Use psi4.MintsHelper to calculate JK matrix.")
        jk_helper = jk.JK_psi4_mints(wfn, occ_val[:nspin], Cocc[:nspin])

    # ==> Set up initial guess <==
    # Build the initial C matrix as the initial guess.
    C_guess = []
    if guess_wfn:
        # SCF for DFA: guess CO from the input guess_wfn.
        local_print(1, "Set initial guess from a reference wavefunction.")
        C_guess = [np.asarray(guess_wfn.Ca()), np.asarray(guess_wfn.Cb())]
    elif losc_ref_wfn:
        # SCF for LOSC-DFA: initial guess has to tbe losc_ref_wfn.
        local_print(1, "Set initial guess from the parent DFA wavefunction.")
        C_guess = [np.asarray(losc_ref_wfn.Ca()),
                   np.asarray(losc_ref_wfn.Cb())]
    else:
        # SCF for DFA: set guess from psi4 guess setting.
        local_print(1, "Set initial guess from psi4 setting.")
        wfn.form_Shalf()
        wfn.guess()
        C_guess = [np.asarray(wfn.Ca()), np.asarray(wfn.Cb())]
    # Use the initial C guess to update Cocc, C, D matrix.
    for s in range(nspin):
        update_C(occ_idx[s], occ_val[s], C_guess[s], C[s], Cocc[s], D[s])

    # ==> SCF <==
    local_print(1, '\n=> SCF iterations <=')
    Eold = 0.0
    Enuc = wfn.molecule().nuclear_repulsion_energy()
    reference_factor = 2.0 if is_rks else 1.0
    Exc = 0.0
    energy_table = {}
    for SCF_ITER in range(1, maxiter + 1):
        # Build Vxc by using psi4.potential object. So we need to feed the
        # psi4.potential object with density matrix in psi4 matrix type.
        if supfunc.needs_xc():
            Vpot.set_D(D_psi)
            Vpot.compute_V(Vxc_psi)

        # Build J and K matrices.
        jk_helper.compute()
        J = jk_helper.J()
        K = jk_helper.K()
        J_tot = 2 * J[0] if nspin == 1 else J[0] + J[1]

        # Build Fock and diis helper object.
        G = []
        E_losc = 0
        for i in range(nspin):
            # build DFA Fock matrix
            G_tmp = J_tot - supfunc.x_alpha() * K[i]
            G.append(G_tmp)
            F[i][:] = H + G_tmp
            if supfunc.needs_xc():
                F[i][:] += Vxc[i]

            # ==> !!! The LOSC contribution in SCF: begin !!! <==
            if losc_ref_wfn:
                # build LOSC local occupation matrix
                local_occ = py_losc.local_occupation(C_lo[i], S, D[i])
                # build LOSC effective Fock matrix
                H_losc = py_losc.ao_hamiltonian_correction(
                    S, C_lo[i], curvature[i], local_occ)
                F[i][:] += H_losc
                # form LOSC energy correction
                E_losc += py_losc.energy_correction(curvature[i], local_occ)
            # ==> !!! The LOSC contribution in SCF: end !!! <==

            # DIIS error build and update
            diis_error[i] = F[i].dot(D[i]).dot(S) - S.dot(D[i]).dot(F[i])
            diis_error[i] = (A.T).dot(diis_error[i]).dot(A)
            diis_helper[i].add(F[i], diis_error[i])

        # Calculate SCF energy.
        if supfunc.needs_xc():
            Exc = Vpot.quadrature_values()['FUNCTIONAL']
        hf_E = 0
        for i in range(nspin):
            hf_E += reference_factor * \
                np.einsum('pq,pq->', D[i], H + 0.5 * G[i], optimize=True)
        SCF_E = Enuc + hf_E + Exc + E_losc
        energy_table['E_total'] = SCF_E
        energy_table['E_xc'] = Exc
        energy_table['E_losc'] = E_losc
        energy_table['E_nuc'] = Enuc
        energy_table['E_dfa'] = Enuc + hf_E + Exc

        # Print iteration steps
        dRMS = 0.5 * sum([np.mean(m ** 2) ** 0.5 for m in diis_error])
        local_print(1, '@%3s iter  %3d: Energy = %4.16f   dE = % 1.5E   dRMS = %1.5E  DIIS'
                    % (reference, SCF_ITER, SCF_E, (SCF_E - Eold), dRMS))

        # Check convergence
        if (abs(SCF_E - Eold) < E_conv) and (dRMS < D_conv):
            break

        # Extrapolate Fock matrices and update C, Cocc, D and eigenvalues.
        for s in range(nspin):
            # Form the extrapolated Fock matrix
            F_tmp = diis_helper[s].extrapolate()
            # Diagonalize Fock matrix and update C related variables including
            # wfn.epsilon, wfn.C, Cocc, wfn.D.
            H_tmp = A.dot(F_tmp).dot(A)
            eig[s][:], C2 = np.linalg.eigh(H_tmp)
            C_new = A.dot(C2)
            update_C(occ_idx[s], occ_val[s], C_new, C[s], Cocc[s], D[s])

        # Save scf energy
        Eold = SCF_E

        if SCF_ITER == maxiter:
            psi4.core.clean()
            raise Exception("Maximum number of SCF cycles exceeded.")

    # print total energies
    local_print(1, "\n=> Total Energy <=")
    local_print(1, "Nuclear potential energy: {:.10f}".format(energy_table['E_nuc']))
    local_print(1, "DFA total energy: {:.10f}".format(energy_table['E_dfa']))
    if supfunc.needs_xc():
        local_print(1, "Exchange-Correlation energy: {:.10f}".format(energy_table['E_xc']))
    if losc_ref_wfn:
        local_print(1, "LOSC correction energy: {:.10f}".format(energy_table['E_losc']))
    local_print(1, "SCF total energy: {:.10f}".format(SCF_E))

    # print orbital energies
    local_print(1, "\n=> Orbital Energy <=")
    occ_idx_val = [dict(zip(occ_idx[s], occ_val[s])) for s in range(nspin)]
    nbf = wfn.basisset().nbf()
    eig_factor = 1.0 if orbital_energy_unit == 'au' else constants.hartree2ev
    for s in range(nspin):
        if not is_rks:
            local_print(1, f"{'Alpha' if s == 0 else 'Beta'} orbital energy:")
        local_print(1, "{:<5s}  {:<8s}  {:>14s}".format(
            "Index", "Occ", f"DFA ({orbital_energy_unit})"))
        for i in range(nbf):
            local_print(1, "{:<5d}  {:<8.5f}  {:>14.6f}"
                        .format(i, occ_idx_val[s].get(i, 0),
                                eig[s][i] * eig_factor))
        local_print(1, "")

    # ==> update wfn <==
    # Update energy.
    wfn.set_energy(SCF_E)
    # Update matrices.
    # S: AO overlap matrix.
    # H: core matrix.
    # F: Fock matrix.
    # D: Density matrix.
    # Vxc: DFA Vxc matrix.
    # eig: orbital energy vector.
    wfn_S = np.asarray(wfn.S())
    wfn_H = np.asarray(wfn.H())
    wfn_F = [np.asarray(wfn.Fa()), np.asarray(wfn.Fb())]
    wfn_C = [np.asarray(wfn.Ca()), np.asarray(wfn.Cb())]
    wfn_D = [np.asarray(wfn.Da()), np.asarray(wfn.Db())]
    wfn_Vxc = [np.asarray(wfn.Va()), np.asarray(wfn.Vb())]
    wfn_eig = [np.asarray(wfn.epsilon_a()), np.asarray(wfn.epsilon_b())]

    wfn_S[:] = S
    wfn_H[:] = H
    for s in range(nspin):
        wfn_F[s][:] = F[s]
        wfn_D[s][:] = D[s]
        wfn_C[s][:] = C[s]
        wfn_eig[s][:] = eig[s]
        if supfunc.needs_xc():
            wfn_Vxc[s][:] = Vxc[s]

    # restore psi4 options.
    optstash.restore()
    return wfn


def scf(name, guess_wfn=None, occ={}, verbose=1, orbital_energy_unit='eV'):
    """
    Perform the SCF calculation for a normal DFA. It supports the systems with
    fractional occupations.

    Parameters
    ----------
    name: str
        The name of the psi4 superfunctional, including DFT functional or HF
        method. The same style as psi4.energy() function.
    guess_wfn: psi4.core.RHF or psi4.core.UHF, default to None.
        The wfn used to set the initial guess of the SCF calculation. Setting
        this variable will copy the coefficient matrix from `guess_wfn` as the
        initial SCF guess. Setting this variable will ignore the psi4 global
        guess setting. If you use this variable, make sure you are passing a
        reasonable guess wfn to the function. This function does not validate
        the input `guess_wfn`.
    occ: dict, default to an empty dict.
        A dictionary that specifies the customized occupation number. The
        occupation is obtained based on the aufbau occupation number of the
        current molecule to give the final set of occupation numbers for the SCF
        calculation. All the str in this dictionary is case-insensitive.

        The structure of the dict:
        {
            "alpha": # for alpha spin
            {
                int or str: float,
            }
            "beta": # for beta spin
            {
                int or str: float,
            }
        }
        The key-value pair (int: float): int refers to orbital index (0-based);
        float refers to the customized occupation number (in the range of [0, 1]).
        The key-value pair (str: float): str can be either ['homo', 'lumo'],
        which refers to HOMO and LUMO orbital respectively; float refers to the
        customized occupation number (in the range of [0, 1]).

        Example:
        >> H2 system: charge=0, mult=1
        The aufbau occupation numbers for H2 are:
        alpha occ: 1, 0, 0, 0, ...
        beta occ: 1, 0, 0, 0, ...
        >> Customized occ dictionary:
        occ = {
            "alpha": {"homo": 0.5}
            "beta": {"3": 0.7}
        }
        >> Resulted occupation numbers:
        alpha occ: 0.5, 0, 0, 0, ...
        beta occ:  1,   0, 0, 0.7, ...

    verbose: int, default to 1
        print level to the psi4 output file.
        0 means print nothing. 1 means normal print level. A larger number means
        more details.
    orbital_energy_unit: str, default to 'eV'
        The units of orbital energies used to print in the output.
        Valid choices are ['au', 'eV'].
        'au': atomic unit, hartree.
        'eV': electronvolt.

    Returns
    -------
    wfn: psi4.core.RHF or psi4.core.UHF
        The SCF wavefunction.
        1. Following members in the returned wfn object are updated:
            wfn.S(): AO overlap matrix.
            wfn.H(): Core matrix.
            wfn.Ca(), wfn.Cb(): CO coefficient matrix.
            wfn.Fa(), wfn.Fb(): Fock matrix.
            wfn.Da(), wfn.Db(): density matrix.
            wfn.Va(), wfn.Vb(): DFA (just the DFA, if it is LOSC-DFA) Vxc matrix.
            wfn.epsilon_a(), wfn.epsilon_b(): orbital energies.
            wfn.energy(): total energy.

        2. Other members are not touched. Accessing and using other members
        in the returned wfn may be a undefined behavior.

        3. A new attributes will be added into the returned wfn object:
            losc_data: dict
                The data that can be used for post-SCF-LOSC and SCF-LOSC
                calculations.
                losc_data['occ']: dict
                    The occupation information that is copied from the input
                    `occ`. See `psi4_losc.scf()` for detailed description of
                    occupation information.

    Notes
    -----
    1. This function only relies on the GLOBAL settings in psi4, NOT any local
    settings for any psi4 modules!

    2. There are double memory cost for all the updated matrices/vectors in
    the returned wfn (such as C, D, F matrices): one is for the allocation in
    the python code, the other one is for the allocation in psi4.core C++ code.
    At return, matrices/vectors will be copied from the python side to the psi4
    C++ side code. This is limited by the psi4.core interface. We have to pay
    the price.
    """
    wfn = _scf(name, guess_wfn=guess_wfn, occ=occ, verbose=verbose,
               orbital_energy_unit=orbital_energy_unit)

    # add new attributes
    wfn.losc_data = {'occ': occ.copy()}
    return wfn


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


def scf_losc(dfa_info, dfa_wfn, orbital_energy_unit='eV', verbose=1):
    """
    Perform the SCF-LOSC (frozen-LO) calculation based on a DFA wavefunction.

    This function supports SCF-LOSC (frozen-LO) calculations for integer/fractional
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
    dfa_wfn: psi4.core.HF like psi4 object
        The converged wavefunction from a parent DFA.
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
    wfn: psi4.core.RHF or psi4.core.UHF
        The wavefunction for the SCF-LOSC (frozen-LO) calculation.
        1. Following members in the returned wfn object are calculated:
            wfn.S(): AO overlap matrix.
            wfn.H(): Core matrix.
            wfn.Ca(), wfn.Cb(): CO coefficient matrix.
            wfn.Fa(), wfn.Fb(): Fock matrix.
            wfn.Da(), wfn.Db(): density matrix.
            wfn.Va(), wfn.Vb(): DFA (just the DFA, if it is LOSC-DFA) Vxc matrix.
            wfn.epsilon_a(), wfn.epsilon_b(): orbital energies.
            wfn.energy(): total energy.

        2. Accessing and using other members in the returned wfn is a undefined
        behavior.

        3. A new attributes will be added into the returned wfn object:
            losc_data: dict
                The data that can be used for post-SCF-LOSC and SCF-LOSC
                calculations.
                losc_data['occ']: dict
                    The occupation information. See `psi4_losc.scf()` for more
                    detailed description of occupation information.

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
    1. This function only relies on the GLOBAL settings in psi4, NOT any local
    settings for any psi4 modules!

    2. There are double memory cost for all the updated matrices/vectors in
    the returned wfn (such as C, D, F matrices): one is for the allocation in
    the python code, the other one is for the allocation in psi4.core C++ code.
    At return, matrices/vectors will be copied from the python side to the psi4
    C++ side code. This is limited by the psi4 core interface. We have to pay
    the price.
    """
    post_scf_losc(dfa_info, dfa_wfn,
                  orbital_energy_unit=orbital_energy_unit, verbose=verbose)
    dfa_name = dfa_wfn.functional().name()
    wfn = _scf(dfa_name, losc_ref_wfn=dfa_wfn, dfa_info=dfa_info,
               orbital_energy_unit=orbital_energy_unit, verbose=verbose)
    wfn.losc_data = {'occ': dfa_wfn.losc_data['occ'].copy()}

    return wfn
