"""
Extended SCF procedure for SCF-DFA and SCF-LOSC-DFA (frozen-LO) with
compability to fractional systems.

Notes
-----
If you are interested in calculations of integer systems for LOSC,
use `psi4_losc.post_scf_losc()` or `psi4_losc.scf_losc()` instead, since their
SCF procedure are handled by psi4 and may be more efficient and reliable.
"""

import psi4
import numpy as np
import py_losc
from psi4_losc import post_scf_losc
from psi4_losc import utils
from psi4_losc import diis
from psi4_losc import jk
from qcelemental import constants


def _scf(name, guess_wfn=None, losc_ref_wfn=None, curvature=None, C_lo=None,
         dfa_info=None, occ={}, verbose=1, orbital_energy_unit='eV'):
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
        guess setting; (2) if you use this variable, make sure you are
        passing a reasonable guess wfn to the function. This function does not
        check the validity the input.
    losc_ref_wfn: psi4.core.RHF or psi4.core.UHF, default to None.
        The wfn used to do SCF-LOSC (frozen-LO) calculation. This variable is
        the flag to trigger SCF-LOSC (frozen-LO) calculation.
    curvature: [np.array, ...]
        The LOSC curvature matrix.
    C_lo: [np.array, ...]
        The LOSC LO coefficient matrix.
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
    1. The option settins are handled by psi4 SCF module.

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

    optstash = psi4.driver.p4util.OptionsState(
        ['SCF', 'REFERENCE'],
        ['SCF', 'PRINT'])

    local_print = utils.init_local_print(verbose)

    # If we do DFT, we the reference should be either 'rks' or 'uks'.
    # Using 'rhf` or `uhf` for DFT calculation will lead to error in
    # psi4.core.HF constructor.
    is_hf = True
    if name.upper() not in ['HF', 'SCF']:
        is_hf = False
        reference = psi4.core.get_option('SCF', 'REFERENCE')
        if reference == 'RHF':
            psi4.core.set_local_option('SCF', 'REFERENCE', 'RKS')
        elif reference == 'UHF':
            psi4.core.set_local_option('SCF', 'REFERENCE', 'UKS')

    mol = psi4.core.get_active_molecule()
    # print out molecule information
    if verbose >= 1:
        mol.print_out()

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
        if not curvature:
            raise Exception("SCF-LOSC miss argument: curvature.")
        if not C_lo:
            raise Exception("SCF-LOSC miss argument: C_lo.")

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

    # Get psi4 settings
    maxiter = psi4.core.get_option('SCF', 'MAXITER')
    E_conv = psi4.core.get_option('SCF', 'E_CONVERGENCE')
    D_conv = psi4.core.get_option('SCF', 'D_CONVERGENCE')
    reference = psi4.core.get_option('SCF', 'REFERENCE')
    is_diis_rms = psi4.core.get_option('SCF', 'DIIS_RMS_ERROR')
    is_dfjk = psi4.core.get_option('SCF', 'SCF_TYPE').endswith('DF')
    basis = psi4.core.get_option('SCF', 'BASIS')
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
    wfn = psi4.core.Wavefunction.build(mol, basis)
    wfn = psi4.driver.scf_wavefunction_factory(name, wfn, reference)
    mintshelper = psi4.core.MintsHelper(wfn.basisset())

    supfunc = wfn.functional()
    nbf = wfn.basisset().nbf()

    # Create the occupation number, including the fractional cases.
    nocc, occ_idx, occ_val = utils.form_occ(wfn, occ)
    nelec = [sum(x) for x in occ_val]
    is_integer = utils.is_integer_system(wfn, occ)
    is_aufbau = utils.is_aufbau_system(wfn, occ)
    local_print(1, "=> Occupation Number <=")
    local_print(1, f"Is integer system: {is_integer}")
    local_print(1, f"Is aufbau occupation: {is_aufbau}")
    local_print(
        1, f"nelec_a: {nelec[0]} {f' ! update and differ to wfn.nalpha(): nalpha={wfn.nalpha()}.' if nelec[0] != wfn.nalpha() else ''}")
    local_print(
        1, f"nelec_b: {nelec[1]} {f' ! update and differ to wfn.nbeta(): nbeta={wfn.nbeta()}.' if nelec[1] != wfn.nbeta() else ''}")
    local_print(1, "")
    for s in range(nspin):
        if not is_rks:
            local_print(
                1, f'{"Alpha" if s == 0 else "Beta"} Occupation Number:')
        else:
            local_print(1, 'Occupation Number:')
        max_idx = -1 if not occ_idx[s] else occ_idx[s][-1]
        occ_idx_val = dict(zip(occ_idx[s], occ_val[s]))
        for i in range(max_idx+1):
            occ_val_t = occ_idx_val.get(i, 0)
            local_print(1, 'index={:<5d} occ={:<10f}'
                        .format(i, occ_val_t))
        local_print(1, "")

    # ==> Set up matrices <==
    # S, H matrix
    if guess_wfn:
        S = np.asarray(guess_wfn.S())
        H = np.asarray(guess_wfn.H())
    else:
        S = np.asarray(mintshelper.ao_overlap())
        V = np.asarray(mintshelper.ao_potential())
        if losc_ref_wfn:
            T = np.asarray(mintshelper.ao_kinetic())
            H = V + T
        else:
            mintshelper.one_electron_integrals()
            wfn.form_H()
            H = np.array(wfn.H())
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
    local_print(1, "%s                        Total Energy        Delta E     %s |[F,P]|\n" %
                ("   " if is_dfjk else "", "RMS" if is_diis_rms else "MAX"))

    Eold = 0.0
    Enuc = wfn.molecule().nuclear_repulsion_energy()
    reference_factor = 2.0 if is_rks else 1.0
    Exc = 0.0
    energy_table = {}
    scf_status = ['DIIS']
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
        E_J = 0
        E_K = 0
        E_H = 0
        for i in range(nspin):
            E_J += reference_factor * \
                np.einsum('pq,pq->', D[i], 0.5 * J_tot, optimize=True)
            E_K += reference_factor * \
                np.einsum('pq,pq->', D[i], -0.5 *
                          supfunc.x_alpha() * K[i], optimize=True)
            E_H += reference_factor * \
                np.einsum('pq,pq->', D[i], H, optimize=True)
        SCF_E = Enuc + E_H + E_J + E_K + Exc + E_losc
        energy_table['E_total'] = SCF_E
        energy_table['E_H'] = E_H
        energy_table['E_J'] = E_J
        energy_table['E_K'] = E_K
        energy_table['E_xc'] = Exc
        energy_table['E_G'] = E_J + E_K
        energy_table['E_hf'] = E_H + E_J + E_K
        energy_table['E_losc'] = E_losc
        energy_table['E_nuc'] = Enuc
        energy_table['E_dfa'] = Enuc + E_H + E_J + E_K + Exc

        # Print iteration steps
        D_error = 0
        if is_diis_rms:  # rms type error
            rms = [np.sqrt(np.mean(np.square(m))) for m in diis_error]
            D_error = np.sqrt(np.mean(np.square(rms)))
        else:  # max type error
            D_error = max([np.max(np.abs(m)) for m in diis_error])
        local_print(1,
                    "   @%s%s iter %3s: %20.14f   %12.5e   %-11.5e %s" %
                    ("DF-" if is_dfjk else "", reference, SCF_ITER,
                     SCF_E, (SCF_E - Eold), D_error, '/'.join(scf_status)))

        # Check convergence
        if (abs(SCF_E - Eold) < E_conv) and (D_error < D_conv):
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
    func_type = "HF" if is_hf else "DFA"
    local_print(1, "\n=> Total Energy <=")
    local_print(1, "Nuclear potential energy: {:.10f}".format(
        energy_table['E_nuc']))
    local_print(1, "{:s} core energy: {:.10f}".format(
        func_type, energy_table['E_H']))
    local_print(1, "{:s} J energy: {:.10f}".format(
        func_type, energy_table['E_J']))
    local_print(1, "{:s} K energy: {:.10f}".format(
        func_type, energy_table['E_K']))
    local_print(1, "{:s} G energy: {:.10f}".format(
        func_type, energy_table['E_G']))
    if supfunc.needs_xc():
        local_print(1, "{:s} XC energy: {:.10f}".format(
            func_type, energy_table['E_xc']))
    local_print(1, "{:s} total energy: {:.10f}".format(
        func_type, energy_table['E_dfa']))
    if losc_ref_wfn:
        local_print(1, "LOSC correction energy: {:.10f}".format(
            energy_table['E_losc']))
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
        initial SCF guess. Setting this variable will ignore the psi4 guess
        setting. If you use this variable, make sure you are passing a
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

        3. The input argument `occ` is added as a new attribute to the returned
        wfn object as `wfn.losc_data['occ'] = occ`.

    Notes
    -----
    1. The option settings are handeld by psi4 SCF module.

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


def scf_losc(dfa_info, dfa_wfn, orbital_energy_unit='eV', verbose=1):
    """
    Perform the SCF-LOSC (frozen-LO) calculation based on a DFA wavefunction.

    This function supports SCF-LOSC (frozen-LO) calculations for integer/fractional
    systems with aufbau/non-aufbau occupations:
    (1) If you want to calculate integer system with aufbau occupations, use
    `psi4.energy()` or `psi4_losc.scf.scf()` to generate the input `dfa_wfn`.
    (2) If you want to calculate integer system with non-aufbau occupation, or
    fractional system with aufbau/non-aufbau occupations, use `psi4_losc.scf.scf()`
    to generate the input `dfa_wfn`.

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

    See Also
    --------
    py_losc.DFAInfo(): constructor of the DFA information class.
    psi4.energy(): return a DFA SCF wavefunction. psi4.energy() only supports
        calculations for integer systems with aufbau occupations.
    psi4_losc.scf.scf(): return a DFA SCF wavefunction. `psi4_losc.scf.scf()`
        supports calculations for integer/fractional systems with aufbau/non-aufbau
        occupations.

    Notes
    -----
    1. The option settings are handled by psi4 SCF module.

    2. There are double memory cost for all the updated matrices/vectors in
    the returned wfn (such as C, D, F matrices): one is for the allocation in
    the python code, the other one is for the allocation in psi4.core C++ code.
    At return, matrices/vectors will be copied from the python side to the psi4
    C++ side code. This is limited by the psi4 core interface. We have to pay
    the price.
    """
    _, _, losc_data = post_scf_losc(dfa_info, dfa_wfn, verbose=verbose,
                                    orbital_energy_unit=orbital_energy_unit,
                                    return_losc_data=True)
    dfa_name = dfa_wfn.functional().name()
    occ = {}
    if hasattr(dfa_wfn, 'losc_data'):
        occ = dfa_wfn.losc_data['occ']
    wfn = _scf(dfa_name, occ=occ, losc_ref_wfn=dfa_wfn, dfa_info=dfa_info,
               curvature=losc_data['curvature'], C_lo=losc_data['C_lo'],
               orbital_energy_unit=orbital_energy_unit, verbose=verbose)

    return wfn
