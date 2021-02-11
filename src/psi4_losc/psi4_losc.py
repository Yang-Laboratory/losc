"""
Integrate py_losc module with psi4 package to perform LOSC calculation.
"""

import numpy as np
from psi4_losc import jk
from psi4_losc import diis

import psi4
from pkg_resources import parse_version
from psi4.driver.p4util.exceptions import ValidationError
from py_losc import py_losc
from qcelemental import constants

if parse_version(psi4.__version__) >= parse_version('1.3a1'):
    build_superfunctional = psi4.driver.dft.build_superfunctional
else:
    build_superfunctional = psi4.driver.dft.build_superfunctional


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


def _print_header(scf_losc=True):
    if scf_losc:
        name = "SCF-LOSC (frozen-LO)"
    else:
        name = "post-SCF-LOSC"
    psi4.core.print_out('---------------------------------------\n')
    psi4.core.print_out(f'          {name}\n')
    psi4.core.print_out(f'          By Yuncai Mei\n')
    psi4.core.print_out('---------------------------------------\n')


def form_df_basis_matrix(wfn):
    """
    Form the three-center integral used in density fitting.

    Parameters
    ----------
    wfn: psi4.core.HF
        A psi4 wavefunction object that is associated with a DFT calculations.

    Returns
    -------
    df_pmn: np.array [nfitbasis, nbasis x nbasis]
        The three-center integral of <fitbasis|basis, basis>.
    df_Vpq_inv: np.array [nfitbasis, nfitbasis]
        The inverse of matrix <fitbasis|1/r|fitbasis>.
    """
    # basis set
    basis = wfn.basisset()
    zero_bas = psi4.core.BasisSet.zero_ao_basis_set()
    # auxiliary basis for density fitting. TODO: let user choose the fitbasis.
    aux_bas = psi4.core.BasisSet.build(
        wfn.molecule(), "ORBITAL", "my_aug-cc-pVTZ")
    aux_bas.print_out()
    # psi4.mintshelper object to help building AO integrals.
    mints = wfn.mintshelper()
    # build three-center integral <fitbasis|lo, lo>
    df_pmn = np.asarray(mints.ao_eri(aux_bas, zero_bas, basis, basis))
    df_pmn = np.squeeze(df_pmn)
    # build density fitting Vpq inverse
    df_Vpq = np.asarray(mints.ao_eri(aux_bas, zero_bas, aux_bas, zero_bas))
    df_Vpq = np.squeeze(df_Vpq)
    df_Vpq_inv = np.linalg.inv(df_Vpq)

    return df_pmn, df_Vpq_inv


def form_grid_lo(wfn, C_lo):
    """
    Form the matrix that stores the values of LOs on grid points.

    Parameters
    ----------
    wfn: psi4.core.HF
        A psi4 wavefunction object that is associated with a DFT calculations.
    C_lo: np.array [nbasis, nlo]
        LO coefficient matrix on AOs.

    Returns
    -------
    grid_lo: np.array [npts, nlo]
        A matrix for the values of LOs on grid points. npts is the number of
        grid points, and nlo is the number of LOs.
    """
    # psi4.VBase object to help building grid.
    Vpot = wfn.V_potential()
    # number of grid points
    npts = Vpot.grid().npoints()
    # number of LOs.
    nlo = C_lo.shape[1]
    # Get the psi4.potential object for help
    # Note:
    # Vpot.properties() returns `nthread` number of pointer functions
    # (psi4 core object). These multiple pointer functions work in parallel
    # in psi4 core.
    # Here, we work in the python layer and do not use parallel construction
    # of grid yet. So we always use the first pointer function object to
    # build grid.
    points_func = Vpot.properties()[0]

    # allocate grid_lo matrices.
    grid_lo = np.zeros((npts, nlo))

    # loop over the blocks to build grid_lo
    npts_count = 0
    for b in range(Vpot.nblocks()):
        # Obtain block information
        block = Vpot.get_block(b)
        points_func.compute_points(block)
        npoints = block.npoints()
        lpos = np.array(block.functions_local_to_global())

        # Compute grid_ao on the fly to build grid_lo.
        # The full matrix of grid_ao is [npts, nbasis]. Since the grid is very
        # sparse on atoms, we do not always need the full matrix. `lpos`
        # represents the AOs that have non-zero grid values for the given grid
        # block.
        grid_ao = np.array(points_func.basis_values()["PHI"])[:npoints, :lpos.shape[0]]
        grid_lo_blk = grid_lo[npts_count:npts_count+npoints, :]
        grid_lo_blk[:] = grid_ao.dot(C_lo[lpos, :])  # npoints x nlo

        # update block starting position
        npts_count += npoints

    return grid_lo


def form_grid_w(wfn):
    """
    The grid points weights

    Parameters
    ----------
    wfn: psi4.core.HF
        A psi4 wavefunction object that is associated with a DFT calculations.

    Returns
    -------
    grid_w: np.array [npts]
        A vector for the weights of grid points.


    References
    ----------
    psi4numpy/Tutorials/04_Density_Functional_Theory/4b_LDA_kernel.ipynb
    """
    # psi4.VBase object to help building grid.
    Vpot = wfn.V_potential()
    # number of grid points
    npts = Vpot.grid().npoints()

    # build grid_w vector
    grid_w = np.zeros((npts,))
    npts_count = 0
    # loop over the blocks to build grid_w
    for b in range(Vpot.nblocks()):
        # Obtain block of weights
        block = Vpot.get_block(b)
        npoints = block.npoints()
        grid_w[npts_count:npts_count+npoints] = np.asarray(block.w())

        # update block starting position
        npts_count += npoints

    return grid_w


def post_scf_losc(dfa_wfn, orbital_energy_unit='eV'):
    """
    Perform the post-SCF-LOSC calculation.

    Parameters
    ----------
    dfa_wfn: psi4.core.HF like psi4 object
        The converged wavefunction from a parent DFA.
    orbital_energy_unit: str, default to 'eV'
        The units of orbital energies used to print in the output.
        Valid choices are ['au', 'eV'].
        'au': atomic unit, hartree.
        'eV': electronvolt.

    Returns
    -------
    energy: float
        The total energy of LOSC-DFA.
    eig: [list, ...]
        The orbital energies from LOSC-DFA. `len(eig)` equals 1 if it is RKS,
        and equals 2 if it is UKS. `len(eig[0])` equals nlo.
    """
    # sanity-check of input dfa wfn.
    _validate_dfa_wfn(dfa_wfn)
    if orbital_energy_unit not in ['eV', 'au']:
        raise Exception(
            'Invalid input for "orbital_energy_unit". Choices are ["eV", "au"].')
    eig_factor = 1.0 if orbital_energy_unit == 'au' else constants.hartree2ev

    _print_header(scf_losc=False)
    is_rks = dfa_wfn.same_a_b_orbs() and dfa_wfn.same_a_b_dens()
    nspin = 1 if is_rks else 2

    # map needed matrices to DFA wfn.
    C_co = [np.asarray(dfa_wfn.Ca()), np.asarray(dfa_wfn.Cb())]
    H_ao = [np.asarray(dfa_wfn.Fa()), np.asarray(dfa_wfn.Fb())]
    D_ao = [np.asarray(m) for m in dfa_wfn.mintshelper().ao_dipole()]
    S = np.asarray(dfa_wfn.S())
    D = [np.asarray(dfa_wfn.Da()), np.asarray(dfa_wfn.Db())]

    # ====> !!! Start to construct LOSC !!! <====
    # ==> LOSC localization <==
    C_lo = [None] * nspin
    for s in range(nspin):
        # create losc localizer object
        localizer = py_losc.LocalizerV2(C_co[s], H_ao[s], D_ao)
        # compute LOs
        C_lo[s] = localizer.lo()

    # ==> LOSC curvature matrix <==
    # build DFA information object for py_losc module
    # TODO: how to get following dfa information from psi4.superfunction?
    gga_x = 0.8
    hf_x = 0.2
    dfa_info = py_losc.DFAInfo(gga_x, hf_x, name='b3lyp')

    # build matrices related to density fitting
    df_pmn, df_Vpq_inv = form_df_basis_matrix(dfa_wfn)
    df_pii = [None] * nspin
    for s in range(nspin):
        # build three-center integral <fitbasis|lo, lo>
        df_pii[s] = np.einsum('pmn,mi,ni->pi', df_pmn,
                              C_lo[s], C_lo[s], optimize=True)

    # build weights of grid points
    grid_w = form_grid_w(dfa_wfn)

    # build values of LOs on grid points
    grid_lo = [form_grid_lo(dfa_wfn, C_lo_) for C_lo_ in C_lo]

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

    # ==> Print energies to output <==
    psi4.core.print_out('\n\n=> Total Energy <=\n')
    psi4.core.print_out('DFA SCF energy: {:.8f} hartree\n'.format(dfa_wfn.energy()))
    psi4.core.print_out('LOSC energy: {:.8f} hartree\n'.format(E_losc_tot))
    psi4.core.print_out('LOSC-DFA energy: {:.8f} hartree\n'.format(dfa_wfn.energy() + E_losc_tot))

    # print orbital energies
    psi4.core.print_out('\n\n=> Orbital Energy <=\n')
    dfa_eigs = [np.asarray(dfa_wfn.epsilon_a()) * eig_factor,
                np.asarray(dfa_wfn.epsilon_b()) * eig_factor]
    for spin in range(nspin):
        psi4.core.print_out(
            f'\n{"Alpha" if spin == 0 else "Beta"} orbital energies in {orbital_energy_unit}\n')
        header = "{:<6s} {:>20s} {:>20s}\n".format(
            "Index", "DFA", "LOSC-DFA")
        psi4.core.print_out(header)
        for i in range(len(dfa_eigs[spin])):
            psi4.core.print_out('{:<6d} {:>20.8f} {:>20.8f}\n'.format(
                i, dfa_eigs[spin][i], losc_eig[spin][i]))

    return dfa_wfn.energy() + E_losc_tot, losc_eig


def scf_losc(dfa_wfn, inplace=True, orbital_energy_unit='eV'):
    """
    Perform the SCF-LOSC (frozen-LO) calculation.

    Parameters
    ----------
    dfa_wfn: psi4.core.HF like psi4 object
        The converged wavefunction from a parent DFA. It will be modified if the
        `inplace` is set to True.
    inplace: bool, default to True.
        Modify the `dfa_wfn` in place.
    orbital_energy_unit: str, default to 'eV'
        The units of orbital energies used to print in the output.
        Valid choices are ['au', 'eV'].
        'au': atomic unit, hartree.
        'eV': electronvolt.

    Returns
    -------
    wfn: psi4.core.HF
        A wavefunction that represents LOSC-DFA functional. The type of `wfn`
        is the same `dfa_wfn`. So it would be either `psi4.core.RHF` or
        `psi4.core.UHF`. You can the same interface as defined in psi4.

        The following internal data of `wfn` are modifed to involve LOSC
        contributions:
        1. Fock matrix: modified by `wfn.Fa()` and `wfn.Fb()` interface.
        2. Density matrix: modified by `wfn.Da()` and `wfn.Db()` interface.
        3. CO coefficient matrix: modified by `wfn.Ca()` and `wfn.Cb()` interface.
        4. CO energy: modified by `wfn.epsilon_a()` and `wfn.epsilon_b()` interface.
        5. total energy: modified by `wfn.set_energy()` interface.

        The following new data will be added into `wfn` object:
        `wnf.losc_E`: float
            the LOSC energy correction.
    """
    # ==> Sanity-check of input dfa_wfn <==
    _validate_dfa_wfn(dfa_wfn)
    if orbital_energy_unit not in ['eV', 'au']:
        raise Exception('Invalid input for "orbital_energy_unit". Choices are ["eV", "au"].')
    eig_factor = 1.0 if orbital_energy_unit == 'au' else constants.hartree2ev

    _print_header(scf_losc=True)
    if inplace:
        wfn = dfa_wfn
    else:
        wfn = dfa_wfn.deep_copy()

    # ==> get settings and basic information <==
    maxiter = psi4.core.get_option('SCF', 'MAXITER')
    E_conv = psi4.core.get_option('SCF', 'E_CONVERGENCE')
    D_conv = psi4.core.get_option('SCF', 'D_CONVERGENCE')
    is_rks = wfn.same_a_b_orbs() and wfn.same_ab_dens()
    nspin = 1 if is_rks else 2
    nbf = wfn.nso()
    nelec = [wfn.nalpha(), wfn.nbeta()]

    psi4.core.print_out('\n=> Settings and basic information <=\n')
    psi4.core.print_out(f'max_iter: {maxiter}\n')
    psi4.core.print_out(f'E_conv: {E_conv}\n')
    psi4.core.print_out(f'D_conv: {D_conv}\n')
    psi4.core.print_out(f'Reference: {"RKS" if is_rks else "UKS"}\n')
    psi4.core.print_out(f'nspin: {nspin}\n')
    psi4.core.print_out(
        f'Number of electrons: alpha={nelec[0]} beta={nelec[1]}\n')
    psi4.core.print_out(
        f'Number of singly occupied orbitals: {nelec[0] - nelec[1]}\n')
    psi4.core.print_out(f'Number of basis functions: {nbf}\n')

    # ==> Set up matrices used in scf <==
    S = np.asarray(wfn.S())  # AO overlap matrix
    H = np.asarray(wfn.H())  # core matrix.
    A = psi4.core.Matrix(nbf, nbf).from_array(S)  # S^(-1/2)
    A.power(-0.5, 1.e-16)
    A = np.asarray(A)
    Cocc_psi = [wfn.Ca_subset('SO', 'OCC'), wfn.Cb_subset('SO', 'OCC')]
    Cocc = [np.asarray(m) for m in Cocc_psi]
    C = [np.asarray(wfn.Ca()), np.asarray(wfn.Cb())]
    D = [np.asarray(wfn.Da()), np.asarray(wfn.Db())]
    F = [np.asarray(wfn.Fa()), np.asarray(wfn.Fb())]
    eig = [np.asarray(wfn.epsilon_a()), np.asarray(wfn.epsilon_b())]
    C_co = C
    H_ao = F
    D_ao = [np.asarray(m) for m in dfa_wfn.mintshelper().ao_dipole()]

    # ==> Set up psi4.JK, Vxc and DIIS helper used in scf <==
    # psi4 JK object to build J and K matrices.
    jk = psi4.core.JK.build(wfn.basisset())
    jk.initialize()
    for i in range(nspin):
        jk.C_left_add(Cocc_psi[i])

    # psi4 DFA potential object to build Vxc matrix
    supfunc = dfa_wfn.functional()
    if supfunc.needs_xc():
        Vpot = wfn.V_potential()
        Vpot.initialize()
        Vxc = [np.asarray(wfn.Va()), np.asarray(wfn.Vb())]

    # DIIS helper object
    diis_helper = [diis.diis() for _ in range(nspin)]
    diis_error = [np.zeros((nbf, nbf)) for _ in range(nspin)]

    # ==> !!! LOSC: build Frozen-LO and LOSC curvature !!! <==
    # => step1: build frozen LO <=
    C_lo = [None] * nspin
    for s in range(nspin):
        # create losc localizer object
        localizer = py_losc.LocalizerV2(C[s], F[s], D_ao)
        # compute LOs
        C_lo[s] = localizer.lo()

    # => step2: build frozen curvature matrix <==
    # build DFA information object for py_losc module
    # TODO: how to get following dfa information from psi4.superfunction?
    gga_x = 0.8
    hf_x = 0.2
    dfa_info = py_losc.DFAInfo(gga_x, hf_x, name='b3lyp')

    # build matrices related to density fitting
    df_pmn, df_Vpq_inv = form_df_basis_matrix(dfa_wfn)
    df_pii = [None] * nspin
    for s in range(nspin):
        # build three-center integral <fitbasis|lo, lo>
        df_pii[s] = np.einsum('pmn,mi,ni->pi', df_pmn, C_lo[s], C_lo[s], optimize=True)

    # build weights of grid points
    grid_w = form_grid_w(dfa_wfn)

    # build values of LOs on grid points
    grid_lo = [form_grid_lo(dfa_wfn, C_lo_) for C_lo_ in C_lo]

    # build LOSC curvature matrices
    curvature = [None] * nspin
    for s in range(nspin):
        # build losc curvature matrix
        curvature_helper = py_losc.CurvatureV2(
            dfa_info, df_pii[s], df_Vpq_inv, grid_lo[s], grid_w)
        curvature[s] = curvature_helper.kappa()

    # ==> SCF-LOSC (frozen-LO) <==
    psi4.core.print_out('\n\n=> Start SCF iterations <=\n\n')
    Eold = 0.0
    Enuc = wfn.molecule().nuclear_repulsion_energy()
    x_alpha = 1.0 if not supfunc.needs_xc() else Vpot.functional().x_alpha()
    reference_factor = 2.0 if is_rks else 1.0
    Exc = 0.0
    for SCF_ITER in range(1, maxiter + 1):
        # build Vxc using psi4 potential object
        if supfunc.needs_xc():
            Vpot.set_D([wfn.Da(), wfn.Db()])
            Vpot.compute_V([wfn.Va(), wfn.Vb()])

        # build JK using psi4.libJK
        jk.compute()
        J = [np.asarray(jk.J()[i]) for i in range(nspin)]
        K = [np.asarray(jk.K()[i]) for i in range(nspin)]
        J_tot = 2 * J[0] if nspin == 1 else J[0] + J[1]

        # build Fock and diis helper object.
        G = []
        E_losc = 0
        for i in range(nspin):
            # build DFA Fock matrix
            G_tmp = J_tot - x_alpha * K[i]
            G.append(G_tmp)
            F[i][:] = H + G_tmp
            if supfunc.needs_xc():
                F[i][:] += Vxc[i]

            # ==> !!! The LOSC contribution in SCF: begin !!! <==
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

        # SCF energy and update: [Szabo:1996], exercise 3.40, pp. 215
        if supfunc.needs_xc():
            Exc = Vpot.quadrature_values()['FUNCTIONAL']
        SCF_E = 0
        for i in range(nspin):
            SCF_E += reference_factor * \
                np.einsum('pq,pq->', D[i], H + 0.5 * G[i])
        SCF_E += Enuc + Exc + E_losc

        # print iteration steps
        dRMS = 0.5 * sum([np.mean(m ** 2) ** 0.5 for m in diis_error])
        psi4.core.print_out('SCF Iteration %3d: Energy = %4.16f   dE = % 1.5E   dRMS = %1.5E\n'
                            % (SCF_ITER, SCF_E, (SCF_E - Eold), dRMS))

        # check convergence
        if (abs(SCF_E - Eold) < E_conv) and (dRMS < D_conv):
            break

        # Extrapolate Fock matrices and update density
        for i in range(nspin):
            # Form the extrapolated Fock matrix
            F_tmp = diis_helper[i].extrapolate()
            # Diagonalize Fock matrix and update C related variables including
            # wfn.epsilon, wfn.C, Cocc, wfn.D.
            H_tmp = A.dot(F_tmp).dot(A)
            eig[i][:], C2 = np.linalg.eigh(H_tmp)
            C[i][:] = A.dot(C2)
            Cocc[i][:] = C[i][:, :nelec[i]]
            D[i][:] = np.einsum('pi,qi->pq', Cocc[i], Cocc[i])

        # save scf energy
        Eold = SCF_E

        if SCF_ITER == maxiter:
            psi4.core.clean()
            raise Exception("Maximum number of SCF cycles exceeded.")

    wfn.set_energy(SCF_E)
    psi4.core.print_out('\n\n=> Total Energy <=\n')
    psi4.core.print_out('Final SCF energy: {:.8f} hartree\n'.format(SCF_E))

    # print orbital energies
    psi4.core.print_out('\n\n=> Orbital Energy <=\n')
    dfa_eigs = [np.asarray(dfa_wfn.epsilon_a()) * eig_factor,
                np.asarray(dfa_wfn.epsilon_b()) * eig_factor]
    losc_eigs = [i * eig_factor for i in eig]
    for spin in range(nspin):
        psi4.core.print_out(
            f'\n{"Alpha" if spin == 0 else "Beta"} orbital energies in {orbital_energy_unit}.\n')
        header = "{:<6s} {:>20s} {:>20s}\n".format(
            "Index", "DFA", "LOSC-DFA")
        psi4.core.print_out(header)
        for i in range(len(dfa_eigs[spin])):
            psi4.core.print_out('{:<6d} {:>20.8f} {:>20.8f}\n'.format(
                i, dfa_eigs[spin][i], losc_eigs[spin][i]))
    return wfn


def _scf(name, guess_wfn=None, losc_ref_wfn=None, occ={}, verbose=1):
    """
    Perform the SCF calculation for systems that may have fractional occupations.

    Parameters
    ----------
    name: str
        The name of the DFT functional or HF method. The same style as
        psi4.energy() function.
    guess_wfn: psi4.core.RHF or psi4.core.UHF, default to None.
        The wfn used to set the initial guess of the SCF calculation. This is
        equivalent to copy guess_wfn.C (coefficient matrix) and guess_wfn.D
        (density matrix) to the scf wavefunction object to initialize the SCF
        calculation. Setting this variable will overwrite the psi4 global guess
        setting.
        !!!
        Note: if you use this variable, make sure you are passing a reasonable
        guess wfn to the function. We do not check the validity the input.
    occ: dict, default to an empty dict.
        A dictionary that specifies the customized occupation number. This
        variable is used to modify the aufbau occupation number of current
        molecule to give the final set of occupation numbers for the SCF
        calculation.

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
        aufbau occupation numbers are:
        alpha occ: 1, 0, 0, 0, ...
        beta occ: 1, 0, 0, 0, ...

        >> customized occ:
        occ = {
            "alpha": {"homo": 0.5}
            "beta": {"3": 0.7}
        }

        >> resulted occ:
        alpha occ: 0.5, 0, 0, 0, ...
        beta occ:  1,   0, 0, 0.7, ...

        Note:
        All the str in this dictionary is case-insensitive.

    verbose: int, default to 1
        print level. 0 means print nothing. 1 means normal print level. A larger
        number means more details.

    Returns
    -------
    wfn: psi4.core.RHF or psi4.core.UHF
        The SCF wavefunction.
        1. Updated data for the returned wfn object:
            wfn.S(): AO overlap matrix.
            wfn.H(): Core matrix.
            wfn.Ca(), wfn.Cb(): CO coefficient matrix.
            wfn.Fa(), wfn.Fb(): Fock matrix.
            wfn.Da(), wfn.Db(): density matrix.
            wfn.Va(), wfn.Vb(): DFA (just the DFA, if it is LOSC-DFA) Vxc matrix.
            wfn.epsilon_a(), wfn.epsilon_b(): orbital energies.
            wfn.energy(): total energy.

        2. Other members may not be initialized. Accessing and using other members
        in the returned wfn is a undefined behavior.

        3. New attributes added to the returned wfn object:
            wfn.elec_a, wfn.elec_b: the total number of electrons.
                Note: if it is a fractional system, these two attributes will
                differ to the returned value of wfn.nalpha() and wfn.nbeta().


    Notes
    -----
    1. Connections to psi4 global setting:
    This function only relies on the GLOBAL settings in psi4, NOT any local
    settings for any psi4 modules!

    2. To build J and K matrices, we directly transforms AO ERI into J and K
    matrix with density matrix. We can not used psi4.JK object to help, because
    psi4.JK object does not provide interface to accept density matrix as input
    (it only accept coefficient matrix so far). Therefore, this function is very
    memory consuming (in N^4 for AO ERI). Only use this function for small
    systems. This is enough, because you will only care the fractional
    calculations for small systems.

    3. AO ERI is constructed by psi4.mintshelper.ao_eri(). So the SCF type is
    always `direct`.

    4. Ignored psi4 global setting:
        scf_type: always ignored.
        guess: ignored when `guess_wfn` is provided and not None.

    5. For the updated matrices/vectors in returned wfn (such as C, D, F matrices),
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
        D[:] = np.einsum('i,ui,vi->uv', np.asarray(occ_val), Cocc, Cocc)


    def update_occ(wfn, occ):
        nbf = wfn.basisset().nbf()
        nelec = [wfn.nalpha(), wfn.nbeta()]
        is_aufbau = True
        is_integer = True
        # Build aufbau occupation.
        rst_occ = [{i: 1 for i in range(n)} for n in nelec]
        for k, v in occ.items():
            k = k.lower()
            spin_chanel = ['alpha', 'beta']
            if k not in spin_chanel:
                raise Exception(f"invalid customized occupation spin chanel: {k}.")
            s = spin_chanel.index(k)
            for orb_i, occ_i in v.items():
                if isinstance(orb_i, str):
                    orb_i = orb_i.lower()
                    if orb_i not in ['homo', 'lumo']:
                        raise Exception(f"unknown customized occupation index: {orb_i}.")
                    if orb_i == 'homo':
                        orb_i = nelec[s] - 1
                    else:
                        orb_i = nelec[s]
                if orb_i >= nbf:
                    raise Exception(f"customized occupation index is out-of-range: {orb_i} (orbital index) > {nbf} (basis size).")
                if not 0 <= occ_i <= 1:
                    raise Exception(f"customized occupation number is invalid: occ={occ_i}.")
                rst_occ[s][orb_i] = occ_i

                if int(occ_i) != float(occ_i):
                    is_integer = False

        occ_idx = []
        occ_val = []
        for d in rst_occ:
            idx_occ = list(d.items())
            idx_occ.sort()
            idx, occ = zip(*idx_occ)
            occ_idx.append(idx)
            occ_val.append(occ)

        nocc = [len(x) for x in occ_idx]

        for s in range(nspin):
            if occ_idx and occ_idx[s][-1] >= nocc[s]:
                is_aufbau = False

        return is_integer, is_aufbau, nocc, occ_idx, occ_val

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
    # TODO:
    # Currenty only support c1 symmetry.
    if mol.schoenflies_symbol() != 'c1':
        raise Exception('Current SCF code only supports C1 symmetry')

    local_print(0, "-----------------------")
    local_print(0, "          SCF")
    local_print(0, "    by Yuncai Mei")
    local_print(0, "-----------------------")

    # Get global settings
    maxiter = psi4.core.get_global_option('MAXITER')
    E_conv = psi4.core.get_global_option('E_CONVERGENCE')
    D_conv = psi4.core.get_global_option('D_CONVERGENCE')
    reference = psi4.core.get_global_option('REFERENCE')
    basis = psi4.core.get_global_option('BASIS')

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
    is_integer, is_aufbau, nocc, occ_idx, occ_val = update_occ(wfn, occ)
    local_print(1, "=> Occupation Number <=")
    local_print(1, f"Is integer system: {is_integer}")
    local_print(1, f"Is aufbau occupation: {is_aufbau}")
    for s in range(nspin):
        if not is_rks:
            local_print(1, f'{"Alpha" if s == 0 else "Beta"} Occupation Number:')
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
    C_psi = [psi4.core.Matrix(nbf, nbf) for s in range(nspin)]
    C = [np.asarray(m) for m in C_psi]
    Cocc = [np.zeros((nbf, nocc[s])) for s in range(nspin)]
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
        local_print(0, "Use psi4.JK to calculate JK matrix.")
        jk_helper = jk.JK_psi4_jk(wfn, C_psi)
    else:
        # To make fractional or non-aufbau calculation available, we need to
        # transform J and K matrices based on density matrix (explicitly
        # constructed based on occupation number). psi4 JK class does not
        # provide interface to accept density matrix as input to build J and K
        # matrices. So here we only the choice to transform from AO ERI to MO
        # ERI for J and K with help of psi4.MinsHelper library.
        local_print(0, "Use psi4.MintsHelper to calculate JK matrix.")
        jk_helper = jk.JK_psi4_mints(wfn, occ_val[:nspin], Cocc[:nspin])

    # ==> Set up initial guess <==
    # Build the initial C matrix as the initial guess.
    C_guess = []
    if guess_wfn:
        local_print(0, "Set initial guess from a reference wavefunction.")
        C_guess = [np.asarray(guess_wfn.Ca()), np.asarray(guess_wfn.Cb())]
    else:
        local_print(0, "Set initial guess from psi4 setting.")
        wfn.form_Shalf()
        wfn.guess()
        C_guess = [np.asarray(wfn.Ca()), np.asarray(wfn.Cb())]
    # Use the initial C guess to update Cocc, C, D matrix.
    for s in range(nspin):
        update_C(occ_idx[s], occ_val[s], C_guess[s], C[s], Cocc[s], D[s])

    # ==> SCF <==
    local_print(0, '\n=> Start SCF iterations <=')
    Eold = 0.0
    Enuc = wfn.molecule().nuclear_repulsion_energy()
    reference_factor = 2.0 if is_rks else 1.0
    Exc = 0.0
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

            # DIIS error build and update
            diis_error[i] = F[i].dot(D[i]).dot(S) - S.dot(D[i]).dot(F[i])
            diis_error[i] = (A.T).dot(diis_error[i]).dot(A)
            diis_helper[i].add(F[i], diis_error[i])

        # Calculate SCF energy.
        if supfunc.needs_xc():
            Exc = Vpot.quadrature_values()['FUNCTIONAL']
        SCF_E = 0
        for i in range(nspin):
            SCF_E += reference_factor * \
                np.einsum('pq,pq->', D[i], H + 0.5 * G[i])
        SCF_E += Enuc + Exc + E_losc

        # Print iteration steps
        dRMS = 0.5 * sum([np.mean(m ** 2) ** 0.5 for m in diis_error])
        local_print(0, '@%3s iter  %3d: Energy = %4.16f   dE = % 1.5E   dRMS = %1.5E  DIIS'
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
    local_print(0, "\n=> SCF Energy <=")
    local_print(0, "SCF energy: {:.10f}".format(SCF_E))

    # print orbital energies
    local_print(0, "\n=> Orbital Energy <=")
    for s in range(nspin):
        if not is_rks:
            local_print(0, f"{'Alpha' if s == 0 else 'Beta'} orbital energy:")
        local_print(0, "{:<5s}  {:<8s}  {:<14s}".format("Index", "Occ", "DFA(eV)"))
        for i in range(nocc[s]):
            local_print(0, "{:<5d}  {:<8.5f}  {:<14.6f}"
                        .format(occ_idx[s][i], occ_val[s][i],
                        eig[s][i] * constants.hartree2ev))
        local_print(0, "")

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

    # Add new attributes to the returned wfn.
    # elec_a: number of alpha electrons.
    # elec_b: number of beta electrons.
    wfn.elec_a, wfn.elec_b = [sum(x) for x in occ_val]

    # restore psi4 options.
    optstash.restore()
    return wfn
