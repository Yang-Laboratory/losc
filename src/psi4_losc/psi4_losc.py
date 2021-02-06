"""
Integrate py_losc module with psi4 package to perform LOSC calculation.
"""

from pkg_resources import parse_version
from psi4_losc import diis
import psi4
import numpy as np
from psi4.driver.p4util.exceptions import ValidationError
from py_losc import py_losc

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
        wfn.molecule(), "DF_BASIS_SCF", "", "JKFIT", "aug-cc-pVDZ")
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
    npts = 0
    # number of LOs.
    nlo = C_lo.shape[1]
    # Get the psi4.potential object for help
    # TODO: check if points_func is related to threads.
    points_func = Vpot.properties()[0]
    for b in range(Vpot.nblocks()):
        npts += Vpot.get_block(b).npoints()

    # allocate grid_lo matrices.
    grid_lo = np.zeros((npts, nlo))

    # loop over the blocks to build grid_lo
    npts_count = 0
    for b in range(Vpot.nblocks()):
        # Obtain block information
        block = Vpot.get_block(b)
        points_func.compute_points(block)
        npoints = block.npoints()

        # compute grid_ao on the fly to build grid_lo
        grid_ao = np.array(points_func.basis_values()["PHI"])[:npoints, :]
        grid_lo_blk = grid_lo[npts_count:npts_count+npoints, :]
        grid_lo_blk[:] = grid_ao.dot(C_lo)

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
    npts = 0
    # Get the psi4.potential object for help
    # TODO: check if points_func is related to threads.
    points_func = Vpot.properties()[0]
    for b in range(Vpot.nblocks()):
        npts += Vpot.get_block(b).npoints()

    # build grid_w vector
    grid_w = np.zeros((npts,))
    npts_count = 0
    # loop over the blocks to build grid_w
    for b in range(Vpot.nblocks()):
        # Obtain block information
        block = Vpot.get_block(b)
        # TODO: check do we need compute points to get weights?
        points_func.compute_points(block)
        npoints = block.npoints()

        # view of the grid_w block
        grid_w_blk = grid_w[npts_count:npts_count+npoints]

        # Obtain the grid weight
        grid_w_blk[:] = np.array(block.w())

    return grid_w


def post_scf_losc(dfa_wfn):
    """
    Perform the post-SCF-LOSC calculation.

    Parameters
    ----------
    dfa_wfn: psi4.core.HF like psi4 object
        The converged wavefunction from a parent DFA.

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

    _print_header(scf_losc=False)
    is_rks = dfa_wfn.same_a_b_orbs() and dfa_wfn.same_ab_dens()
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
        # TODO: check if I use np.einsum right.
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
        losc_eig[s] = py_losc.orbital_energy_post_scf(
            H_ao[s], H_losc[s], C_co[s])

    return dfa_wfn.energy() + sum(E_losc), losc_eig


def scf_losc(dfa_wfn, inplace=True):
    """
    Perform the SCF-LOSC (frozen-LO) calculation.

    Parameters
    ----------
    dfa_wfn: psi4.core.HF like psi4 object
        The converged wavefunction from a parent DFA. It will be modified if the
        `inplace` is set to True.
    inplace: bool, default to True.
        Modify the `dfa_wfn` in place.

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
        C_lo[s] = localizer.lo("identity")

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
        # TODO: check if I use np.einsum right.
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
    dfa_eigs = [np.asarray(dfa_wfn.epsilon_a()),
                np.asarray(dfa_wfn.epsilon_b())]
    for spin in range(nspin):
        psi4.core.print_out(
            f'\n{"Alpha" if spin == 0 else "Beta"} orbital energies in a.u.\n')
        header = "{:<6s} {:>20s} {:>20s}\n".format(
            "Index", "DFA", "LOSC-DFA")
        psi4.core.print_out(header)
        for i in range(len(dfa_eigs[spin])):
            psi4.core.print_out('{:<6d} {:>20.8f} {:>20.8f}\n'.format(
                i, dfa_eigs[spin][i], eig[spin][i]))
    return wfn
