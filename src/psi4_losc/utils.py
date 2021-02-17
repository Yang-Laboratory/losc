"""
Psi4_losc package utils, a collection of helper functions to help to integrate
LOSC implementation with psi4.
"""

import psi4
import numpy as np


def _local_print(print_level, verbose_level, *args):
    """
    Print arguments to psi4 output file when verbose level >= print level.

    Parameters
    ----------
    print_level: int
        The print level.
    verbose_level: int
        The verbose level.
    """
    if verbose_level >= print_level:
        t = [f'{i}' for i in args]
        psi4.core.print_out(f'{" ".join(t)}\n')


def init_local_print(verbose_level):
    """
    Create a local printer and set its verbose level.

    Parameters
    ----------
    verbose_level: int
        The verbose level that is compared to print level.

    Returns
    -------
    local_print: function(print_level, *args)
        A printer with setting verbose level to be `verbose_level`.
    """
    def local_print(print_level, *args):
        """
        Print arguments into psi4 output file when verbose level >= print level.
        The verbose level is initialized at the call of
        `psi4_losc.utils.init_local_print()`.
        """
        _local_print(print_level, verbose_level, *args)
    return local_print


def print_total_energies(verbose_level, losc_data, print_level=1):
    """
    Print total energies to pis4 output file. Default print level is 1.
    """
    local_print = init_local_print(verbose_level)
    l = print_level
    local_print(l, f'\n=> {losc_data["losc_type"]} Energetics <=')
    local_print(l, 'DFA energy: {:.16f}'.format(losc_data['dfa_energy']))
    local_print(l, 'LOSC correction energy: {:.16f}'.format(
        losc_data['losc_energy']))
    local_print(
        l, 'LOSC-DFA total energy: {:.16f}'.format(losc_data['losc_dfa_energy']))


def print_orbital_energies(verbose_level, wfn, losc_data, print_level=1):
    """
    Print orbital energies to pis4 output file. Default print level is 1.
    """
    local_print = init_local_print(verbose_level)

    nspin = losc_data['nspin']
    orbital_energy_unit = losc_data['orbital_energy_unit']
    dfa_eigs = losc_data['dfa_orbital_energy']
    losc_eigs = losc_data['losc_dfa_orbital_energy']

    nbf = wfn.basisset().nbf()
    is_rks = True if nspin == 1 else False
    _, occ_idx, occ_val = form_occ(wfn, losc_data['occ'])
    occ_idx_val = [dict(zip(occ_idx[s], occ_val[s])) for s in range(nspin)]

    local_print(
        print_level, f'\n=> {losc_data["losc_type"]} Orbital Energies <=')
    for s in range(nspin):
        if not is_rks:
            local_print(
                print_level, f"{'Alpha' if s == 0 else 'Beta'} orbital energy:")
        local_print(print_level, "{:<5s}  {:<8s}  {:>14s} {:>14s}"
                    .format("Index", "Occ", f"DFA ({orbital_energy_unit})",
                            f"LOSC ({orbital_energy_unit})"))
        for i in range(nbf):
            local_print(print_level, "{:<5d}  {:<8.5f}  {:>14.6f} {:>14.6f}"
                        .format(i, occ_idx_val[s].get(i, 0),
                                dfa_eigs[s][i], losc_eigs[s][i]))
        local_print(print_level, "")


def form_df_basis_matrix(wfn):
    """
    Form the three-center integral used in density fitting.

    Parameters
    ----------
    wfn: psi4.core.wavefunction
        A psi4 wavefunction object.

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
    aux_bas_name = psi4.core.get_global_option('DF_BASIS_SCF').lower()
    aux_bas = psi4.core.BasisSet.build(wfn.molecule(), "ORBITAL", aux_bas_name)
    aux_bas.print_out()
    # psi4.mintshelper object to help building AO integrals.
    mints = psi4.core.MintsHelper(basis)
    # build three-center integral <fitbasis|ao, ao>
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
        A psi4 HF object that is associated with a DFT calculations.
    C_lo: np.array [nbasis, nlo]
        LO coefficient matrix on AOs.

    Returns
    -------
    grid_lo: np.array [npts, nlo]
        A matrix for the values of LOs on grid points. npts is the number of
        grid points, and nlo is the number of LOs.
    """
    if not isinstance(wfn, psi4.core.HF):
        raise Exception(
            "Unknown type of argument wfn. It has to be the type of psi4.core.HF.")
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
        grid_ao = np.array(points_func.basis_values()["PHI"])[
            :npoints, :lpos.shape[0]]
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
        A psi4 HF object that is associated with a DFT calculations.

    Returns
    -------
    grid_w: np.array [npts]
        A vector for the weights of grid points.


    References
    ----------
    psi4numpy/Tutorials/04_Density_Functional_Theory/4b_LDA_kernel.ipynb
    """
    if not isinstance(wfn, psi4.core.HF):
        raise Exception(
            "Unknown type of argument wfn. It has to be the type of psi4.core.HF.")
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


def form_occ(wfn, occ={}):
    """
    Form the occupation number for a given psi4 wavefunction object.

    Parameters
    ----------
    wfn: psi4.core.Wavefunction
        A psi4 wavefunction object.
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
        customized occupation number (in the range of [0, 1]). Using `homo` and
        the orbital index of HOMO as the key will cause ambiguity to specify
        the HOMO occupation number, so it is not allowed. The same rule applies
        to LUMO.

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

        Returns
        -------
        nocc: [int, int]
            Number of non-zero occupied orbitals.
        occ_idx: [[int, ...], [int, ...]]
            The index (0-based) of non-zero occupied orbitals. If the alpha or
            beta spin chanel has zero number of non-zero occupied orbitals,
            the index list is an empty list.
        occ_val: [[float, ...], [float, ...]]
            The occupation number of non-zero occupied orbitals. If the alpha or
            beta spin chanel has zero number of non-zero occupied orbitals,
            the occupation number list is an empty list.
    """
    nbf = wfn.basisset().nbf()
    nelec = [wfn.nalpha(), wfn.nbeta()]
    # Build aufbau occupation.
    rst_occ = [{i: 1 for i in range(n)} for n in nelec]
    for k, v in occ.items():
        k = k.lower()
        spin_chanel = ['alpha', 'beta']
        if k not in spin_chanel:
            raise Exception(f"invalid customized occupation spin chanel: {k}.")
        s = spin_chanel.index(k)
        homo_status = []
        lumo_status = []
        for orb_i, occ_i in v.items():
            # check orbital index
            if isinstance(orb_i, str):
                orb_i = orb_i.lower()
                if orb_i not in ['homo', 'lumo']:
                    raise Exception(
                        f"unknown customized occupation index: {orb_i}.")
                if orb_i == 'homo':
                    orb_i = nelec[s] - 1
                else:
                    orb_i = nelec[s]
            elif not isinstance(orb_i, int):
                raise Exception(
                    'Orbital index should be either int, "homo" or "lumo".')
            if not 0 <= orb_i < nbf:
                raise Exception(
                    f"customized occupation index is out-of-range [0, {nbf}): {orb_i}")
            # check occupation number
            if not 0 <= occ_i <= 1:
                raise Exception(
                    f"customized occupation number is out-of-range [0, 1]: occ={occ_i}.")
            # check if homo and homo_idx appear at the same time:
            if orb_i == nelec[s] - 1:
                homo_status.append(orb_i)
            if len(homo_status) == 2:
                raise ValueError(
                    "HOMO and HOMO index appear at the same time, which causes ambiguity.")
            # check if lumo and lumo_idx appear at the same time:
            if orb_i == nelec[s]:
                lumo_status.append(orb_i)
            if len(lumo_status) == 2:
                raise ValueError(
                    "LUMO and LUMO index appear at the same time, which causes ambiguity.")

            # update occupation number
            rst_occ[s][orb_i] = occ_i

    occ_idx = []
    occ_val = []
    for d in rst_occ:
        if d:
            idx_occ = list(d.items())
            # remove 0 occupied orbital
            idx_occ = list(filter(lambda x: abs(x[1]) > 0, idx_occ))
            idx_occ.sort()
            if idx_occ:
                idx, occ = zip(*idx_occ)
            else:
                idx, occ = tuple(), tuple()
            occ_idx.append(idx)
            occ_val.append(occ)
        else:
            occ_idx.append(tuple())
            occ_val.append(tuple())

    nocc = tuple(len(x) for x in occ_idx)

    return nocc, occ_idx, occ_val


def is_integer_system(wfn, occ={}):
    """
    Check if the wavefunction is an integer system.

    Parameters
    ----------
    wfn: psi4.core.Wavefunction
        The wavefunction object to be check.
    occ: dict
        The occupation information dictionary. See `form_occ()`.

    Returns
    -------
    out: bool
        If the wavefunction is an integer system.
    """
    _, _, occ_val = form_occ(wfn, occ=occ)
    for occ_val_i in occ_val:
        for o in occ_val_i:
            if int(o) != float(o):
                return False
    return True


def is_aufbau_system(wfn, occ={}):
    """
    Check if the wavefunction is an aufbau system.

    Parameters
    ----------
    wfn: psi4.core.Wavefunction
        The wavefunction object to be check.
    occ: dict
        The occupation information dictionary. See `form_occ()`.

    Returns
    -------
    out: bool
        If the wavefunction is an aufbau system.
    """
    norbs, occ_idx, occ_val = form_occ(wfn, occ=occ)
    nelec = [sum(x) for x in occ_val]
    occ_idx_val = [list(zip(occ_idx[s], occ_val[s]))
                   for s in range(len(norbs))]
    for s in range(len(norbs)):
        n_orbitals = norbs[s]
        if n_orbitals > 0:
            highest_orbital_idx = occ_idx_val[s][-1][0]
            electron_number = nelec[s]
            is_aufbau = (n_orbitals - 1 <= electron_number <= n_orbitals and
                         highest_orbital_idx == n_orbitals - 1)
            if not is_aufbau:
                return False
    return True
