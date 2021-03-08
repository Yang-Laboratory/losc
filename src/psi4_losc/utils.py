"""
Psi4_losc package utils, a collection of helper functions to help to integrate
LOSC implementation with psi4.
"""

import psi4
import numpy as np
from random import shuffle
from qcelemental import constants
from psi4_losc import options


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


def print_orbital_energies(verbose_level, wfn, losc_data, print_level=1,
                           window=None):
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
            orbE = dfa_eigs[s][i]
            if orbital_energy_unit != 'eV':
                orbE *= constants.hartree2ev
            if not window or (window[0] <= orbE <= window[1]):
                local_print(print_level, "{:<5d}  {:<8.5f}  {:>14.6f} {:>14.6f}"
                            .format(i, occ_idx_val[s].get(i, 0),
                                    dfa_eigs[s][i], losc_eigs[s][i]))
        local_print(print_level, "")


def print_full_matrix(mat, line_limit=5):
    """
    Print matrix into psi4 output.
    """
    psi4.core.print_out(f'Dimension: {mat.shape}\n')
    row, col = mat.shape
    for i in range(row):
        psi4.core.print_out(f'    {i}:\n')
        for j in range(col):
            psi4.core.print_out(' {:>15.8e}'.format(mat[i, j]))
            if j != col - 1:
                if (j+1) % line_limit == 0:
                    psi4.core.print_out('\n')
        psi4.core.print_out('\n')


def print_sym_matrix(mat, line_limit=5):
    """
    Print matrix into psi4 output.
    """
    psi4.core.print_out(f'Dimension: {mat.shape}\n')
    row, col = mat.shape
    for i in range(row):
        psi4.core.print_out(f'    {i}:\n')
        for j in range(0, i + 1):
            psi4.core.print_out(' {:>15.8e}'.format(mat[i, j]))
            if j != i:
                if (j+1) % line_limit == 0:
                    psi4.core.print_out('\n')
        psi4.core.print_out('\n')


def form_df_matrix(wfn, C_lo):
    """
    Build density fitting related matrices. DF_pii, DF_Vpq_inv.

    Parameters
    ----------
    wfn: psi4.core.HF
        The wavefunction object.
    C_lo: [np.array, ...]
        LOSC LO coefficient matrix.

    Returns
    -------
    df_pii: [np.array, ...]
        The three center <fitbasis|lo, lo> matrix.
    df_Vpq_inv: np.array
        The inverse of <fitbasis|1/r|fitbasis> matrix.
    """
    # basis set
    basis = wfn.basisset()
    zero_bas = psi4.core.BasisSet.zero_ao_basis_set()
    # psi4.mintshelper object to help building AO integrals.
    mints = psi4.core.MintsHelper(basis)

    # build auxiliary fitbasis for each fragments.
    nspin = len(C_lo)
    mol = wfn.molecule()
    frag_mol, whole_mol = split_molecule(mol, return_whole_mol=True,
                                         frag_size=options.get_param('curvature', 'df_molecular_fragment_size'))
    aux_bas_name = psi4.core.get_global_option('DF_BASIS_SCF').lower()
    aux_bas = []
    nfitbasis = 0
    for frag in frag_mol:
        aux_bas_t = psi4.core.BasisSet.build(frag, "ORBITAL", aux_bas_name)
        nfitbasis += aux_bas_t.nbf()
        aux_bas.append(aux_bas_t)

    # ==> build df_pii <==
    df_pii = [np.zeros((nfitbasis, C_lo[s].shape[1])) for s in range(nspin)]
    nfitbasis_count = 0
    for i in range(len(aux_bas)):
        # build three-center integral <fitbasis|ao, ao>
        df_pmn = np.asarray(mints.ao_eri(aux_bas[i], zero_bas, basis, basis))
        df_pmn = np.squeeze(df_pmn)

        # transform from AO to LO.
        nbf = aux_bas[i].nbf()
        for s in range(nspin):
            # build three-center integral <fitbasis|lo, lo>
            df_pii[s][nfitbasis_count:nfitbasis_count + nbf, :] = (
                np.einsum('pmn,mi,ni->pi', df_pmn,
                          C_lo[s], C_lo[s], optimize=True))

        # update row counter
        nfitbasis_count += nbf

    # ==> build df_Vpq_inv <==
    aux_bas_all = psi4.core.BasisSet.build(whole_mol, "ORBITAL", aux_bas_name)
    df_Vpq = np.asarray(mints.ao_eri(
        aux_bas_all, zero_bas, aux_bas_all, zero_bas))
    df_Vpq = np.squeeze(df_Vpq)
    df_Vpq_inv = np.linalg.inv(df_Vpq)

    return df_pii, df_Vpq_inv


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


def split_molecule(mol, frag_size=5, return_whole_mol=False):
    """
    Split a molecule into small fragments based on number of atoms.

    This function is based on the `psi4.core.Molecule.to_dict()` functions to
    extract the information of molecular xyz and atom labels. At the return,
    a list of molecular fragments are created by calling `psi4.geometry()`
    function.

    Parameters
    ----------
    mol: psi4.core.Molecule
        The input molecule to be splitted.
    frag_size: int, default to 5.
        The maximum size of fragments in atom numbers.
    return_whole_mol: bool, default to False.
        Return the whole molecule which has the same order of returned fragments.
        The order needs consideration because the generated fragments are
        shuffled.

    Returns
    -------
    mol_frags: [psi4.core.Molecule, ...]
        A list of small molecular fragments.
    whole_mol: psi4.core.Molecule
        The whole molecule with the corresponding order for splitted fragmetns.
    """
    mol_dict = mol.to_dict()
    unit = mol_dict['units']
    xyz = mol_dict['geom'].reshape(-1, 3)
    atoms = mol_dict['elem']

    def make_mol(geom_str_list):
        # specify no re-orientation, no re-centering and use the same units.
        geom = geom_str_list.copy()
        geom.append('no_reorient')
        geom.append('no_com')
        geom.append(f'units {unit}')
        return psi4.geometry('\n'.join(geom))

    if len(xyz) != len(atoms):
        raise Exception('Miss some atoms.')

    # form xyz string to represent the molecule.
    geom_str = []
    for i in range(len(atoms)):
        atom_xyz = ['{:.16f}'.format(x) for x in xyz[i, :]]
        geom_str.append(' '.join([str(atoms[i])] + atom_xyz))

    # shuffle atoms before spliting into chunks
    shuffle(geom_str)
    whole_mol = make_mol(geom_str)

    # split into molecular fragments
    mol_frags = []
    for i in range(0, len(geom_str), frag_size):
        frag = geom_str[i:i+frag_size]
        mol_frags.append(make_mol(frag))

    # we need to reset the wfn.molecule to be the active one.
    psi4.core.set_active_molecule(mol)

    if return_whole_mol:
        return mol_frags, whole_mol
    else:
        return mol_frags
