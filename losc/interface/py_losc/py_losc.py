"""
Python interface for localized orbital scaling correction (LOSC) library.

Notations
---------
LO refers to the localized orbital.
AO refers to the atomic orbitals.
CO refers to the canonical orbital.
`nlo` refers to the number of LOs.
`nbasis` refers to the number of AOs.
`nfitbasis` refers to the number of fitting basis for density fitting.
`npts` refers to the number of grids.

Notes
-----
1. All the matrices are represented and stored in 2-dimensional
`np.ndarray` object. The storage order of matrices is required to be
Fortran-style (column major). If the input matrix is not Fortran-style,
it will be converted into Fortran-style with calling `np.asfortranarray()`
function interanally.
"""

import numpy as np
import py_losc_core as core
import functools
from typing import List


def _np2eigen(func):
    """Decorator function that implicitly convert arguments of `func` that
    are `np.ndarray` type into fortran-style, in order to make the storage
    type in python side compatible with `Eigen::MatrixXd` used in LOSC library.
    """
    def to_fortran(t):
        "Convert the variable into a fortran-style np.array if it is np.ndarray."
        if isinstance(t, np.ndarray):
            if not t.flags['F_CONTIGUOUS']:
                t = np.asfortranarray(t)
        elif isinstance(t, list):
            for i in range(len(t)):
                if not t[i].flags['F_CONTIGUOUS']:
                    t[i] = np.asfortranarray(t[i])
        return t

    @functools.wraps(func)
    def wrapper(*args, **kargs):
        n_args = []
        for i in range(len(args)):
            print(f'arg: {i}, {args[i]}')
            n_args.append(to_fortran(args[i]))
        n_kargs = {}
        for k in kargs:
            n_kargs[k] = to_fortran(kargs[k])
        return func(*n_args, **n_kargs)
    return wrapper


# Interface for curvature matrix.
class DFAInfo(core.DFAInfo):
    def __init__(self, gga_x, hf_x, name=''):
        """
        Constructor of DFAInfo class that represents a DFA.

        Parameters
        ----------
        name: str
           The description of the DFA. Default to an empty string.
        hf_x: float
           The weight of HF exchange in the DFA.
        gga_x: float
           The total weights of ALL GGA and LDA type exchange in the DFA.

        Examples
        --------
        B3LYP functional is:
        E_B3LYP = E_LDA_x + a0 * (E_HF_x - E_LDA_x) + ax * (E_GGA_x - E_LDA_x)
        + E_LDA_c + ac * (E_GGA_c - E_LDA_c),
        in which exchanges end with suffix "_x" and correlations end with
        suffix "_c", and (a0, ax, ac) are coefficients.

        The total weights of GGA and LDA exchanges are:
        gga_x = 1 + a0 * (-1) + ax * (1 - 1) = 1 - a0

        The total weights of HF exchanges is clear:
        hf_x = a0
        >>> b3lyp = DFAInfo(gga_x, hf_x, "B3LYP")
        """
        core.DFAInfo.__init__(self, gga_x, hf_x, name)

    def __repr__(self):
        "Representation of DFAInfo object."
        return ("<py_losc.DFAInfo> object: {"
                f"name: {None if self.name() == '' else self.name()}, "
                f"gga_x: {self.gga_x()}, "
                f"hf_x: {self.hf_x()}"
                "}")


class CurvatureV1(core.CurvatureV1):
    """LOSC curvature matrix version 1.

    Curvature version 1 is defined as $\kappa$ in Eq. 10 of the original
    LOSC paper (https://doi.org/10.1093/nsr/nwx11).
    """

    @_np2eigen
    def __init__(self, dfa_info: DFAInfo,
                 C_lo: np.ndarray,
                 df_pii: np.ndarray,
                 df_vpq_inv: np.ndarray,
                 grid_basis_val: np.ndarray,
                 grid_weight: list):
        """
        Curvature version 1 class constructor.

        Parameters
        ----------
        dfa_info: DFAInfo
            The information of the associated DFA.
        C_lo: np.ndarray [nbasis, nlo]
            LOs' coefficient matrix under AOs. The relation between LOs
            ($\phi_i$) and AOs ($\psi_\mu$) via the coefficient matrix
            ($C_{\mu i}$) is the following:
            $\phi_i = \psi_{\mu} C_{\mu i}$.
        df_pii: np.ndarray [nfitbasis, nlo]
            Three-body integral <p|ii> for density fitting, in which index
            p is for fitting basis and index i is for LOs.
        df_Vpq_inv: np.ndarray [nfitbasis, nfitbasis]
            Inverse matrix of integral <p|1/r|q> for density fitting, in
            which index p and q are for fitting basis.
        grid_basis_val: np.ndarray [npts, nbasis]
            Grid values of AOs. grid_basis_val[ip, i] is the i-th AO's value
            on ip-th grid point.
        grid_weight: list or np.ndarray [npts, 1]
            Weights of grids.

        See Also
        --------
        DFAInfo: the DFA information class.
        LoscLocalizerV2: it can build the LOs.
        """
        # Pybind11 has to call the base.__init__() explicitly, instead of
        # using super().__init__() to initialize the base class.
        core.CurvatureV1.__init__(self, dfa_info, C_lo, df_pii, df_vpq_inv,
                                  grid_basis_val, grid_weight)


class CurvatureV2(core.CurvatureV2):
    """LOSC curvature matrix version 2.

    Curvature version 2 is defined as $\tilde{\kappa}$ in Eq. 8 of the
    LOSC2 paper (J. Phys. Chem. Lett. 2020, 11, 4, 1528-1535).
    """

    @_np2eigen
    def __init__(self, dfa_info: core.DFAInfo,
                 C_lo: np.ndarray,
                 df_pii: np.ndarray,
                 df_vpq_inv: np.ndarray,
                 grid_basis_val: np.ndarray,
                 grid_weight: list):
        """
        Curvature version 2 class constructor.

        Notes
        -----
        Same interface for input parameters as described in CurvatureV1 class.

        See Also
        --------
        CurvatureV1: constructor of LOSC curvature version 1.
        LoscLocalizerV2: it can build the LOs.
        """
        core.CurvatureV2.__init__(self, dfa_info, C_lo, df_pii, df_vpq_inv,
                                  grid_basis_val, grid_weight)

# Interface for localization matrix.
class LoscLocalizerV2(core.LoscLocalizerV2):
    """LOSC localization version 2.

    The LOSC localization version 2 is defined in Eq. 7 of the
    LOSC2 paper (J. Phys. Chem. Lett. 2020, 11, 4, 1528-1535).
    """

    @_np2eigen
    def __init__(self, C_lo_basis: np.ndarray,
                 H_ao: np.ndarray,
                 D_ao: List[np.ndarray]):
        """
        LOSC localization version 2 class constructor.

        Parameters
        ----------
        C_lo_basis: np.ndarray [nbasis, nlo]
            The coefficient matrix under AO for the molecular orbitals that
            serves as the basis to expand the LOs. These molecular orbitals
            are called LO basis. The LO basis $\psi_i$ is expanded under AOs
            $\phi_\mu$ via the coefficient matrix $C_{\mu i}$ as
            $\psi_i = \phi_\mu C_{\mu i}$.
        H_ao: np.ndarray [nbasis, nbasis]
            The Halmiltonian matrix under AOs. The LOSC localization v2,
            the Halmiltonian matrix is defined as just the DFA's
            Halmiltonian (note, not the LOSC-DFA Hamiltonian).
        D_ao: list(np.ndarray [nbasis, nbasis])
            The dipole matrix under AOs in the order of x, y and z directions.
        """
        core.LoscLocalizerV2.__init__(self, C_lo_basis, H_ao, D_ao)
