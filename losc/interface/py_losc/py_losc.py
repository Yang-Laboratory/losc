"""
Python interface for localized orbital scaling correction (LOSC) library.

Notation
-----------
LO refers to the localized orbital.
AO refers to the atomic orbitals.
CO refers to the canonical orbital.
`nlo` refers to the number of LOs.
`nbasis` refers to the number of AOs.
`nfitbasis` refers to the number of fitting basis for density fitting.
`npts` refers to the number of grids.

Note
--------
1. All the matrices are represented and stored in 2-dimensional
`np.ndarray` object. The storage order of matrices is required to be
Fortran-style (column major). If the input matrix is not Fortran-style,
it will be converted into Fortran-style with calling `np.asfortranarray()`
function interanally.
"""

from py_losc_core import DFAInfo
import numpy as np
import py_losc_core as core
import functools


def _np2eigen(func):
    """Decorator function that implicitly convert arguments of `func` that
    are `np.ndarray` type into fortran-style, in order to make the storage
    type in python side compatible with `Eigen::MatrixXd` used in LOSC library.
    """
    def to_fortran(t):
        "Convert the variable into a fortran-style np.array if it is np.ndarray."
        if isinstance(t, np.ndarray):
            if not t.flags['F_CONTIGUOUS']:
                return np.asfortranarray(t)
        return t

    @functools.wraps(func)
    def wrapper(*args, **kargs):
        n_args = []
        for i in range(len(args)):
            n_args.append(to_fortran(args[i]))
        n_kargs = {}
        for k in kargs:
            n_kargs[k] = to_fortran(kargs[k])
        return func(*n_args, **n_kargs)
    return wrapper


def _add_common_attr(cls):
    @property
    def nlo(self):
        "Number of LOs"
        return self._calc_man.nlo()

    @property
    def nbasis(self):
        "Number of AOs"
        return self._calc_man.nbasis()

    @property
    def nfitbasis(self):
        "Number of fit basis"
        return self._calc_man.nfitbasis()

    @property
    def npts(self):
        "Number of grid points"
        return self._calc_man.npts()

    setattr(cls, 'nlo', nlo)
    setattr(cls, 'nbasis', nbasis)
    setattr(cls, 'nfitbasis', nfitbasis)
    setattr(cls, 'npts', npts)
    return cls


# Interface for curvature matrix.
@_add_common_attr
class CurvatureV1():
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
        ------------
        dfa_info: DFAInfo
            The information of the associated DFA.
        C_lo: np.ndarray [nbasis, nlo]
            LO coefficient matrix.
        df_pii: np.ndarray [nfitbasis, nlo]
            Three-body integral <p|ii> for density fitting.
        df_Vpq_inv: np.ndarray [nfitbasis, nfitbasis]
            Inverse matrix of integral <p|1/r|q> for density fitting.
        grid_basis_val: np.ndarray [npts, nbasis]
            Grid values of AOs.
        grid_weight: list or np.ndarray [npts, 1]
            Weights of grids.
        """
        self._calc_man = core.CurvatureV1(dfa_info, C_lo, df_pii,
                                          df_vpq_inv, grid_basis_val,
                                          grid_weight)

    def kappa(self) -> np.ndarray:
        """
        Compute the curvature matrix version 1.

        Return
        -------
        np.ndarray [nlo, nlo]
            The curvature matrix.
        """
        return self._calc_man.kappa()


@_add_common_attr
class CurvatureV2():
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

        Parameters
        ------------
        Same interface for input parameters as described in CurvatureV1 class.
        """
        self._calc_man = core.CurvatureV2(dfa_info, C_lo, df_pii,
                                          df_vpq_inv, grid_basis_val,
                                          grid_weight)

    def kappa(self):
        """
        Compute the curvature matrix version 2.

        Return
        -------
        np.ndarray [nlo, nlo]
            The curvature matrix.
        """
        return self._calc_man.kappa()

# TODO: Interface for localization matrix.
