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
from py_losc import py_losc_core as core
import functools
from typing import List

# ==> Interface for LOSC corrections <==
from py_losc.py_losc_core import ao_hamiltonian_correction
from py_losc.py_losc_core import orbital_energy_post_scf
from py_losc.py_losc_core import energy_correction

# ==> Interface for LOSC local occupation <==
from py_losc.py_losc_core import local_occupation

# ==> Interface for LOSC curvature matrix <==


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

    def __init__(self, dfa_info: DFAInfo,
                 df_pii: np.ndarray,
                 df_vpq_inv: np.ndarray,
                 grid_lo: np.ndarray,
                 grid_weight: np.ndarray):
        """
        Curvature version 1 class constructor.

        Parameters
        ----------
        dfa_info: DFAInfo
            The information of the associated DFA.
        df_pii: np.ndarray [nfitbasis, nlo]
            Three-body integral <p|ii> for density fitting, in which index
            p is for fitting basis and index i is for LOs.
        df_Vpq_inv: np.ndarray [nfitbasis, nfitbasis]
            Inverse matrix of integral <p|1/r|q> for density fitting, in
            which index p and q are for fitting basis.
        grid_lo: np.ndarray [npts, nbasis]
            Grid values of LOs. grid_lo[ip, i] is the i-th LO's value
            on ip-th grid point.
        grid_weight: list or np.ndarray [npts, 1]
            Weights of grids.

        See Also
        --------
        DFAInfo: the DFA information class.
        LocalizerV2: it can build the LOs.
        """
        # Pybind11 has to call the base.__init__() explicitly, instead of
        # using super().__init__() to initialize the base class.
        self._df_pii = np.asfortranarray(df_pii)
        self._df_vpq_inv = np.asfortranarray(df_vpq_inv)
        self._grid_lo = np.asfortranarray(grid_lo)
        self._grid_wt = np.asarray(grid_weight)
        self._grid_wt = np.asfortranarray(self._grid_wt)
        core.CurvatureV1.__init__(self, dfa_info, self._df_pii,
                                  self._df_vpq_inv, self._grid_lo,
                                  self._grid_weight)


class CurvatureV2(core.CurvatureV2):
    """LOSC curvature matrix version 2.

    Curvature version 2 is defined as $\tilde{\kappa}$ in Eq. 8 of the
    LOSC2 paper (J. Phys. Chem. Lett. 2020, 11, 4, 1528-1535).
    """

    # @_np2eigen
    def __init__(self, dfa_info: core.DFAInfo,
                 df_pii: np.ndarray,
                 df_vpq_inv: np.ndarray,
                 grid_lo: np.ndarray,
                 grid_weight: np.ndarray):
        """
        Curvature version 2 class constructor.

        Notes
        -----
        Same interface for input parameters as described in CurvatureV1 class.

        See Also
        --------
        CurvatureV1: constructor of LOSC curvature version 1.
        LocalizerV2: it can build the LOs.
        """
        self._df_pii = np.asfortranarray(df_pii)
        self._df_vpq_inv = np.asfortranarray(df_vpq_inv)
        self._grid_lo = np.asfortranarray(grid_lo)
        self._grid_wt = np.asarray(grid_weight)
        self._grid_wt = np.asfortranarray(self._grid_wt)
        core.CurvatureV2.__init__(self, dfa_info, self._df_pii,
                                  self._df_vpq_inv, self._grid_lo,
                                  self._grid_wt)


# ==> Interface for LOSC localization <==
class LocalizerV2(core.LocalizerV2):
    """LOSC localization version 2.

    The LOSC localization version 2 is defined in Eq. 7 of the
    LOSC2 paper (J. Phys. Chem. Lett. 2020, 11, 4, 1528-1535).
    """

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
        self._C_lo_basis = np.asfortranarray(C_lo_basis)
        self._H_ao = np.asfortranarray(H_ao)
        self._D_ao = [np.asfortranarray(x) for x in D_ao]
        core.LocalizerV2.__init__(
            self, self._C_lo_basis, self._H_ao, self._D_ao)
