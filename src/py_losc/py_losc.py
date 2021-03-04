import numpy as np
from py_losc import py_losc_core as core
from typing import List

# ==> Interface for LOSC curvature matrix <==

def _convert_mat(m, order='C'):
    """
    Convert input np.array into the specified order.

    Parameters
    ----------
    m: np.array
        The array to convert.
    order: str, choices = ['C', 'F'], default to 'C'.
        'C': C-style row major order.
        'F': Fortran-style column major order.

    Returns
    -------
    out: np.array
        The array with storage order equaling to `order`.
    """
    if order == 'C':
        return np.ascontiguousarray(m)
    elif order == 'F':
        return np.asfortranarray(m)
    else:
        raise Exception(f'Unknown storage order: {order}')

class DFAInfo(core.DFAInfo):
    """The information of a density functional approximation.
    """

    def __init__(self, gga_x, hf_x, name=''):
        """
        Constructor of DFAInfo class that represents a DFA.

        Parameters
        ----------
        name : str
           The description of the DFA. Default to an empty string.
        hf_x : float
           The weight of HF exchange in the DFA.
        gga_x : float
           The total weights of ALL GGA and LDA type exchange in the DFA.

        Examples
        --------
        B3LYP functional is

        .. math:: E^{\\rm{B3LYP}}_{\\rm{xc}} = E^{\\rm{LDA}}_{\\rm{x}} +
           a_0 (E^{\\rm{HF}}_{\\rm{x}} - E^{\\rm{LDA}}_{\\rm{x}}) +
           a_x (E^{\\rm{GGA}}_{\\rm{x}} - E^{\\rm{LDA}}_{\\rm{x}}) +
           E^{\\rm{LDA}}_{\\rm{c}} + a_c (E^{\\rm{GGA}}_{\\rm{c}} -
           E^{\\rm{LDA}}_{\\rm{c}}),

        in which exchanges end with suffix ``_x`` and correlations end with
        suffix ``_c``, :math:`a_0=0.20`, :math:`a_x=0.72` and :math:`a_c=0.81`.
        The total weights of GGA and LDA exchanges are
        :math:`1 - a_0 + a_x \times (1 - 1) = 1 - a_0 = 0.80`.
        The total weights of HF exchanges is
        :math:`a_0 = 0.20`. To construct a ``DFAInfo`` object for B3LYP
        functional, one should do the following

        .. code-block:: python

            >>> b3lyp = DFAInfo(0.80, 0.20, "B3LYP")
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

    References
    ----------
    Eq. 10 in the original LOSC paper https://doi.org/10.1093/nsr/nwx11.
    """

    def __init__(self, dfa_info: DFAInfo,
                 df_pii: np.ndarray,
                 df_vpq_inv: np.ndarray,
                 grid_lo: np.ndarray,
                 grid_weight: np.ndarray):
        """
        Constructor of LOSC curvature version1.

        Parameters
        ----------
        dfa_info : DFAInfo
            The information for the DFA.
        df_pii : numpy.array
            The three-center integral :math:`\langle p | ii \\rangle`
            used in density fitting, in which index `p` refers to the =
            fitting basis and index `i` refers to the LO. The dimension of
            `df_pii` is [nfitbasis, nlo].
        df_Vpq_inv : numpy.array
            The inverse of integral matrix
            :math:`\langle p | 1/\mathbf{r} | q \\rangle` used in density
            fitting, in which indices `p` and `q` refer to the fitting basis.
            The dimension of `df_Vpq_inv` is [nfitbasis, nfitbasis].
        grid_lo : numpy.array
            The LOs values on grid points. The dimension of `grid_lo` is
            [npts, nlo].
        grid_weight : numpy.array
            The weights for numerical integral over grid points. The
            dimension of `grid_lo` is [npts, ].

        See Also
        --------
        DFAInfo
        LocalizerV2
        """
        # Pybind11 has to call the base.__init__() explicitly, instead of
        # using super().__init__() to initialize the base class.
        self._df_pii = _convert_mat(df_pii)
        self._df_vpq_inv = _convert_mat(df_vpq_inv)
        self._grid_lo = _convert_mat(grid_lo)
        self._grid_wt = _convert_mat(grid_weight)
        core.CurvatureV1.__init__(self, dfa_info, self._df_pii,
                                  self._df_vpq_inv, self._grid_lo,
                                  self._grid_wt)


class CurvatureV2(core.CurvatureV2):
    """LOSC curvature matrix version 2.

    References
    ----------
    Eq. 8 in the LOSC2 paper (J. Phys. Chem. Lett. 2020, 11, 4, 1528-1535).
    """

    # @_np2eigen
    def __init__(self, dfa_info: core.DFAInfo,
                 df_pii: np.ndarray,
                 df_vpq_inv: np.ndarray,
                 grid_lo: np.ndarray,
                 grid_weight: np.ndarray):
        """

        Curvature version 2 class constructor.

        See Also
        --------
        CurvatureV1, LocalizerV2

        Notes
        -----
        Same interface for input parameters as `CurvatureV1.__init__`.

        """
        self._df_pii = _convert_mat(df_pii)
        self._df_vpq_inv = _convert_mat(df_vpq_inv)
        self._grid_lo = _convert_mat(grid_lo)
        self._grid_wt = _convert_mat(grid_weight)
        core.CurvatureV2.__init__(self, dfa_info, self._df_pii,
                                  self._df_vpq_inv, self._grid_lo,
                                  self._grid_wt)


# ==> Interface for LOSC localization <==
class LocalizerV2(core.LocalizerV2):
    """LOSC localization version 2.

    References
    ----------
    Eq. 7 in the LOSC2 paper (J. Phys. Chem. Lett. 2020, 11, 4, 1528-1535).
    """

    def __init__(self, C_lo_basis: np.ndarray,
                 H_ao: np.ndarray,
                 D_ao: List[np.ndarray]):
        """
        Constructor of LOSC localizerV2.

        Parameters
        ----------
        C_lo_basis : numpy.array
            The coefficient matrix of LO basis that is expanded under AO with
            dimension [nbasis, nlo]. ``C_lo_basis[:, i]`` is the coefficient of
            the i-th LO basis. The LO basis refers to a set of MOs that
            will be unitary transformed to be the LOs. The basis is usually
            being the COs from the associated DFA.

        H_ao : numpy.array
            The Hamiltonian matrix under AO representation for the associated
            DFA. The dimension is [nbasis, nbasis].

        D_ao : list of numpy.array
            The dipole matrix under AO representation for x, y and z directions.
        """
        self._C_lo_basis = _convert_mat(C_lo_basis)
        self._H_ao = _convert_mat(H_ao)
        self._D_ao = [_convert_mat(x) for x in D_ao]
        core.LocalizerV2.__init__(
            self, self._C_lo_basis, self._H_ao, self._D_ao)
