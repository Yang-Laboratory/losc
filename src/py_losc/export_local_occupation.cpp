#include "export_local_occupation.hpp"
#include <losc/local_occupation.hpp>
#include <pybind11/eigen.h>

void export_local_occupation(py::module &m)
{
    // local_occupation
    m.def("local_occupation", &losc::local_occupation,
          R"pddoc(
    Calculate the LOSC local occupation matrix.

    Parameters
    ----------
    C_lo : numpy.array
        The LO coefficient matrix. ``C_lo[:, i]`` is the coefficients of
        the i-th LO represented on AO.
    S: numpy.ndarray
        The AO overlap matrix with dimension [nbasis, nbasis].
    D: numpy.ndarray
        The spin density matrix under AO with dimension [nbasis, nbasis].

    Returns
    -------
    numpy.ndarray
        The local occupation matrix with dimension [nlo, nlo].

    See Also
    --------
    LocalizerV2.lo

    Notes
    -----
    The local occupation matrix is defined as

    .. math:: \lambda_{ij} = \langle \phi_i|\rho|\phi_j \rangle,

    where :math:`\lambda_{ij}` is the local occupation matrix element,
    :math:`\phi_i` and :math:`\phi_j` are the i-th and j-th localized orbital
    and :math:`\rho` is the spin density operator.
    )pddoc");
}
