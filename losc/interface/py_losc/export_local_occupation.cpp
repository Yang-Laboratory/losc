#include "export_local_occupation.h"
#include "../../src/local_occupation.h"
#include <pybind11/eigen.h>

void export_local_occupation(py::module &m)
{
    // local_occupation
    m.def("local_occupation", &losc::local_occupation,
          R"pddoc(
    Calculate the local occupation matrix.

    The local occupation matrix is defined as
    $\lambda_{ij} = \langle \phi_i|\rho|\phi_j \rangle$,
    where $\lambda_{ij}$ is the local occupation matrix element,
    $\phi_i$ and $\phi_j$ are the i-th and j-th localized orbital and
    $\rho$ is the spin density operator. Note that $\lambda$, $\phi$
    and $\rho$ are for the same spin. In matrix form, local occupation
    matrix $L$ is expressed as $L = C_{LO}^T * S * D * S * C_{LO}$.

    Parameters
    ----------
    C_lo: np.ndarray [nbasis, nlo]
        The LO coefficient matrix.
    S: np.ndarray [nbasis, nbasis]
        The AO overlap matrix.
    D: np.ndarray [nbasis, nbasis]
        The spin density matrix under AO.

    Returns
    -------
    out: np.ndarray [nlo, nlo]
        The local occupation matrix.

    See Also
    --------
    LocalizerV2::lo(), ...: Obtain the LOs' coefficient matrix.
    )pddoc");
}
