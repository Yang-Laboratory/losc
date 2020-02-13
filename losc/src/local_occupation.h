#ifndef _LOSC_LOCAL_OCCUPATION_H_
#define _LOSC_LOCAL_OCCUPATION_H_

#include <matrix/matrix.h>
#include <memory>

namespace losc {

using matrix::Matrix;
using SharedMatrix = std::shared_ptr<Matrix>;

/**
 * Calculate the local occupation matrix.
 *
 * The local occupation matrix is defined as
 * $\lambda_{ij} = <\phi_i|\rho|\phi_j>$, where $\lambda_{ij}$ is the local occupation matrix
 * element, $\phi_i$ and $\phi_j$ are the i-th and j-th localized orbital and $\rho$
 * is the spin density operator. Note that $\lambda$, $\phi$ and $\rho$ are for
 * the same spin. In matrix form, local occupation matrix $L$ is expressed as
 * $L = C_{LO} * S * D * S * C_{LO}^T$.
 *
 * @ param [in] C_lo: LO coefficient matrix which can be generated from `localization::LocalizerBase`.
 * @ param [in] S: atomic orbital (basis) overlap matrix. It is a full matrix with
 *      dimension (nbasis, nbasis).
 * @ param [in] D: spin density matrix under AO. It is a full matrix with dimension (nbasis, nbasis).
 * @ return L: the local occupation matrix.
 */
inline SharedMatrix local_occupation_matrix(SharedMatrix C_lo, SharedMatrix S, SharedMatrix D)
{
    size_t nlo = C_lo->row();
    size_t nbasis = C_lo->col();
    Matrix SDS(nbasis, nbasis);
    SharedMatrix L = std::make_shared<Matrix> (nlo, nlo);
    matrix::mult_dgemm_ABAT(*S, *D, SDS);
    matrix::mult_dgemm_ABAT(*C_lo, SDS, *L);
    return L;
}

}

#endif // _LOSC_LOCAL_OCCUPATION_H_
