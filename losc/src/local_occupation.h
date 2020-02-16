#ifndef _LOSC_SRC_LOCAL_OCCUPATION_H_
#define _LOSC_SRC_LOCAL_OCCUPATION_H_

#include <matrix/matrix.h>
#include <memory>
#include "exception.h"

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
 * @ param [in] C_lo: LO coefficient matrix with dimension [nlo, nbasis],
 *   which can be generated from `losc::LocalizerBase::compute()`.
 * @ param [in] S: AO overlap matrix with dimension [nbasis, nbasis].
 * @ param [in] D: spin density matrix under AO with dimension [nbasis, nbasis].
 * @ return : the local occupation matrix with dimension [nlo, nlo].
 */
inline SharedMatrix local_occupation_matrix(const Matrix& C_lo, const Matrix& S, const Matrix& D)
{
    size_t nlo = C_lo.row();
    size_t nbasis = C_lo.col();

    if (nlo > nbasis) {
        throw exception::DimensionError("wrong dimension for LO coefficient matrix: number of LO is larger than the number of AO.");
    }
    if (!S.is_square() || nbasis != S.row()) {
        throw exception::DimensionError(S, nbasis, nbasis, "wrong dimension for AO overlap matrix.");
    }
    if (!D.is_square() || nbasis != D.row()) {
        throw exception::DimensionError(D, nbasis, nbasis, "wrong dimension for density matrix.");
    }

    Matrix SDS(nbasis, nbasis);
    SharedMatrix L = std::make_shared<Matrix> (nlo, nlo);
    matrix::mult_dgemm_ABAT(S, D, SDS);
    matrix::mult_dgemm_ABAT(C_lo, SDS, *L);
    return L;
}

}   // namespace losc
#endif // _LOSC_SRC_LOCAL_OCCUPATION_H_
