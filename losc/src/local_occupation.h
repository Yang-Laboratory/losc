/**
 * @file
 * @brief Losc localization occupation matrix related declarations.
 */
#ifndef _LOSC_SRC_LOCAL_OCCUPATION_H_
#define _LOSC_SRC_LOCAL_OCCUPATION_H_

#include "exception.h"
#include <matrix/matrix.h>
#include <memory>

namespace losc {

using matrix::Matrix;
using SharedMatrix = std::shared_ptr<Matrix>;

/**
 * @brief Calculate the local occupation matrix.
 *
 * @details The local occupation matrix is defined as
 * \f$\lambda_{ij} = \langle \phi_i|\rho|\phi_j \rangle,\f$
 * where \f$\lambda_{ij}\f$ is the local occupation matrix
 * element, \f$\phi_i\f$ and \f$\phi_j\f$ are the i-th and j-th localized
 * orbital and \f$\rho\f$ is the spin density operator. Note that \f$\lambda\f$,
 * \f$\phi\f$ and \f$\rho\f$ are for the same spin. In matrix form,
 * local occupation matrix \f$L\f$ is expressed as
 * \f$L = C_{LO} * S * D * S * C_{LO}^T\f$.
 *
 * @param [in] C_lo: LO coefficient matrix with dimension [nlo, nbasis].
 * See LocalizerBase::compute().
 * @param [in] S: AO overlap matrix with dimension [nbasis, nbasis].
 * @param [in] D: Spin density matrix under AO with dimension [nbasis, nbasis].
 * @return SharedMatrix: The local occupation matrix with dimension [nlo, nlo].
 */
inline SharedMatrix local_occupation_matrix(const Matrix &C_lo, const Matrix &S,
                                            const Matrix &D)
{
    size_t nlo = C_lo.row();
    size_t nbasis = C_lo.col();

    if (nlo > nbasis) {
        throw exception::DimensionError(
            "wrong dimension for LO coefficient matrix: number of LO is larger "
            "than the number of AO.");
    }
    if (!S.is_square() || nbasis != S.row()) {
        throw exception::DimensionError(
            S, nbasis, nbasis, "wrong dimension for AO overlap matrix.");
    }
    if (!D.is_square() || nbasis != D.row()) {
        throw exception::DimensionError(D, nbasis, nbasis,
                                        "wrong dimension for density matrix.");
    }

    Matrix SDS(nbasis, nbasis);
    SharedMatrix L = std::make_shared<Matrix>(nlo, nlo);
    matrix::mult_dgemm_ABAT(S, D, SDS);
    matrix::mult_dgemm_ABAT(C_lo, SDS, *L);
    return L;
}

} // namespace losc
#endif // _LOSC_SRC_LOCAL_OCCUPATION_H_
