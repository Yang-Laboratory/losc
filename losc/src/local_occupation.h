/**
 * @file
 * @brief Losc localization occupation matrix related declarations.
 */
#ifndef _LOSC_SRC_LOCAL_OCCUPATION_H_
#define _LOSC_SRC_LOCAL_OCCUPATION_H_

#include "exception.h"
#include "matrix.h"
#include <memory>

namespace losc {

using losc::Matrix;
using std::shared_ptr;

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
 * \f$L = C_{LO}^T * S * D * S * C_{LO}\f$.
 *
 * @param [in] C_lo: LO coefficient matrix with dimension [nbasis, nlo].
 * See LocalizerBase::compute().
 * @param [in] S: AO overlap matrix with dimension [nbasis, nbasis].
 * @param [in] D: Spin density matrix under AO with dimension [nbasis, nbasis].
 * @return shared_ptr<Matrix>: The local occupation matrix with dimension [nlo,
 * nlo].
 */
inline shared_ptr<Matrix>
local_occupation_matrix(const Matrix &C_lo, const Matrix &S, const Matrix &D)
{
    size_t nlo = C_lo.cols();
    size_t nbasis = C_lo.rows();

    if (nlo > nbasis) {
        throw exception::DimensionError(
            "wrong dimension for LO coefficient matrix: number of LO is larger "
            "than the number of AO.");
    }
    if (!S.is_square() || nbasis != S.rows()) {
        throw exception::DimensionError(
            S, nbasis, nbasis, "wrong dimension for AO overlap matrix.");
    }
    if (!D.is_square() || nbasis != D.rows()) {
        throw exception::DimensionError(D, nbasis, nbasis,
                                        "wrong dimension for density matrix.");
    }

    shared_ptr<Matrix> L = std::make_shared<Matrix>(nlo, nlo);
    (*L).noalias() = C_lo.transpose() * S * D * S * C_lo;
    return L;
}

} // namespace losc
#endif // _LOSC_SRC_LOCAL_OCCUPATION_H_
