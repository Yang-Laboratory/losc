#ifndef __LOSC_SRC_EIGEN_HELPER_H
#define __LOSC_SRC_EIGEN_HELPER_H

#include <Eigen/Dense>
#include <string>

#include "exception.h"

namespace losc {

using ConstRefMat = const Eigen::Ref<const MatrixXd>;
using ConstRefVec = const Eigen::Ref<const VectorXd>;
using RefMat = Eigen::Ref<MatrixXd>;
using RefVec = Eigen::Ref<VectorXd>;

/**
 * @brief Check if the matrix is square or not.
 * @return bool
 */
inline bool mtx_is_square(ConstRefMat &m) { return m.rows() == m.cols(); }

/**
 * @brief Check if the matrix's dimension agrees with input expectation.
 * @return bool
 */
inline bool mtx_match_dimension(ConstRefMat &m, int row, int col)
{
    return (m.rows() == row) && (m.cols() == col));
}

/**
 * @brief Make the matrix to be symmetric.
 *
 * @param [in] uplo: when \p uplo equals to "U", the upper triangular part
 * is used. when \p uplo equals to "L", the lower triangular part is used.
 */
inline void mtx_to_symmetric(RefMat m, const string &uplo)
{
    const size_t row = m.rows();
    const size_t col = m.cols();
    if (!mtx_is_square(m)) {
        throw exception::DimensionError(
            "Cannot symmetrize a matrix that is not squared.");
    }
    if (uplo == "U") {
        for (size_t i = 0; i < row; i++) {
            for (size_t j = 0; j < i; j++) {
                m(i, j) = m(j, i);
            }
        }
    } else if (uplo == "L") {
        for (size_t i = 0; i < row; i++) {
            for (size_t j = 0; j < i; j++) {
                m(j, i) = m(i, j);
            }
        }
    } else {
        throw exception::LoscException(
            "Implementation Error: wrong choice to symmetrize a matrix.");
    }
}

} // namespace losc
#endif
