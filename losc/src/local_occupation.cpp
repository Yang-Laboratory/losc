#include "local_occupation.h"
#include "eigen_helper.h"
#include "exception.h"

namespace losc {

MatrixXd local_occupation(ConstRefMat &C_lo, ConstRefMat &S, ConstRefMat &D)
{
    const size_t nlo = C_lo.cols();
    const size_t nbasis = C_lo.rows();

    if (nlo > nbasis) {
        throw exception::DimensionError(
            "wrong dimension for LO coefficient matrix: number of LO is larger "
            "than the number of AO.");
    }
    if (!mtx_match_dimension(S, nbasis, nbasis)) {
        throw exception::DimensionError(
            S, nbasis, nbasis, "wrong dimension for AO overlap matrix.");
    }
    if (!mtx_match_dimension(D, nbasis, nbasis)) {
        throw exception::DimensionError(D, nbasis, nbasis,
                                        "wrong dimension for density matrix.");
    }

    return C_lo.transpose() * S * D * S * C_lo;
}

} // namespace losc
