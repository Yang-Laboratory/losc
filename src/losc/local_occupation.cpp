#include "eigen_helper.hpp"
#include <losc/exception.hpp>
#include <losc/local_occupation.hpp>
#include <utility>

namespace losc {

MatrixXd local_occupation(ConstRefMat &C_lo, ConstRefMat &S, ConstRefMat &D)
{
    const size_t nlo = C_lo.cols();
    MatrixXd local_occ(nlo, nlo);
    C_API_local_occupation(C_lo, S, D, local_occ);
    return std::move(local_occ);
}

void C_API_local_occupation(ConstRefMat &C_lo, ConstRefMat &S, ConstRefMat &D,
                            RefMat LocalOcc)
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
    if (!mtx_match_dimension(LocalOcc, nlo, nlo)) {
        throw exception::DimensionError(
            LocalOcc, nlo, nlo, "wrong dimension for local occupation matrix.");
    }

    LocalOcc = C_lo.transpose() * S * D * S * C_lo;
}

} // namespace losc
