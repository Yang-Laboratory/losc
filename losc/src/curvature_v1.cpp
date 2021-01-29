/**
 * @file
 * @brief LOSC curvature version 1 implementation.
 */
#include "curvature.h"
#include "eigen_helper.h"
#include "exception.h"
#include <cmath>
#include <cstdlib>

namespace losc {

MatrixXd CurvatureV1::compute_kappa_J() const
{
    // K_ij = <ii|p> V_{pq}^-1 <q|jj>
    return df_pii_.transpose() * df_Vpq_inverse_ * df_pii_;
}

MatrixXd CurvatureV1::compute_kappa_xc() const
{
    MatrixXd kappa_xc(nlo_, nlo_);
    kappa_xc.setZero();

    // The LO grid value matrix has dimension of (npts, nlo), which could be
    // very large. To limit the memory usage on LO grid value, we build it by
    // blocks. For now we limit the memory usage to be 1GB.
    const size_t block_size =
        1000ULL * 1000ULL * 1000ULL / sizeof(double) / nlo_;
    size_t nBLK = npts_ / block_size;
    const size_t res = npts_ % block_size;
    if (res != 0) {
        nBLK += 1;
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    // loop over all the blocks of grid.
    for (size_t n = 0; n < nBLK; ++n) {
        size_t size = block_size;
        if (n == nBLK - 1 && res != 0)
            size = res;
        // Calculate LO grid value for each block.
        // Formula: grid_lo = grid_ao * C_lo
        // grid_lo: [block_size, nlo]
        // grid_ao: [block_size, nbasis]
        // C_lo: [nbasis, nlo]
        MatrixXd grid_lo_block(size, nlo_);
        grid_lo_block.noalias() =
            grid_basis_value_.block(n * block_size, 0, size, nbasis_) * C_lo_;

        // calculate LO grid value to the power of 4/3.
        grid_lo_block = grid_lo_block.array().square().matrix();
        grid_lo_block = grid_lo_block.array().pow(2.0 / 3.0).matrix();

        // sum over current block contribution to kappa xc.
        auto wt = grid_weight_.data() + n * block_size;
        for (size_t ip = 0; ip < size; ++ip) {
            const double wt_value = wt[ip];
            for (size_t i = 0; i < nlo_; ++i) {
                for (size_t j = 0; j <= i; ++j) {
                    const double pi = grid_lo_block(ip, i);
                    const double pj = grid_lo_block(ip, j);
                    kappa_xc(i, j) += wt_value * pi * pj;
                }
            }
        }
    }
    mtx_to_symmetric(kappa_xc, "L");

    return kappa_xc;
}

void CurvatureV1::kappa(RefMat K) const
{
    if (!mtx_match_dimension(K, nlo_, nlo_)) {
        throw exception::DimensionError(
            K, nlo_, nlo_,
            "CurvatureV1::kappa(): wrong dimension of the input kappa matrix.");
    }
    MatrixXd kappa_J = compute_kappa_J();
    MatrixXd kappa_xc = compute_kappa_xc();
    K.noalias() = (1 - dfa_info_.hf_x()) * kappa_J -
                  dfa_info_.gga_x() * (2.0 * tau_ * cx_ / 3.0) * kappa_xc;
}

} // namespace losc
