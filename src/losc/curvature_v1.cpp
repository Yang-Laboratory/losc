#include "eigen_helper.hpp"
#include <cmath>
#include <cstdlib>
#include <losc/curvature.hpp>
#include <losc/exception.hpp>

namespace losc {

LOSCMatrix CurvatureV1::compute_kappa_J() const
{
    // K_ij = <ii|p> V_{pq}^-1 <q|jj>
    return df_pii_.transpose() * df_Vpq_inverse_ * df_pii_;
}

LOSCMatrix CurvatureV1::compute_kappa_xc() const
{
    LOSCMatrix kappa_xc(nlo_, nlo_);
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

        const auto grid_lo_block =
            grid_lo_.block(n * block_size, 0, size, nlo_);
        const auto wt = grid_weight_.segment(n * block_size, size);

        // sum over current block contribution to kappa xc.
        for (size_t ip = 0; ip < size; ++ip) {
            const double wt_value = wt[ip];
            for (size_t i = 0; i < nlo_; ++i) {
                for (size_t j = 0; j <= i; ++j) {
                    const double pi = grid_lo_block(ip, i);
                    const double pj = grid_lo_block(ip, j);
                    kappa_xc(i, j) += wt_value * std::pow(pi * pi, 2.0 / 3.0) *
                                      std::pow(pj * pj, 2.0 / 3.0);
                }
            }
        }
    }
    mtx_to_symmetric(kappa_xc, "L");

    return kappa_xc;
}

void CurvatureV1::C_API_kappa(RefMat K) const
{
    if (!mtx_match_dimension(K, nlo_, nlo_)) {
        throw exception::DimensionError(
            K, nlo_, nlo_,
            "CurvatureV1::kappa(): wrong dimension of the input kappa matrix.");
    }
    LOSCMatrix kappa_J = compute_kappa_J();
    LOSCMatrix kappa_xc = compute_kappa_xc();
    K.noalias() = (1 - dfa_info_.hf_x()) * kappa_J -
                  dfa_info_.gga_x() * (2.0 * tau_ * cx_ / 3.0) * kappa_xc;
}

} // namespace losc
