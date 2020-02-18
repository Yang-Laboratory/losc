/**
 * @file
 * @brief Losc curvature version 1 implementation.
 */
#include "curvature.h"

#include "blas_base.h"
#include "exception.h"
#include <cmath>
#include <cstdlib>
#include <matrix/matrix.h>

namespace losc {

SharedMatrix CurvatureV1::compute_kappa_J()
{
    // i, j: LO index.
    // p, q: fitbasis index.
    // m, n: AO basis index.

    // (p|ii): [nfitbasis, nlo].
    SharedMatrix df_pii = std::make_shared<Matrix>(nfitbasis_, nlo_);
    for (size_t p = 0; p < nfitbasis_; ++p) {
        // at each p: (p|mn).
        // dimension: nbasis x nbasis.
        SharedMatrix df_pmn_p_mn = std::make_shared<Matrix>(nbasis_, nbasis_);
        for (size_t m = 0; m < nbasis_; ++m) {
            for (size_t n = 0; n <= m; ++n) {
                const size_t mn = m * (m + 1) / 2 + n;
                (*df_pmn_p_mn)(m, n) = (*df_pmn_)(p, mn);
            }
        }
        df_pmn_p_mn->to_symmetric("L");

        // at each p: (p|in)
        // dimension: nlo x nbasis.
        // (p|in) = C_lo * (p|mn).
        Matrix df_pmn_p_in(nlo_, nbasis_);
        matrix::mult_dgemm(1.0, *C_lo_, "N", *df_pmn_p_mn, "N", 0.0,
                           df_pmn_p_in);

        // calculate element (p|ii)
        for (size_t i = 0; i < nlo_; ++i) {
            int nbasis = nbasis_;
            (*df_pii)(p, i) = losc::blas::ddot_(
                &nbasis, df_pmn_p_in.data() + i * nbasis, losc::blas::ione,
                C_lo_->data() + i * nbasis, losc::blas::ione);
        }
    }

    // kappa_J_ij = \sum_{pq} (\rho_i | p) V^{-1}_{pq} (q | \rho_j)
    // kappa_J = df_pii^T * V^{-1} * df_pii.
    auto kappa_J = std::make_shared<Matrix>(nlo_, nlo_);
    SharedMatrix V_pii = std::make_shared<Matrix>(nfitbasis_, nlo_);
    matrix::mult_dgemm(1.0, *df_Vpq_inverse_, "N", *df_pii, "N", 0.0, *V_pii);
    matrix::mult_dgemm(1.0, *df_pii, "T", *V_pii, "N", 0.0, *kappa_J);

    return kappa_J;
}

SharedMatrix CurvatureV1::compute_kappa_xc()
{
    auto kappa_xc = std::make_shared<Matrix>(nlo_, nlo_);

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
        // calculate LO grid value for each block, which
        // has dimension of (size, nlo).
        auto grid_lo_block = std::make_shared<Matrix>(size, nlo_);
        auto grid_basis_block = std::make_shared<Matrix>(
            size, nbasis_, grid_basis_value_->data() + n * block_size * nbasis_,
            matrix::Matrix::kShallowCopy);
        matrix::mult_dgemm(1.0, *grid_basis_block, "N", *C_lo_, "T", 0.0,
                           *grid_lo_block);
        grid_basis_block.reset();

        // calculate LO grid value to the power of 2/3.
        for (size_t i = 0; i < grid_lo_block->size(); ++i) {
            grid_lo_block->data()[i] *= grid_lo_block->data()[i];
            grid_lo_block->data()[i] =
                std::pow(grid_lo_block->data()[i], 2.0 / 3.0);
        }

        // sum over current block contribution to kappa xc.
        const double *wt = grid_weight_->data() + n * block_size;
        for (size_t ip = 0; ip < size; ++ip) {
            const double wt_value = wt[ip];
            for (size_t i = 0; i < nlo_; ++i) {
                for (size_t j = 0; j <= i; ++j) {
                    const double pi = (*grid_lo_block)(ip, i);
                    const double pj = (*grid_lo_block)(ip, j);
                    (*kappa_xc)(i, j) += (*grid_weight_)[ip] * pi * pj;
                }
            }
        }
    }
    kappa_xc->to_symmetric("L");

    return kappa_xc;
}

SharedMatrix CurvatureV1::compute()
{
    SharedMatrix kappa_J = compute_kappa_J();
    SharedMatrix kappa_xc = compute_kappa_xc();

    // combine J and xc
    auto kappa = std::make_shared<Matrix>(nlo_, nlo_);
    const double xc_factor =
        -para_tau_ * para_cx_ * 2.0 / 3.0 * (1.0 - para_alpha_);
    const double j_factor = 1.0 - para_alpha_ - para_beta_;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (size_t i = 0; i < kappa->size(); ++i) {
        kappa->data()[i] =
            j_factor * kappa_J->data()[i] + xc_factor * kappa_xc->data()[i];
    }

    return kappa;
}

} // namespace losc
