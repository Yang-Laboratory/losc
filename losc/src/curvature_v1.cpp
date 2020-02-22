/**
 * @file
 * @brief Losc curvature version 1 implementation.
 */
#include "curvature.h"

#include "blas_base.h"
#include "exception.h"
#include "matrix.h"
#include <cmath>
#include <cstdlib>

namespace losc {

shared_ptr<Matrix> CurvatureV1::compute_kappa_J()
{
    // i, j: LO index.
    // p, q: fitbasis index.
    // m, n: AO basis index.

    // (p|ii): [nfitbasis, nlo].
    Matrix df_pii(nfitbasis_, nlo_);
    const Matrix &C_lo = *C_lo_;
    const Matrix &df_pmn = *df_pmn_;
    for (size_t p = 0; p < nfitbasis_; ++p) {
        // For each p0 (p0|mn), the (m, n) block has dimension:
        // [nbasis x nbasis].
        Matrix df_pmn_p0_mn(nbasis_, nbasis_);
        for (size_t m = 0; m < nbasis_; ++m) {
            for (size_t n = 0; n <= m; ++n) {
                const size_t mn = m * (m + 1) / 2 + n;
                df_pmn_p0_mn(m, n) = df_pmn(p, mn);
            }
        }
        df_pmn_p0_mn.to_symmetric("L");

        // For each p0 (p0|in), the (i, n) block dimension:
        // [nlo x nbasis].
        // Formula: (p0|in) = C_lo^T * (p|mn).
        Matrix df_pmn_p0_in(nlo_, nbasis_);
        df_pmn_p0_in.noalias() = C_lo.transpose() * df_pmn_p0_mn;

        // calculate element (p|ii)
        for (size_t i = 0; i < nlo_; ++i) {
            df_pii(p, i) = df_pmn_p0_in.row(i).dot(C_lo.col(i));
        }
    }

    // kappa_J_ij = \sum_{pq} (\rho_i | p) V^{-1}_{pq} (q | \rho_j)
    // kappa_J = df_pii^T * V^{-1} * df_pii.
    auto kappa_J = std::make_shared<Matrix>(nlo_, nlo_);
    const Matrix &df_Vpq_inverse = *df_Vpq_inverse_;
    (*kappa_J).noalias() = df_pii.transpose() * df_Vpq_inverse * df_pii;

    return kappa_J;
}

shared_ptr<Matrix> CurvatureV1::compute_kappa_xc()
{
    auto kappa_xc = std::make_shared<Matrix>(nlo_, nlo_);
    kappa_xc->setZero();

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
        Matrix grid_lo_block(size, nlo_);
        // const Eigen::Block<Matrix, > &grid_ao_block =
        auto grid_ao_block =
            (*grid_basis_value_).block(n * block_size, 0, size, nbasis_);
        grid_lo_block.noalias() = grid_ao_block * (*C_lo_);

        // calculate LO grid value to the power of 4/3.
        grid_lo_block = grid_lo_block.array().square().matrix();
        grid_lo_block = grid_lo_block.array().pow(2.0 / 3.0).matrix();

        // sum over current block contribution to kappa xc.
        vector<double>::const_iterator wt =
            grid_weight_->begin() + n * block_size;
        for (size_t ip = 0; ip < size; ++ip) {
            const double wt_value = wt[ip];
            for (size_t i = 0; i < nlo_; ++i) {
                for (size_t j = 0; j <= i; ++j) {
                    const double pi = grid_lo_block(ip, i);
                    const double pj = grid_lo_block(ip, j);
                    (*kappa_xc)(i, j) += wt_value * pi * pj;
                }
            }
        }
    }
    kappa_xc->to_symmetric("L");

    return kappa_xc;
}

shared_ptr<Matrix> CurvatureV1::compute()
{
    shared_ptr<Matrix> p_kappa_J = compute_kappa_J();
    shared_ptr<Matrix> p_kappa_xc = compute_kappa_xc();

    // combine J and xc
    auto kappa = std::make_shared<Matrix>(nlo_, nlo_);
    const double xc_factor =
        -para_tau_ * para_cx_ * 2.0 / 3.0 * (1.0 - para_alpha_);
    const double j_factor = 1.0 - para_alpha_ - para_beta_;

    (*kappa).noalias() = j_factor * (*p_kappa_J) + xc_factor * (*p_kappa_xc);

    return kappa;
}

} // namespace losc
