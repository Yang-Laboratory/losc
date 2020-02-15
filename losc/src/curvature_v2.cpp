#include <cmath>

#include "curvature.h"
#include "exception.h"

namespace losc {

SharedMatrix CurvatureV2::compute()
{
    // construct absolute overlap under LO.
    auto S_lo = std::make_shared<Matrix>(nlo_, nlo_);
    // The LO grid value matrix has dimension of (npts, nlo), which could be very large.
    // To limit the memory usage on LO grid value, we build it by
    // blocks. For now we limit the memory usage to be 1GB.
    const size_t block_size = 1000ULL * 1000ULL * 1000ULL / sizeof(double) / nlo_;
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
        auto grid_basis_block =
            std::make_shared<Matrix>(size, nbasis_,
                                     grid_basis_value_->data() + n * block_size * nbasis_,
                                     matrix::Matrix::kShallowCopy);
        matrix::mult_dgemm(1.0, *grid_basis_block, "N", *C_lo_, "T", 0.0, *grid_lo_block);
        grid_basis_block.reset();
        // sum over current block for the LO overlap.
        const double *wt = grid_weight_->data() + n * block_size;
        for (size_t ip = 0; ip < size; ++ip) {
            const double wt_value = wt[ip];
            for (size_t i = 0; i < nlo_; ++i) {
                for (size_t j = 0; j <= i; ++j) {
                    const double pi = (*grid_lo_block)(ip, i);
                    const double pj = (*grid_lo_block)(ip, j);
                    (*S_lo)(i, j) += wt_value * std::fabs(pi * pj);
                }
            }
        }
    }
    S_lo->to_symmetric("L");

    // build the curvature version 1.
    CurvatureV1 kappa1_man(dfa_type_, C_lo_, df_pmn_, df_Vpq_inverse_, grid_basis_value_, grid_weight_);
    auto kappa1 = kappa1_man.compute();

    // build LOSC2 kappa matrix:
    // K2[ij] = erf(tau * S[ij]) * sqrt(abs(K1[ii] * K1[jj])) + erfc(tau * S[ij]) * K][ij]
    using std::erf;
    using std::erfc;
    using std::fabs;
    using std::sqrt;
    auto kappa2 = std::make_shared<Matrix>(nlo_, nlo_);
    for (size_t i = 0; i < nlo_; ++i) {
        const double K1_ii = (*kappa1)(i, i);
        (*kappa2)(i, i) = K1_ii;
        for (size_t j = 0; j < i; ++j) {
            const double S_ij = (*S_lo)(i, j);
            const double K1_ij = (*kappa1)(i, j);
            const double K1_jj = (*kappa1)(j, j);
            const double f = para_zeta_ * S_ij;
            (*kappa2)(i, j) = erf(f) * sqrt(fabs(K1_ii * K1_jj)) + erfc(f) * K1_ij;
        }
    }
    kappa2->to_symmetric("L");

    return kappa2;
}


}
