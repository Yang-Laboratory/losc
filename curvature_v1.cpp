#include "curvature.h"
#include <matrix/matrix.h>
#include <cmath>
#include <cstdlib>

namespace losc {

CurvatureV1::CurvatureV1(SharedMatrix C_lo, SharedMatrix df_pmn, SharedMatrix df_Vpq_inverse)
    : C_lo_{C_lo}, df_pmn_{df_pmn}, df_Vpq_inverse_{df_Vpq_inverse},
    nlo_{C_lo_->row()}, nbasis_{C_lo_->col()}, nfitbasis_{df_Vpq_inverse_->row()}
{
    if (df_pmn_->col() != nbasis_ * (nbasis_ + 1) / 2) {
        printf("Dimension error: df_mnp matrix.\n");
        std::exit(EXIT_FAILURE);
    }
    if (!df_Vpq_inverse_->is_square()) {
        printf("Dimension error: df_Vpq_inverse matrix.\n");
        std::exit(EXIT_FAILURE);
    }
}

void CurvatureV1::compute_kappa_J()
{
    // i, j: LO index.
    // p, q: fitbasis index.
    // m, n: AO basis index.

    // (p|ii): [nfitbasis, nlo].
    SharedMatrix df_pii = std::make_shared<Matrix> (nlo_, nfitbasis_);
    for (size_t p = 0; p < nfitbasis_; ++p) {
        // at each p: (p|mn).
        // dimension: nbasis x nbasis.
        SharedMatrix df_pmn_p_mn = std::make_shared<Matrix> (nbasis_, nbasis_);
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
        matrix::mult_dgemm(1.0, *C_lo_, "N", *df_pmn_p_mn, "N", 0.0, df_pmn_p_in);

        // calculate element (p|ii)
        for (size_t i = 0; i < nlo_; ++i) {
            int nbasis = nbasis_;
            (*df_pii)(p, i) = matrix::blas::ddot_(&nbasis, df_pmn_p_in.data() + i * nbasis,
                                                  matrix::blas::ione,
                                                  C_lo_->data() + i * nbasis,
                                                  matrix::blas::ione);
        }
    }

    // kappa_J_ij = \sum_{pq} (\rho_i | p) V^{-1}_{pq} (q | \rho_j)
    // kappa_J = df_pii^T * V^{-1} * df_pii.
    kappa_J_ = std::make_shared<Matrix> (nlo_, nlo_);
    SharedMatrix V_pii = std::make_shared<Matrix> (nfitbasis_, nlo_);
    matrix::mult_dgemm(1.0, *df_Vpq_inverse_, "N", *df_pii, "N", 0.0, *V_pii);
    matrix::mult_dgemm(1.0, *df_pii, "T", *V_pii, "N", 0.0, *kappa_J_);
}

}
