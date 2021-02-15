#include "eigen_helper.hpp"
#include <losc/curvature.hpp>
#include <losc/exception.hpp>

namespace losc {

using exception::DimensionError;

CurvatureBase::CurvatureBase(const DFAInfo &dfa_info, ConstRefMat &df_pii,
                             ConstRefMat &df_Vpq_inverse, ConstRefMat &grid_lo,
                             ConstRefVec &grid_weight)
    : npts_{grid_weight.size()}, nlo_{df_pii.cols()},
      nfitbasis_{df_pii.rows()}, dfa_info_{dfa_info}, df_pii_{df_pii},
      df_Vpq_inverse_{df_Vpq_inverse}, grid_lo_{grid_lo}, grid_weight_{
                                                              grid_weight}
{
    if (!mtx_match_dimension(df_Vpq_inverse_, nfitbasis_, nfitbasis_)) {
        throw DimensionError(
            df_Vpq_inverse_, nfitbasis_, nfitbasis_,
            "wrong dimension for density fitting Vpq inverse matrix.");
    }
    if (!mtx_match_dimension(grid_lo_, npts_, nlo_)) {
        throw DimensionError(grid_lo_, npts_, nlo_,
                             "wrong dimension for grid value of LOs.");
    }
    if (grid_weight_.size() != npts_) {
        throw DimensionError(grid_weight_, npts_, 1,
                             "wrong dimension for grid weights.");
    }
}

namespace utils {

/**
 * @note This function will not track if the matrix `df_pii` is fully filled.
 * The user should be carefully traverse all the blocks of `df_pmn_block`.
 */
void convert_df_pmn2pii_blockwise(const vector<size_t> &p_index,
                                  ConstRefMat &df_pmn_block, ConstRefMat &C_lo,
                                  RefMat df_pii)
{
    // i, j: LO index.
    // p, q: fitbasis index.
    // m, n: AO basis index.

    // (p|ii): [nfitbasis, nlo].
    const size_t nfitbasis = df_pii.rows();
    const size_t nbasis = C_lo.rows();
    const size_t nlo = C_lo.cols();

    // check input dimension.
    if (df_pmn_block.cols() != nbasis * (nbasis + 1) / 2)
        throw exception::DimensionError("Wrong dimension of df_pmn_block.");
    if (df_pii.rows() != nfitbasis && df_pii.cols() != nlo)
        throw exception::DimensionError("Wrong dimension of df_pii.");
    if (p_index.size() != df_pmn_block.rows())
        throw exception::DimensionError(
            "number of rows in matrix df_pii dose not match with p index.");

    for (size_t pi = 0; pi < p_index.size(); ++pi) {
        const size_t p = p_index[pi];
        // For each p0 (p0|mn), the (m, n) block has dimension:
        // [nbasis x nbasis].
        LOSCMatrix df_pmn_p0_mn(nbasis, nbasis);
        for (size_t m = 0; m < nbasis; ++m) {
            for (size_t n = 0; n <= m; ++n) {
                const size_t mn = m * (m + 1) / 2 + n;
                df_pmn_p0_mn(m, n) = df_pmn_block(pi, mn);
            }
        }
        mtx_to_symmetric(df_pmn_p0_mn, "L");

        // For each p0 (p0|in), the (i, n) block dimension:
        // [nlo x nbasis].
        // Formula: (p0|in) = C_lo^T * (p|mn).
        LOSCMatrix df_pmn_p0_in(nlo, nbasis);
        df_pmn_p0_in.noalias() = C_lo.transpose() * df_pmn_p0_mn;

        // calculate element (p|ii)
        for (size_t i = 0; i < nlo; ++i) {
            df_pii(p, i) = df_pmn_p0_in.row(i).dot(C_lo.col(i));
        }
    }
}

} // namespace utils

} // namespace losc
