/**
 * @file
 * @brief definition relates to Losc curvature.
 */
#include "curvature.h"

#include "exception.h"
#include "matrix.h"

namespace losc {

using exception::DimensionError;

CurvatureBase::CurvatureBase(DFAType dfa, const shared_ptr<Matrix> &C_lo,
                             const shared_ptr<Matrix> &df_pmn,
                             const shared_ptr<Matrix> &df_Vpq_inverse,
                             const shared_ptr<Matrix> &grid_basis_value,
                             const shared_ptr<vector<double>> &grid_weight)
    : dfa_type_{dfa}, npts_{grid_weight->size()}, nlo_{C_lo->cols()},
      nbasis_{C_lo->rows()}, nfitbasis_{df_pmn->rows()}, C_lo_{C_lo},
      df_pmn_{df_pmn}, df_Vpq_inverse_{df_Vpq_inverse},
      grid_basis_value_{grid_basis_value}, grid_weight_{grid_weight}
{
    if (df_pmn_->cols() != nbasis_ * (nbasis_ + 1) / 2) {
        throw DimensionError(*df_pmn_, nfitbasis_, nbasis_,
                             "wrong dimension for density fitting three-body "
                             "integral matrix <p|mn>.");
    }
    if (!df_Vpq_inverse_->is_square() ||
        df_Vpq_inverse_->rows() != nfitbasis_) {
        throw DimensionError(
            *df_Vpq_inverse_, nfitbasis_, nfitbasis_,
            "wrong dimension for density fitting Vpq inverse matrix.");
    }
    if (npts_ != grid_basis_value_->rows() ||
        nbasis_ != grid_basis_value_->cols()) {
        throw DimensionError(*grid_basis_value_, npts_, nbasis_,
                             "wrong dimension for grid value of AO basis.");
    }
    switch (dfa) {
    case losc::LDA: {
        para_alpha_ = para_beta_ = 0.0;
        break;
    }
    case losc::GGA: {
        para_alpha_ = para_beta_ = 0.0;
        break;
    }
    case losc::B3LYP: {
        para_alpha_ = 0.2;
        para_beta_ = 0.0;
        break;
    }
    default: {
        throw exception::LoscException("Unknown DFA choice.");
    }
    }
}

} // namespace losc
