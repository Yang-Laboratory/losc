/**
 * @file
 * @brief definition relates to Losc curvature.
 */
#include "curvature.h"

#include "exception.h"

namespace losc {

using exception::DimensionError;

CurvatureBase::CurvatureBase(enum DFAType dfa, const SharedMatrix &C_lo,
                             const SharedMatrix &df_pmn,
                             const SharedMatrix &df_Vpq_inverse,
                             const SharedMatrix &grid_basis_value,
                             const SharedDoubleVector &grid_weight)
    : dfa_type_{dfa}, npts_{grid_weight->size()}, nlo_{C_lo->row()},
      nbasis_{C_lo->col()}, nfitbasis_{df_pmn->row()}, C_lo_{C_lo},
      df_pmn_{df_pmn}, df_Vpq_inverse_{df_Vpq_inverse},
      grid_basis_value_{grid_basis_value}, grid_weight_{grid_weight}
{
    if (df_pmn_->col() != nbasis_ * (nbasis_ + 1) / 2) {
        throw DimensionError(*df_pmn_, nfitbasis_, nbasis_,
                             "wrong dimension for density fitting three-body "
                             "integral matrix <p|mn>.");
    }
    if (!df_Vpq_inverse_->is_square() || df_Vpq_inverse_->row() != nfitbasis_) {
        throw DimensionError(
            *df_Vpq_inverse_, nfitbasis_, nfitbasis_,
            "wrong dimension for density fitting Vpq inverse matrix.");
    }
    if (npts_ != grid_basis_value_->row() ||
        nbasis_ != grid_basis_value_->col()) {
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
