#include "localization.h"
#include "eigen_helper.h"
#include "exception.h"

namespace losc {

void LocalizerBase::set_u_guess(RefMat U, const string &guess) const
{
    if (guess == "identity")
        U.setIdentity();
    // TODO
    else if (guess == "random") {
    } else if (guess == "random_fixed_seed") {
    } else {
        throw exception::LoscException(
            "Unknown initial guess of U matrix for localization.");
    }
}

void LocalizerBase::set_u_guess(RefMat U, ConstRefMat &U_guess,
                                double threshold) const
{
    if (!mtx_match_dimension(U, nlo_, nlo_)) {
        throw exception::DimensionError(
            U, nlo_, nlo_, "wrong dimension for localization U matrix.");
    }
    if (!U.isUnitary(threshold)) {
        throw exception::LoscException("Invalid U matrix: input U matrix "
                                       "is not unitary");
    }
    U = U_guess;
}

LocalizerBase::LocalizerBase(ConstRefMat &C_lo_basis)
    : nbasis_{C_lo_basis.rows()}, nlo_{C_lo_basis.cols()}, C_lo_basis_{
                                                               C_lo_basis}
{
    if (!mtx_match_dimension(C_lo_basis_, nbasis_, nlo_)) {
        throw exception::DimensionError(C_lo_basis_, nbasis_, nlo_,
                                        "LocalizerBase: wrong dimension for "
                                        "LO basis matrix under AO.");
    }
}

LocalizerBase::~LocalizerBase() {}

} // namespace losc
