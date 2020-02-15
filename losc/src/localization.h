/**
 * file: localization.h
 */
#ifndef _LOSC_LOCALIZATION_H_
#define _LOSC_LOCALIZATION_H_

#include <vector>
#include <iostream>
#include <cstddef>
#include <memory>
#include <string>
#include <matrix/matrix.h>

#include "exception.h"

namespace losc {

using matrix::Matrix;
using std::vector;
using std::shared_ptr;

using SharedMatrix = shared_ptr<Matrix>;

enum PrintLevel {
    kPrintLevelNo,
    kPrintLevelNormal,
    kPrintLevelDebug1,
    kPrintLevelDebug2,
};

class LocalizerBase {
    protected:
    size_t nlo_;    // number of LO.
    size_t nbasis_; // number of AO basis.

    enum PrintLevel print_level_ = kPrintLevelNo;

    /**
     * Unitary matrix that transfer LO basis coefficient matrix into LO coefficient matrix.
     * Dimension: nlo x nlo.
     * phi_i = \sum_j U_{ij} psi_j, where phi_i is the i-th LO and psi_j is the j-th LO basis.
     */
    SharedMatrix U_;

    /**
     * LO basis coefficient matrix under AO.
     * Dimension: nlo x nbasis.
     */
    SharedMatrix C_lo_basis_;

    public:
    LocalizerBase(const SharedMatrix& C_lo_basis)
        : nlo_{C_lo_basis->row()}, nbasis_{C_lo_basis->col()},
        print_level_{kPrintLevelNo}, C_lo_basis_{C_lo_basis}
    {
        if (nlo_ > nbasis_) {
            throw exception::DimensionError("wrong dimension for LO coefficient matrix: number of LO is larger than the number of AO.");
        }
        U_ = std::make_shared<Matrix> (nlo_, nlo_);
        U_->set_identity();
    }

    /**
     * get the U matrix.
     */
    SharedMatrix get_u() const { return U_; }

    /**
     * Set up the initial U matrix.
     */
    void set_initial_u_matrix(const SharedMatrix& U)
    {
        if (U) {
            throw exception::LoscException("invalid U matrix: input U matrix is null.");
        }
        if (!U->is_square() && U->row() != nlo_) {
            throw exception::DimensionError(*U, nlo_, nlo_, "wrong dimension for localization U matrix.");
        }
        U_ = U;
    }

    /**
     * Set up the print level.
     */
    void set_print(enum PrintLevel level) {print_level_ = level;}

    /**
     * do localization and compute the LO coefficient matrix.
     *
     * Calling this function will calculate the `C_lo_` and `U_`.
     * @ return: LO coefficient matrix under AO with dimension [nlo, nbasis].
     */
    virtual SharedMatrix compute() = 0;
};

/**
 * Losc2 Localization
 */
class Losc2Localizer : public LocalizerBase {
    private:
    SharedMatrix H_ao_;
    vector<SharedMatrix> Dipole_ao_;

    bool js_random_permutation_ = true;
    size_t js_max_iter_ = 1000;
    double js_tol_ = 1e-10;
    double para_c_ = 1000;
    double para_gamma_ = 0.707;

    void message(std::string t, ...) const;

    public:
    Losc2Localizer(const SharedMatrix& C_lo_basis, const SharedMatrix& H_ao, const vector<SharedMatrix>& Dipole_ao);

    void set_tolerance(double tol) {js_tol_ = tol;}
    void set_max_iteration(size_t max_iter) {js_max_iter_ = max_iter;}
    void set_random_permutation(bool true_or_false) {js_random_permutation_ = true_or_false;}

    SharedMatrix compute() override;
};

}

#endif //_LOSC_LOCALIZATION_H_
