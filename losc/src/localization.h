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
    size_t nlo_;    /* number of LO. */
    size_t nbasis_; /* number of AO basis. */

    enum PrintLevel print_level_ = kPrintLevelNo;

    /**
     * Unitary matrix that transfer LO basis coefficient matrix into LO
     * coefficient matrix.
     * Dimension: [nlo, nlo]
     * $C_{lo} = U * C_{lo_basis}$.
     */
    SharedMatrix U_;

    /**
     * LO basis coefficient matrix under AO.
     * Dimension: [nlo, nbasis]
     */
    SharedMatrix C_lo_basis_;

    public:
    /**
     * @ param [in] C_lo_basis: LO basis coefficient matrix under Ao.
     *
     * The U matrix will be intialized as an identity matrix.
     */
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
     * @ return: the U matrix.
     */
    SharedMatrix get_u() const { return U_; }

    /**
     * Set up the U matrix.
     *
     * @ param [in] threshold: the threshold to check if the input matrix is unitary.
     */
    void set_u_matrix(const SharedMatrix& U, double threshold = 1e-8)
    {
        if (U) {
            throw exception::LoscException("invalid U matrix: input U matrix is null.");
        }
        if (!U->is_square() && U->row() != nlo_) {
            throw exception::DimensionError(*U, nlo_, nlo_, "wrong dimension for localization U matrix.");
        }
        Matrix UU(nlo_, nlo_);
        matrix::mult_dgemm(1.0, *U, "N", *U, "T", 0.0, UU);
        if (UU.is_identity(threshold)) {
            throw exception::LoscException("Invalid U matrix: input U matrix is not unitary, U * U^T != I.");
        }
        matrix::mult_dgemm(1.0, *U, "T", *U, "N", 0.0, UU);
        if (UU.is_identity(threshold)) {
            throw exception::LoscException("Invalid U matrix: input U matrix is not unitary, U^T * U != I.");
        }
        U_ = U;
    }

    /**
     * Set up the print level.
     */
    void set_print(enum PrintLevel level) {print_level_ = level;}

    /**
     * Compute the LO coefficient matrix under AO.
     *
     * The current U matrix stored in `U_` will be used as the initial guess.
     * After calling this function, `U_` matrix is updated.
     *
     * @ return: the LO coefficient matrix under AO.
     */
    virtual SharedMatrix compute() = 0;
};

/**
 * Losc localization version 2. This is the version used in the Losc2
 * paper (xxx).
 */
class LoscLocalizerV2 : public LocalizerBase {
    private:
    /**
     * Hamiltonian matrix under AO used in Losc2 localization.
     * In Losc2 method, it is the DFA Hamiltonian.
     * Dimension: [nbasis, nbasis].
     */
    SharedMatrix H_ao_;

    /**
     * Dipole matrix under AO in order of x, y and z directions.
     * Dimension: [nbasis, nbasis].
     */
    vector<SharedMatrix> Dipole_ao_;

    /**
     * Maximum iteration number of Jacobi-sweep algorithm.
     * Default: 1000
     */
    size_t js_max_iter_ = 1000;

    /**
     * If use random permutation in Jacobi-sweep algorithm or not.
     * Default: true
     */
    bool js_random_permutation_ = true;

    /**
     * Convergence tolerance tolerance in Jacobi-sweep algorithm.
     * Default: 1e-10
     */
    double js_tol_ = 1e-10;

    /**
     * Parameter $C$ in Losc2 localization cost function.
     * See Eq. (xxx) in Losc2 paper.
     * Default: 1000
     */
    double para_c_ = 1000;

    /**
     * Parameter $\gamma$ in Losc2 localization cost function.
     * See Eq. (xxx) in Losc2 paper.
     * Default: 0.707
     */
    double para_gamma_ = 0.707;

    void message(std::string t, ...) const;

    public:
    /**
     * @ param [in] C_lo_basis: LO basis coefficient matrix under AO.
     * @ param [in] H_ao: Hamiltonian matrix under AO used in localization.
     * @ param [in] Dipole_ao: dipole matrix under AO.
     */
    LoscLocalizerV2(const SharedMatrix& C_lo_basis, const SharedMatrix& H_ao, const vector<SharedMatrix>& Dipole_ao);

    void set_tolerance(double tol) {js_tol_ = tol;}
    void set_max_iteration(size_t max_iter) {js_max_iter_ = max_iter;}
    void set_random_permutation(bool true_or_false) {js_random_permutation_ = true_or_false;}

    /**
     * Compute the LO coefficient matrix under AO.
     *
     * The current U matrix stored in `U_` will be used as the initial guess.
     * After calling this function, `U_` matrix is updated.
     *
     * @ return: the LO coefficient matrix under AO.
     */
    SharedMatrix compute() override;
};

}

#endif //_LOSC_LOCALIZATION_H_
