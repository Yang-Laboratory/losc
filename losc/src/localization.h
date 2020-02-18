/**
 * @file
 * @brief Declaration relates to orbital localization.
 */
#ifndef _LOSC_SRC_LOCALIZATION_H_
#define _LOSC_SRC_LOCALIZATION_H_

#include "exception.h"
#include <cstddef>
#include <iostream>
#include <matrix/matrix.h>
#include <memory>
#include <string>
#include <vector>

namespace losc {

using matrix::Matrix;
using std::shared_ptr;
using std::vector;

using SharedMatrix = shared_ptr<Matrix>;

enum PrintLevel {
    kPrintLevelNo,
    kPrintLevelNormal,
    kPrintLevelDebug1,
    kPrintLevelDebug2,
};

/**
 * @brief Base class for Losc localization.
 */
class LocalizerBase {
  protected:
    size_t nlo_;    /**<number of LO. */
    size_t nbasis_; /**< number of AO basis. */

    enum PrintLevel print_level_ = kPrintLevelNo;

    /**
     * @brief The unitary matrix that transfer LO basis coefficient matrix into
     * LO coefficient matrix.
     * @details Dimension: [nlo, nlo]. Relation: C_lo = U * C_lo_basis
     */
    SharedMatrix U_;

    /**
     * @brief LO basis coefficient matrix under AO.
     * @details
     * Dimension: [nlo, nbasis]
     *
     * Relation between LO basis and AO via the coefficient matrix \f$C\f$ is
     * \f[
     * \psi_i = \sum_\mu C_{i\mu} \phi_{\mu},
     * \f]
     * where \f$C\f$ matrix is stored in `C_lo_basis_`.
     *
     * @par Availability
     * You have to build the matrix by yourself. Usually, Losc localization
     * use the CO from DFA as the LO basis.
     */
    SharedMatrix C_lo_basis_;

  public:
    /**
     * @brief LocalizerBase class constructor
     * @details After creation, the matrix (LocalizerBase::U_) will be
     * intialized as an identity matrix.
     * @param [in] C_lo_basis: LO basis coefficient matrix under Ao. See
     * LocalizerBase::C_lo_basis_.
     */
    LocalizerBase(const SharedMatrix &C_lo_basis)
        : nlo_{C_lo_basis->row()}, nbasis_{C_lo_basis->col()},
          print_level_{kPrintLevelNo}, C_lo_basis_{C_lo_basis}
    {
        if (nlo_ > nbasis_) {
            throw exception::DimensionError(
                "wrong dimension for LO coefficient matrix: number of LO is "
                "larger than the number of AO.");
        }
        U_ = std::make_shared<Matrix>(nlo_, nlo_);
        U_->set_identity();
    }

    /**
     * @return SharedMatrix: the U matrix.
     */
    SharedMatrix get_u() const { return U_; }

    /**
     * @brief Set up the U matrix.
     * @param [in] U: The localization U matrix.
     * @param [in] threshold: the threshold to check if the input matrix is
     * unitary.
     */
    void set_u_matrix(const SharedMatrix &U, double threshold = 1e-8)
    {
        if (U) {
            throw exception::LoscException(
                "invalid U matrix: input U matrix is null.");
        }
        if (!U->is_square() && U->row() != nlo_) {
            throw exception::DimensionError(
                *U, nlo_, nlo_, "wrong dimension for localization U matrix.");
        }
        Matrix UU(nlo_, nlo_);
        matrix::mult_dgemm(1.0, *U, "N", *U, "T", 0.0, UU);
        if (UU.is_identity(threshold)) {
            throw exception::LoscException("Invalid U matrix: input U matrix "
                                           "is not unitary, U * U^T != I.");
        }
        matrix::mult_dgemm(1.0, *U, "T", *U, "N", 0.0, UU);
        if (UU.is_identity(threshold)) {
            throw exception::LoscException("Invalid U matrix: input U matrix "
                                           "is not unitary, U^T * U != I.");
        }
        U_ = U;
    }

    /**
     * @brief Set up the print level.
     */
    void set_print(enum PrintLevel level) { print_level_ = level; }

    /**
     * @brief Compute the LO coefficient matrix under AO.
     *
     * @details The current U matrix stored in LocalizerBase::U_ will
     * be used as the initial guess. After calling this function,
     * LocalizerBase::U_ matrix is updated.
     *
     * LO coefficient matrix dimension: [nlo, nbasis]
     *
     * Relation between LO and AO via the LO coefficient matrix \f$C\f$ is
     * \f[
     * \psi_i = \sum_\mu C_{i\mu} \phi_{\mu},
     * \f]
     *
     * @return SharedMatrix: the LO coefficient matrix under AO.
     */
    virtual SharedMatrix compute() = 0;
};

/**
 * @brief Losc localization version 2.
 * @details This is the version used in the Losc2 paper (xxx).
 */
class LoscLocalizerV2 : public LocalizerBase {
  private:
    /**
     * @brief Hamiltonian matrix under AO used in Losc2 localization.
     * @details In Losc2 method, it is the DFA Hamiltonian.
     *
     * Dimension: [nbasis, nbasis].
     */
    SharedMatrix H_ao_;

    /**
     * @brief Dipole matrix under AO in order of x, y and z directions.
     * @details Dimension: [nbasis, nbasis].
     */
    vector<SharedMatrix> Dipole_ao_;

    /**
     * @brief Maximum iteration number of Jacobi-sweep algorithm for
     * localization. Default: 1000
     */
    size_t js_max_iter_ = 1000;

    /**
     * @brief If use random permutation in Jacobi-sweep algorithm or not.
     * Default: true
     */
    bool js_random_permutation_ = true;

    /**
     * @brief Convergence tolerance tolerance in Jacobi-sweep algorithm.
     * Default: 1e-10
     */
    double js_tol_ = 1e-10;

    /**
     * @brief Parameter \f$C\f$ in Losc2 localization cost function.
     * @details See Eq. (xxx) in Losc2 paper. Default: 1000
     */
    double para_c_ = 1000;

    /**
     * @brief Parameter \f$\gamma\f$ in Losc2 localization cost function.
     * @details See Eq. (xxx) in Losc2 paper. Default: 0.707
     */
    double para_gamma_ = 0.707;

    void message(std::string t, ...) const;

  public:
    /**
     * @param [in] C_lo_basis: LO basis coefficient matrix under AO. See
     * LoscLocalizerV2::C_lo_basis_.
     * @param [in] H_ao: Hamiltonian matrix under AO used in localization with
     * dimension [nbasis, nbasis].
     * @param [in] Dipole_ao: dipole matrix under AO in x, y and z order with
     * dimension [nbasis, nbasis].
     */
    LoscLocalizerV2(const SharedMatrix &C_lo_basis, const SharedMatrix &H_ao,
                    const vector<SharedMatrix> &Dipole_ao);

    /**
     * @param [in] tol: Jacobi-sweep algorithm convergence tolerance.
     */
    void set_tolerance(double tol) { js_tol_ = tol; }

    /**
     * @param [in] max_iter: Jacobi-sweep algorithm maximum iteration number.
     */
    void set_max_iteration(size_t max_iter) { js_max_iter_ = max_iter; }

    /**
     * @param [in] true_or_false: flag to do random permutation or not
     * Jacobi-sweep algorithm.
     */
    void set_random_permutation(bool true_or_false)
    {
        js_random_permutation_ = true_or_false;
    }

    /**
     * @brief Compute the LO coefficient matrix under AO for Losc localization
     * version 2.
     *
     * @return SharedMatrix: the LO coefficient matrix under AO.
     * @see LocalizerBase::compute()
     */
    SharedMatrix compute() override;
};

} // namespace losc
#endif //_LOSC_SRC_LOCALIZATION_H_
