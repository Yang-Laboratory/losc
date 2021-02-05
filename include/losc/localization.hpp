/**
 * @file localization.hpp
 * @brief C++ interface for the LOSC localization.
 */

#ifndef _LOSC_INCLUDE_LOSC_LOCALIZATION_HPP_
#define _LOSC_INCLUDE_LOSC_LOCALIZATION_HPP_

#include <cstddef>
#include <losc/eigen_def.hpp>
#include <string>
#include <vector>

namespace losc {

using std::string;
using std::vector;

/**
 * @brief Base class for LOSC localization.
 */
class LocalizerBase {
  protected:
    size_t nlo_;    /**<number of LO. */
    size_t nbasis_; /**< number of AO basis. */

    /**
     * @brief LO basis coefficient matrix under AOs.
     * @details Dimension: [nbasis, nlo]
     *
     * Relation between LO basis and AO via the coefficient matrix \f$C\f$ is
     * \f[
     * \psi_i = \sum_\mu C_{\mu i} \phi_{\mu},
     * \f]
     * where \f$ \psi_i \f$ is the LO basis and \f$ \phi_\mu \f$ is the AO.
     *
     * @note
     * Usually, LOSC localization uses the converged DFA COs as the LO basis.
     */
    ConstRefMat C_lo_basis_;

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
     * Internal function to set the initial U matrix.
     * @param U [in, out]: a matrix with dimension of [nlo, nlo]. It stores a
     * unitary matrix based on `guess` at the exit, and serves as the initial
     * guess for the localization.
     */
    void set_u_guess(RefMat U, const string &guess) const;

    /**
     * Internal function to set the initial U matrix.
     * @param U [in, out]: the initial U matrix is copied from `U_guess`.
     */
    void set_u_guess(RefMat U, ConstRefMat &U_guess, double threshold) const;

  public:
    /**
     * @brief LocalizerBase class constructor
     * @param [in] C_lo_basis: LO basis coefficient matrix under AOs with
     * dimension of [nbasis, nlo]. See LocalizerBase::C_lo_basis_.
     *
     * @note
     * 1. The initial guess is to set an identity U matrix.
     */
    LocalizerBase(ConstRefMat &C_lo_basis);
    virtual ~LocalizerBase();

    /**
     * @brief Set the maximum iteration number for localization.
     */
    void set_max_iter(size_t max_iter) { js_max_iter_ = max_iter; }

    /**
     * @brief Set the convergence of localization.
     */
    void set_convergence(double tol) { js_tol_ = tol; }

    /**
     * @param [in] true_or_false: flag to do random permutation or not
     * Jacobi-sweep algorithm.
     */
    void set_random_permutation(bool flag) { js_random_permutation_ = flag; }

    /**
     * @brief Calculate the LOs and the unitary transformation matrix.
     *
     * The calculated LOs are expanded under AOs via the following relation
     * \f[
     * \psi_i = \sum_\mu C_{\mu i} \phi_{\mu},
     * \f]
     * in which \f$ \psi_i \f$ is the LO, \f$ C_{\mu i} \f$ is the LO
     * coefficient matrix and \f$ \phi_i \f$ is the AO.
     * Dimension of LO coefficient matrix: [nbasis, nlo].

     * The unitary matrix is associated with the LOs and LOs basis via the
     * following relation.
     * \f[
     * \psi_i = \sum_j U_{i\mu} \phi_\mu,
     * \f]
     * in which \f$ psi_i\f$ is the LO, \f$ U_{i\mu} \f$ is the unitary matrix
     * and \f$ \phi_\mu \f$ is the LO basis that expands the LO.
     * Dimension of the U matrix: [nlo, nlo].
     *
     * @param guess: Valid choices are
     * ["identity", "random", "random_fixed_seed"]. Default to "identity"
     * "identity": initial U matrix is set as an identity matrix.
     * "random": initial U matrix is set as a random unitary matrix.
     * "random_fixed_seed": initial U matrix is set as a random unitary matrix
     * with fixed random seed.
     *
     * @return a size 2 vector of matrix `rst`. `rst[0]` is the LO coefficient
     * matrix. `rst[1]` is the corresponding U matrix.
     */
    virtual vector<MatrixXd> lo_U(const string &guess = "identity") const = 0;

    /**
     * @brief Calculate the LOs and the unitary transformation matrix.
     *
     * @param U_guess: the initial guess of U matrix. Its data will be copied
     * for localization. Its unitarity will be verified and throw an exception
     * if the validation fails.
     * @param threshold: the threshold used to check the unitarity. Default to
     * 1e-8.
     */
    virtual vector<MatrixXd> lo_U(ConstRefMat &U_guess,
                                  double threshold = 1e-8) const = 0;

    /**
     * @brief The C interface to calculate the LOs and the unitary
     * @param L [in, out]: Dimension of [nbasis, nlo]. The content of L is
     * ignored. At exit, it is stores to be the LO coefficient matrix.
     * @param U [in, out]: Dimension of [nlo, nlo]. At entrance, it stores a
     * unitary matrix as the initial guess for the localization. The unitarity
     * of the input U matrix is not verified. At exit, it is updated by the
     * localization process.
     */
    virtual void C_API_lo_U(RefMat L, RefMat U) const = 0;

    /**
     * @brief Calculate the LOs' coefficient matrix under AO.
     *
     * @param guess: a string represents the initial guess. See `lo_U`.
     * @return the LOs' coefficient matrix under AO.
     */
    virtual MatrixXd lo(const string &guess = "identity") const
    {
        return lo_U(guess)[0];
    }

    /**
     * @brief Calculate the LOs' coefficient matrix under AO.
     *
     * @param U_guess: the initial guess of U matrix. See `lo_U`.
     * @return the LOs' coefficient matrix under AO.
     */
    virtual MatrixXd lo(ConstRefMat &U_guess, double threshold = 1e-8) const
    {
        return lo_U(U_guess, threshold)[0];
    }
};

/**
 * @brief LOSC localization version 2.
 * @note This is the version used in the LOSC2 paper.
 * @see J. Phys. Chem. Lett. 2020, 11, 4, 1528–1535.
 */
class LoscLocalizerV2 : public LocalizerBase {
  private:
    /**
     * @brief Hamiltonian matrix under AO used in LOSC localization v2.
     * @details In LOSC2 method, it is just the Hamiltonian of the parent DFA.
     * Dimension: [nbasis, nbasis].
     */
    ConstRefMat H_ao_;

    /**
     * @brief Dipole matrix under AO in order of x, y and z directions.
     * @details Dimension: [nbasis, nbasis].
     */
    const vector<RefConstMat> Dipole_ao_;

    /**
     * @brief Parameter \f$ C \f$ in the cost function of localization v2.
     * @see Eq. 7 in the LOSC2 paper. Default to 1000.
     */
    double c_ = 1000;

    /**
     * @brief Parameter \f$\gamma\f$ in the cost function of localization v2.
     * @see Eq. 7 in the LOSC2 paper. Default to 0.707.
     */
    double gamma_ = 0.707;

    /**
     * Internal function to calculate the optimal rotation angle.
     */
    void js_optimize_one_pair(const size_t i, const size_t j,
                              const vector<MatrixXd> &D_lo,
                              const MatrixXd &H_lo, double &theta_val,
                              double &delta_val) const;
    /**
     * Internal function to rotate two orbitals: rotate D_lo, H_lo and U
     * matrices.
     */
    void js_rotate_one_pair(const size_t i, const size_t j, const double theta,
                            RefMat U, vector<MatrixXd> &D_lo,
                            MatrixXd &H_lo) const;

  public:
    /**
     * @param [in] C_lo_basis: LO basis coefficient matrix under AO with
     * dimension [nbasis, nlo]. See LoscLocalizerV2::C_lo_basis_.
     * @param [in] H_ao: Hamiltonian matrix under AO used in localization with
     * dimension [nbasis, nbasis].
     * @param [in] Dipole_ao: dipole matrix under AO in x, y and z order with
     * dimension [nbasis, nbasis].
     */
    LoscLocalizerV2(ConstRefMat &C_lo_basis, ConstRefMat &H_ao,
                    const vector<RefConstMat> &Dipole_ao);

    /**
     * @brief Calculate the LOs and the unitary transformation matrix from the
     * LOSC localization v2.
     *
     * @param guess: a string represents the initial guess.
     * See LocalizerBase::lo_U.
     * @return the LOs' coefficient matrix under AO.
     * @see LocalizerBase::lo_U()
     */
    virtual vector<MatrixXd>
    lo_U(const string &guess = "identity") const override;

    /**
     * @brief Calculate the LOs and the unitary transformation matrix from the
     * LOSC localization v2.
     *
     * @param U_guess: the initial guess of U matrix. See LocalizerBase::lo_U.
     * @param threshold: the threshold used to check the unitarity. Default to
     * 1e-8.
     * @see LocalizerBase::lo_U()
     */
    virtual vector<MatrixXd> lo_U(ConstRefMat &U_guess,
                                  double threshold = 1e-8) const override;

    /**
     * @brief Calculate the LOs and the unitary transformation matrix in an
     * in-out way from the LOSC localization v2.
     *
     * @see
     * LocalizerBase::lo_U
     */
    virtual void C_API_lo_U(RefMat L, RefMat U) const override;
};

} // namespace losc
#endif //_LOSC_SRC_LOCALIZATION_H_
