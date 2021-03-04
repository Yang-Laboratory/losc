/**
 * @file
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
 * @class __param__C_lo
 * @param [in] C_lo LOs coefficient matrix under AOs with dimension of
 * `[nbasis, nlo]`. The `i`-th column in `C_lo` matrix is the `i`-th LO.
 */

/**
 * @brief Base class for LOSC localization.
 */
class LocalizerBase {
  protected:
    /// @cond
    size_t nlo_;     /**< The number of LO. */
    size_t nbasis_;  /**< The number of AO basis. */
    size_t nsteps_;  /**< The number of iteration steps for the localization. */
    bool converged_; /**< The convergence of the localization. */

    // LO basis coefficient matrix under AOs.
    // Dimension: [nbasis, nlo]
    ConstRefMat C_lo_basis_;

    // Maximum iteration number of Jacobi-sweep algorithm for
    size_t js_max_iter_ = 1000;

    // If use random permutation in Jacobi-sweep algorithm or not.
    bool js_random_permutation_ = true;

    // Convergence tolerance tolerance in Jacobi-sweep algorithm.
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
    /// @endcond

  public:
    /**
     * @class __param__C_lo_basis
     * @param [in] C_lo_basis: LO basis coefficient matrix under AOs with
     * dimension of [nbasis, nlo]. The LO basis relates to the LOs via a unitary
     * transformation matrix. Usually, the converged COs from the associated
     * DFA are used as the LO basis. The `i`-th column of `C_lo_basis` is the
     * `i`-th LO basis.
     */

    /**
     * @brief Constructor of LocalizerBase
     * @copydoc __param__C_lo_basis
     */
    LocalizerBase(ConstRefMat &C_lo_basis);

    /**
     * @brief Deconstructor of LocalizerBase
     */
    virtual ~LocalizerBase();

    /**
     * Return the number of iteration steps for the most recent localization
     * performed by calling the `lo_U()` and `lo()` functions.
     */
    size_t steps() const { return nsteps_; }

    /**
     * Return the convergence of the most recent localization performed by
     * calling the `lo_U()` and `lo()` functions.
     */
    bool is_converged() const { return converged_; }

    /**
     * Return the localization cost function value for the given LOs.
     */
    virtual double cost_func(ConstRefMat &lo) const = 0;

    /**
     * @brief Set the maximum iteration number for localization.
     */
    void set_max_iter(size_t max_iter) { js_max_iter_ = max_iter; }

    /**
     * @brief Set the convergence of localization.
     */
    void set_convergence(double tol) { js_tol_ = tol; }

    /**
     * @brief Set to use random permutation for Jacobi-Sweep algorithm or not.
     * @param [in] flag `True` for doing random permutation and `False`
     * otherwise.
     */
    void set_random_permutation(bool flag) { js_random_permutation_ = flag; }

    /**
     * @class __param__guess
     * @param [in] guess Valid choices are
     * `["identity", "random", "random_fixed_seed"]`. Default to "identity"
     * - "identity": initial U matrix is set as an identity matrix.
     * - "random": initial U matrix is set as a random unitary matrix.
     * - "random_fixed_seed": initial U matrix is set as a random unitary matrix
     * with fixed random seed.
     */

    /**
     * @class __return__lo_U
     * @return a vector of matrix `rst` with size of 2. `rst[0]` is the LO
     * coefficient matrix. `rst[1]` is the corresponding U matrix.
     * The relation between the LOs and AOs is
     * \f$ \displaystyle
     * \psi_i = \sum_\mu C_{\mu i} \phi_{\mu},
     * \f$
     * in which \f$ \psi_i \f$ is the LO, \f$ C_{\mu i} \f$ is the LO
     * coefficient matrix and \f$ \phi_i \f$ is the AO.
     * The dimension of LO coefficient matrix is `[nbasis, nlo]`.
     * The relation between the LOs and LO basis is
     * \f$ \displaystyle
     * \psi_i = \sum_j U_{\mu i} \phi_\mu,
     * \f$
     * in which \f$ \psi_i \f$ is the LO, \f$ U_{\mu i} \f$ is the U matrix
     * and \f$ \phi_\mu \f$ is the LO basis.
     * The dimension of the U matrix is `[nlo, nlo]`.
     */

    /**
     * @brief Calculate the coefficient matrix of LOs on AOs and the unitary
     * transformation matrix.
     * @copydoc __param__guess
     * @copydoc __return__lo_U
     */
    virtual vector<LOSCMatrix> lo_U(const string &guess = "identity") = 0;

    /**
     * @class __param__U_guess__threshold
     * @param [in] U_guess the initial guess of U matrix. Its data will be
     * copied. Its unitarity will be verified and throw an exception
     * if the validation fails.
     * @param [in] threshold the threshold used to check the unitarity. Default
     * to 1e-8.
     */

    /**
     * @brief Calculate the LOs and the unitary transformation matrix.
     * @copydoc __param__U_guess__threshold
     * @copydoc __return__lo_U
     */
    virtual vector<LOSCMatrix> lo_U(ConstRefMat &U_guess,
                                    double threshold = 1e-8) = 0;

    /**
     * @class __param__L
     * @param [in, out] L A matrix with dimension of `[nbasis, nlo]`. The
     * content of L is ignored. At exit, it stores the LO coefficient matrix
     * under AO.
     */

    /**
     * @class __param__U
     * @param [in, out] U A matrix with dimension of `[nlo, nlo]`. At entrance,
     * it stores a unitary matrix as the initial guess for the localization.
     * The unitarity of the input U matrix is not verified. At exit, it is
     * updated by the localization process.
     */

    /**
     * @brief The C interface to calculate the LOs and the unitary
     * transformation matrix
     * @copydoc __param__L
     * @copydoc __param__U
     */
    virtual void C_API_lo_U(RefMat L, RefMat U) = 0;

    /**
     * @class __return__C_lo
     * @return The LO coefficient matrix `C_lo` under AO with dimension
     * `[nbasis, nlo]`. The `i`-th column of `C_lo` is the `i`-th LO.
     */

    /**
     * @brief Calculate the LOs' coefficient matrix under AO.
     * @copydoc __param__guess
     * @copydoc __return__C_lo
     * @see losc::LocalizerBase::lo_U()
     */
    virtual LOSCMatrix lo(const string &guess = "identity")
    {
        return lo_U(guess)[0];
    }

    /**
     * @brief Calculate the LOs' coefficient matrix under AO.
     * @copydoc __param__U_guess__threshold
     * @copydoc __return__C_lo
     * @see losc::LocalizerBase::lo_U()
     */
    virtual LOSCMatrix lo(ConstRefMat &U_guess, double threshold = 1e-8)
    {
        return lo_U(U_guess, threshold)[0];
    }
};

/**
 * @brief LOSC localization version 2.
 * @note This is the version used in the LOSC2 paper.
 * @see J. Phys. Chem. Lett. 2020, 11, 4, 1528–1535.
 */
class LocalizerV2 : public LocalizerBase {
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
                              const vector<LOSCMatrix> &D_lo,
                              const LOSCMatrix &H_lo, double &theta_val,
                              double &delta_val) const;
    /**
     * Internal function to rotate two orbitals: rotate D_lo, H_lo and U
     * matrices.
     */
    void js_rotate_one_pair(const size_t i, const size_t j, const double theta,
                            RefMat U, vector<LOSCMatrix> &D_lo,
                            LOSCMatrix &H_lo) const;

  public:
    /**
     * @brief Constructor of LocalizerV2
     * @copydoc __param__C_lo_basis
     * @param [in] H_ao Hamiltonian matrix under AO used in localization with
     * dimension `[nbasis, nbasis]`.
     * @param [in] Dipole_ao a vector of dipole matrices under AO (in the order
     * of x, y and z). The dimension of each dipole matrix is
     * `[nbasis, nbasis]`.
     */
    LocalizerV2(ConstRefMat &C_lo_basis, ConstRefMat &H_ao,
                const vector<RefConstMat> &Dipole_ao);

    /**
     * set the localization parameter gamma.
     */
    void set_gamma(double gamma) { gamma_ = gamma; }

    /**
     * set the localization parameter c.
     */
    void set_c(double c) { c_ = c; }

    /**
     * Return the relative cost function of LOSC localization v2.
     * @note The cost function is defined as:
     * \f[
     *  F = \sum_i (1 - \gamma)
     *      \Big (
     *          \langle i | \mathbf{r}^2| i \rangle - \langle i | \mathbf{r} |
     *          i \rangle ^2
     *      \Big) +
     *      \gamma C
     *      \Big (
     *      \langle i | H^2| i \rangle - \langle i | H | i \rangle ^2
     *      \Big).
     * \f]
     * Since the unitary transformation does not change the trace of a matrix,
     * The relative cost function is evaluated:
     * \f[
     *  F_{\rm{rel}} = \sum_i (1 - \gamma) \langle i | \mathbf{r} | i \rangle ^2
     *                 + \gamma C \langle i | H | i \rangle ^2.
     * \f]
     */
    double cost_func(ConstRefMat &lo) const override;

    /**
     * @brief Calculate the LOs and the unitary transformation matrix from the
     * LOSC localization v2.
     * @copydoc __param__guess
     * @copydoc __return__C_lo
     */
    virtual vector<LOSCMatrix> lo_U(const string &guess = "identity") override;

    /**
     * @brief Calculate the LOs and the unitary transformation matrix from the
     * LOSC localization v2.
     * @copydoc __param__U_guess__threshold
     * @copydoc __return__lo_U
     * @see losc::LocalizerV2::lo_U(const string&)
     */
    virtual vector<LOSCMatrix> lo_U(ConstRefMat &U_guess,
                                    double threshold = 1e-8) override;

    /**
     * @brief Calculate the LOs and the unitary transformation matrix in an
     * in-out way from the LOSC localization v2.
     * @copydoc __param__L
     * @copydoc __param__U
     */
    virtual void C_API_lo_U(RefMat L, RefMat U) override;
};

} // namespace losc
#endif //_LOSC_SRC_LOCALIZATION_H_
