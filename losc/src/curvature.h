/**
 * @file
 * @brief declaration relates to Losc curvature.
 */
#ifndef _LOSC_SRC_CURVATURE_H_
#define _LOSC_SRC_CURVATURE_H_

#include <matrix/matrix.h>
#include <memory>
#include <vector>

namespace losc {

using matrix::Matrix;
using std::vector;
using SharedMatrix = std::shared_ptr<Matrix>;
using SharedDoubleVector = std::shared_ptr<vector<double>>;

/**
 * @brief Type of density function approximation.
 */
enum DFAType {
    LDA,   /**< Local density approximation. */
    GGA,   /**< General gradient approximation. */
    B3LYP, /**< Hybrid type: B3LYP function. */
};

/**
 * @brief Base class for Losc curvature.
 */
class CurvatureBase {
  protected:
    DFAType dfa_type_; /**< type of DFA the curvature is related to. */

    size_t nlo_;       /**< number of LOs */
    size_t nbasis_;    /**< number of AOs basis. */
    size_t nfitbasis_; /**< number of fitbasis for density fitting in curvature
                          matrix construction. */
    size_t npts_;      /**< number of grid for numerical integration. */

    double para_alpha_ = 0.0;
    double para_beta_ = 0.0;

    /**
     * @brief LO coefficient matrix under AO.
     * @details Dimension: [nlo, nbasis]
     *
     * @par Availability
     * This matrix can be obtained from the losc::LocalizerBase::compute() or
     * any functions overwrite it from the derived localizer class.
     * @see losc::LoscLocalizerV2::compute()
     */
    SharedMatrix C_lo_;

    /**
     * @brief Three-body integral of <p|mn> matrix used in density fitting.
     * @details
     * Dimension: [nfitbasis, nbasis * (nbasis + 1) / 2]
     *
     * <p|mn> integral is defined as
     * \f[
     * \langle \phi_p | \phi_m \phi_n \rangle
     * = \int \frac{\phi_p(\rm{r}) \phi_m(\rm{r'}) \phi_n(\rm{r'})}{|\rm{r}
     *   - \rm{r'}|} \mbox{d} \rm{r} \mbox{d} \rm{r'},
     * \f]
     * where index p is for the fibbasis and m, n are for AO basis.
     *
     * @par Availability
     * You have to construct the matrix by yourself. To set each element,
     * you can do as the following,
     * @code
     * for (size_t p = 0; p < nfitbasis; ++p) {
     *     for (size_t m = 0; m < nbasis; ++m) {
     *         for (size_t n = 0; n <= m; ++n) {
     *             // set each element <p|mn> as zero.
     *             (*C_lo_)(p, m * (m + 1)/2 + n) = 0;
     *         }
     *     }
     * }
     * @endcode
     */
    SharedMatrix df_pmn_;

    /**
     * @brief Inverse of Vpq matrix used in density fitting.
     * @details
     * Dimension: [nfitbasis, nfitbasis]. Symmetric matrix.
     *
     * Element \f$Vpq\f$ is defined as
     * \f$\langle \phi_p| \frac{1}{\rm{r}_{12}} |\phi_q \rangle\f$
     * where p and q are the fitbasis index.
     *
     * @par Availability
     * You have to construct the matrix by yourself. To set each element,
     * you can do as the following,
     * @code
     * for (size_t p = 0; p < nfitbasis; ++p) {
     *     for (size_t q = 0; q < nfitbasis; ++q) {
     *         // set matrix element (p, q) as zero.
     *         (*df_Vpq_inverse_)(p, q) = 0.0;
     *     }
     * }
     * @endcode
     */
    SharedMatrix df_Vpq_inverse_;

    /**
     * @brief AO basis value on grid.
     * @details
     * Dimension: [npts, nbasis].
     *
     * @par Availability
     * You have to build the matrix by yourself. To construct the matrix,
     * you can do as the following,
     * @code
     * for (size_t ip = 0; ip < npts; ++ip) {
     *     for (size_t i = 0; i < nbasis; ++i) {
     *         // set the value of i-th AO
     *         // basis on grid point ip as zero.
     *         (*grid_basis_value_)(ip, i) = 0;
     *     }
     * }
     * @endcode
     */
    SharedMatrix grid_basis_value_;

    /**
     * @brief Coefficient for grid points used in numerical integral.
     * @details Dimension: npts
     *
     * @par Availability
     * You have to build the vector by yourself. To construct the vector,
     * you can do as the following,
     * @code
     * for (size_t ip = 0; ip < npts; ++ip) {
     *     (*grid_weight_)[ip] = 0.0; // set each grid coefficient as zero.
     * }
     * @endcode
     */
    SharedDoubleVector grid_weight_;

  public:
    /**
     * @brief Curvature base class constructor.
     * @param [in] dfa: Type of DFA.
     * @param [in] C_lo: LO coefficient matrix under AO. See
     * CurvatureBase::C_lo_.
     * @param [in] df_pmn: Three-body integral <p|mn> used in density fitting.
     * See CurvatureBase::df_pmn_.
     * @param [in] df_Vpq_inverse: Inverse of Vpq matrix used in density
     * fitting. See CurvatureBase::df_Vpq_inverse_.
     * @param [in] grid_basis_value: Value of AO basis on grid. See
     * CurvatureBase::grid_basis_value_.
     * @param [in] grid_weight: Coefficient vector for numerical integral on
     * grid. See CurvatureBase::grid_weight_.
     */
    CurvatureBase(enum DFAType dfa, const SharedMatrix &C_lo,
                  const SharedMatrix &df_pmn,
                  const SharedMatrix &df_Vpq_inverse,
                  const SharedMatrix &grid_basis_value,
                  const SharedDoubleVector &grid_weight);

    /**
     * @brief Compute the Losc curvature matrix.
     * @return SharedMatrix: the Losc curvature matrix.
     */
    virtual SharedMatrix compute() = 0;
};

/**
 * @brief Losc curvature class for version 1.
 * @details This class take the responsibility to generate curvature
 * version 1 matrix. Curvature version 1 is used in the original
 * Losc paper (https://doi.org/10.1093/nsr/nwx11).
 */
class CurvatureV1 : public CurvatureBase {
  private:
    /**
     * @brief Paramerter \f$C_x\f$ in curvature.
     * @details See Eq. (10) in the original Losc paper
     * (https://doi.org/10.1093/nsr/nwx111).
     */
    double para_cx_ = 0.930526;

    /**
     * @brief Paramerter \f$\tau\f$ in curvature.
     * @details See Eq. (10) in the original Losc paper
     * (https://doi.org/10.1093/nsr/nwx111).
     */
    double para_tau_ = 1.2378;

    SharedMatrix compute_kappa_J();
    SharedMatrix compute_kappa_xc();

  public:
    /**
     * @brief Class constructor for curvature version 1.
     * @details See CurvatureV1::compute() for more details of curvature version
     * 1 matrix.
     * @param [in] dfa: Type of DFA.
     * @param [in] C_lo: LO coefficient matrix under AO. See
     * CurvatureBase::C_lo_.
     * @param [in] df_pmn: Three-body integral <p|mn> used in density fitting.
     * See CurvatureBase::df_pmn_.
     * @param [in] df_Vpq_inverse: Inverse of Vpq matrix used in density
     * fitting. See CurvatureBase::df_Vpq_inverse_.
     * @param [in] grid_basis_value: Value of AO basis on grid. See
     * CurvatureBase::grid_basis_value_.
     * @param [in] grid_weight: Coefficient vector for numerical integral on
     * grid. See CurvatureBase::grid_weight_.
     */
    CurvatureV1(enum DFAType dfa, const SharedMatrix &C_lo,
                const SharedMatrix &df_pmn, const SharedMatrix &df_Vpq_inverse,
                const SharedMatrix &grid_basis_value,
                const SharedDoubleVector &grid_weight)
        : CurvatureBase(dfa, C_lo, df_pmn, df_Vpq_inverse, grid_basis_value,
                        grid_weight)
    {
    }

    /**
     * @brief Compute the Losc curvature version 1 matrix.
     * @details Dimension of curvature matrix is [nlo, nlo].
     *
     * Curvature version 1 matrix is defined as
     * \f[
     * \kappa_{ij} = (1 - \alpha -\beta) \int \int \frac{\rho_i(\rm{r})
     * \rho_j(\rm{r'})}{|\rm{r} - \rm{r'}|} \mbox{d} \rm{r} \mbox{d} \rm{r'}
     * - \frac{2\tau C_x}{3} (1 - \alpha) \int [\rho_i(\rm{r})]^{\frac{2}{3}}
     * [\rho_j(\rm{r})]^{\frac{2}{3}} \mbox{d} \rm{r},
     * \f]
     * where \f$\rho_i\f$ is the LO density \f$|\phi_i|^2\f$.
     *
     * @return SharedMatrix: The Losc curvature version 1 matrix.
     * @see The original Losc paper (https://doi.org/10.1093/nsr/nwx11)
     */
    virtual SharedMatrix compute() override;
};

/**
 * @brief Losc curvature class for version 2.
 * @details This class take the responsibility to generate curvature version 2
 * matrix. Curvature version 2 is used in the Losc2 paper (xxx). Check it out
 * for more details.
 */
class CurvatureV2 : public CurvatureBase {
  private:
    /**
     * @brief Paramerter \f$\zeta\f$ in curvature.
     * @details See Eq. (xxx) in the Losc2 paper.
     */
    double para_zeta_ = 8.0;

    /**
     * @brief Paramerter \f$C_x\f$ in curvature.
     * @details See Eq. (xxx) in the Losc2 paper.
     */
    double para_cx_ = 0.930526;

    /**
     * @brief Paramerter \f$\tau\f$ in curvature.
     * @details See Eq. (xxx) in the Losc2 paper.
     */
    double para_tau_ = 1.2378;

  public:
    /**
     * @brief Class constructor for curvature version 2.
     * @param [in] dfa: Type of DFA.
     * @param [in] C_lo: LO coefficient matrix under AO. See
     * CurvatureBase::C_lo_.
     * @param [in] df_pmn: Three-body integral <p|mn> used in density fitting.
     * See CurvatureBase::df_pmn_.
     * @param [in] df_Vpq_inverse: Inverse of Vpq matrix used in density
     * fitting. See CurvatureBase::df_Vpq_inverse_.
     * @param [in] grid_basis_value: Value of AO basis on grid. See
     * CurvatureBase::grid_basis_value_.
     * @param [in] grid_weight: Coefficient vector for numerical integral on
     * grid. See CurvatureBase::grid_weight_.
     */
    CurvatureV2(enum DFAType dfa, const SharedMatrix &C_lo,
                const SharedMatrix &df_pmn, const SharedMatrix &df_Vpq_inverse,
                const SharedMatrix &grid_basis_value,
                const SharedDoubleVector &grid_weight)
        : CurvatureBase(dfa, C_lo, df_pmn, df_Vpq_inverse, grid_basis_value,
                        grid_weight)
    {
    }

    /**
     * @brief Compute the Losc curvature version 2 matrix.
     *
     * @details Dimension of curvature matrix is [nlo, nlo].
     *
     * Curvature version 2 matrix is defined as
     * \f[
     * \kappa^2_{ij} = \mbox{erf}(\zeta S_{ij}) \sqrt{\kappa^1_{i, i}
     * \kappa^1_{j, j}} + \mbox{erfc}(\zeta S_{ij}) \kappa^1_{i, j},
     * \f]
     * where \f$\kappa^1\f$ is the Losc curvature version 1 matrix,
     * \f$\kappa^2\f$ is the Losc curvature version 2 matrix,
     * \f$ S_{ij} = \int \sqrt{|\rho_i(\rm{r}) \rho_j(\rm{r})|}
     * \mbox{d}\rm{r}\f$, and \f$\rho_i\f$ is the LO density \f$|\phi_i|^2\f$.
     *
     * @return SharedMatrix: the Losc curvature version 2 matrix.
     * @see The Losc2 paper (xxx) for more details.
     */
    virtual SharedMatrix compute() override;
};

} // namespace losc

#endif // _LOSC_SRC_CURVATURE_H_
