/**
 * @file
 * @brief declaration relates to Losc curvature.
 */
#ifndef _LOSC_SRC_CURVATURE_H_
#define _LOSC_SRC_CURVATURE_H_

#include "matrix.h"
#include <memory>
#include <vector>

namespace losc {

using losc::Matrix;
using std::shared_ptr;
using std::vector;

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
     * @details Dimension: [nbasis, nlo]
     *
     * @par Availability
     * This matrix can be obtained from the losc::LocalizerBase::compute() or
     * any functions overwrite it from the derived localizer class.
     * @see losc::LoscLocalizerV2::compute()
     */
    shared_ptr<Matrix> C_lo_;

    /**
     * @brief The three-body integral <p|ii> matrix used in density fitting.
     * @details p is for fitbasis index and i is for LO index.
     *
     * <p|ii> matrix has dimension of [nfitbasis, nlo].
     * The integral is defined as
     * \f[
     * \langle \phi_p | \phi_i \phi_i \rangle
     * = \int \frac{\phi_p(\rm{r}) \phi_i(\rm{r'}) \phi_i(\rm{r'})}{|\rm{r}
     *   - \rm{r'}|} \mbox{d} \rm{r} \mbox{d} \rm{r'},
     * \f]
     *
     * @see losc::utils::convert_df_pmn2pii_blockwise().
     */
    shared_ptr<Matrix> df_pii_;

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
    shared_ptr<Matrix> df_Vpq_inverse_;

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
    shared_ptr<Matrix> grid_basis_value_;

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
    shared_ptr<vector<double>> grid_weight_;

  public:
    /**
     * @brief Curvature base class constructor.
     * @param [in] dfa: Type of DFA.
     * @param [in] C_lo: LO coefficient matrix under AO. See
     * CurvatureBase::C_lo_.
     * @param [in] df_pii: Three-body integral <p|ii> used in density fitting.
     * See CurvatureBase::df_pii_.
     * @param [in] df_Vpq_inverse: Inverse of Vpq matrix used in density
     * fitting. See CurvatureBase::df_Vpq_inverse_.
     * @param [in] grid_basis_value: Value of AO basis on grid. See
     * CurvatureBase::grid_basis_value_.
     * @param [in] grid_weight: Coefficient vector for numerical integral on
     * grid. See CurvatureBase::grid_weight_.
     */
    CurvatureBase(DFAType dfa, const shared_ptr<Matrix> &C_lo,
                  const shared_ptr<Matrix> &df_pii,
                  const shared_ptr<Matrix> &df_Vpq_inverse,
                  const shared_ptr<Matrix> &grid_basis_value,
                  const shared_ptr<vector<double>> &grid_weight);

    /**
     * @brief Compute the Losc curvature matrix.
     * @return shared_ptr<Matrix>: the Losc curvature matrix.
     */
    virtual shared_ptr<Matrix> compute() = 0;
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

    shared_ptr<Matrix> compute_kappa_J();
    shared_ptr<Matrix> compute_kappa_xc();

  public:
    /**
     * @brief Class constructor for curvature version 1.
     * @details See CurvatureV1::compute() for more details of curvature version
     * 1 matrix.
     * @param [in] dfa: Type of DFA.
     * @param [in] C_lo: LO coefficient matrix under AO. See
     * CurvatureBase::C_lo_.
     * @param [in] df_pii: Three-body integral <p|ii> used in density fitting.
     * See CurvatureBase::df_pii_.
     * @param [in] df_Vpq_inverse: Inverse of Vpq matrix used in density
     * fitting. See CurvatureBase::df_Vpq_inverse_.
     * @param [in] grid_basis_value: Value of AO basis on grid. See
     * CurvatureBase::grid_basis_value_.
     * @param [in] grid_weight: Coefficient vector for numerical integral on
     * grid. See CurvatureBase::grid_weight_.
     */
    CurvatureV1(DFAType dfa, const shared_ptr<Matrix> &C_lo,
                const shared_ptr<Matrix> &df_pii,
                const shared_ptr<Matrix> &df_Vpq_inverse,
                const shared_ptr<Matrix> &grid_basis_value,
                const shared_ptr<vector<double>> &grid_weight)
        : CurvatureBase(dfa, C_lo, df_pii, df_Vpq_inverse, grid_basis_value,
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
     * @return shared_ptr<Matrix>: The Losc curvature version 1 matrix.
     * @see The original Losc paper (https://doi.org/10.1093/nsr/nwx11)
     */
    virtual shared_ptr<Matrix> compute() override;
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
     * @param [in] df_pii: Three-body integral <p|ii> used in density fitting.
     * See CurvatureBase::df_pii_.
     * @param [in] df_Vpq_inverse: Inverse of Vpq matrix used in density
     * fitting. See CurvatureBase::df_Vpq_inverse_.
     * @param [in] grid_basis_value: Value of AO basis on grid. See
     * CurvatureBase::grid_basis_value_.
     * @param [in] grid_weight: Coefficient vector for numerical integral on
     * grid. See CurvatureBase::grid_weight_.
     */
    CurvatureV2(DFAType dfa, const shared_ptr<Matrix> &C_lo,
                const shared_ptr<Matrix> &df_pii,
                const shared_ptr<Matrix> &df_Vpq_inverse,
                const shared_ptr<Matrix> &grid_basis_value,
                const shared_ptr<vector<double>> &grid_weight)
        : CurvatureBase(dfa, C_lo, df_pii, df_Vpq_inverse, grid_basis_value,
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
     * @return shared_ptr<Matrix>: the Losc curvature version 2 matrix.
     * @see The Losc2 paper (xxx) for more details.
     */
    virtual shared_ptr<Matrix> compute() override;
};

namespace utils {

/**
 * @brief Convert density fitting matrix three-body integral <p|mn> into <p|ii>
 * matrix block by block.
 *
 * @details
 * i, j: LO index.\n
 * p, q: fitbasis index.\n
 * m, n: AO basis index.\n
 *
 * <p|mn> matrix has dimension of [nfitbasis, nbasis * (nbasis + 1) / 2].
 * <p|mn> integral is associated with on fitbasis and two AO basis, and it
 * is defined as
 * \f[
 * \langle \phi_p | \phi_m \phi_n \rangle
 * = \int \frac{\phi_p(\rm{r}) \phi_m(\rm{r'}) \phi_n(\rm{r'})}{|\rm{r}
 *   - \rm{r'}|} \mbox{d} \rm{r} \mbox{d} \rm{r'},
 * \f]
 *
 * <p|ii> matrix has dimension of [nfitbasis, nlo].
 * <p|ii> integral is associated with on fitbasis and one LO, and it
 * is defined as
 * \f[
 * \langle \phi_p | \phi_i \phi_i \rangle
 * = \int \frac{\phi_p(\rm{r}) \phi_i(\rm{r'}) \phi_i(\rm{r'})}{|\rm{r}
 *   - \rm{r'}|} \mbox{d} \rm{r} \mbox{d} \rm{r'},
 * \f]
 *
 * This function convert a block (block of fitbasis indices) <p|mn> matrix
 * into <p|ii> matrix. The format of <p|mn> block is required and illustrated as
 * the following code.
 * @code
 * for (size_t p = 0; p < block_size; ++p) {
 *     for (size_t m = 0; m < nbasis; ++m) {
 *         for (size_t n = 0; n <= m; ++n) {
 *             // set each element <p|mn> in the block as zero.
 *             df_pmn_block(p, m * (m + 1)/2 + n) = 0;
 *         }
 *     }
 * }
 * @endcode
 *
 * @note To build the whole <p|ii> matrix, you need to traverse all the <p|mn>
 * blocks.

 * @param [in] p_index: the corresponding fitbasis index in `df_pmn_block`.
 * For example, the i-th row in `df_pmn_block` matrix coresponding to
 * `p_index[i]`-th fitbasis.
 * @param [in] df_pmn_block: density fitting <p|mn> matrix block.
 * @param [in] C_lo: LO coefficient matrix under AO. See losc::LocalizerBase().
 * @param [in, out] df_pii: density fitting <p|ii> matrix. On exit, the rows
 * whose index store in `p_index` vector are updated.
 */
void convert_df_pmn2pii_blockwise(const vector<size_t> &p_index,
                                  const Matrix &df_pmn_block,
                                  const Matrix &C_lo, Matrix &df_pii);

} // namespace utils

} // namespace losc

#endif // _LOSC_SRC_CURVATURE_H_
