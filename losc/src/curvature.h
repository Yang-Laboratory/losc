/**
 * @file curvature.h
 * @brief C++ interface for the curvature matrix of LOSC.
 *
 * @note
 * 1. We use the famous Eigen template library to share data between the
 * construction of curvature matrix and users. The `Eigen::MatrixXd` is used to
 * represent matrices. To avoid the unnecessary copy of data from the user, we
 * use `Eigen::Ref<MatrixXd>` in the interface to map the needed matrices to the
 * existing matrices that are provided by the user. Thus, the life of all the
 * input matrices are controlled by the users, not us! Make sure these input
 * matrices are alive during computation of curvature matrix. We do not perform
 * any check to these conditions.
 * 2. Caution: a common pitfall related to the life of input matrices is passing
 * a temporary matrix to the curvature constructor. This could happen if you
 * pass a matrix expression, supported by Eigen, which is finally evalated into
 * a temporary matrix object (note, matrix and matrix expression are different
 * concept in Eigen). So keep in mind that always prepared the matrices at
 * first, then call the curvature constructor.
 * 3. All the input matrices from users are stored in column-wise, which follows
 * the default behavior of the Eigen library.
 * 4. We never change the input data from the user.
 */

#ifndef _LOSC_SRC_CURVATURE_H_
#define _LOSC_SRC_CURVATURE_H_

#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

namespace losc {

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::string;
using std::vector;
using ConstRefMat = const Eigen::Ref<const MatrixXd>;
using ConstRefVec = const Eigen::Ref<const VectorXd>;
using RefMat = Eigen::Ref<MatrixXd>;
using RefVec = Eigen::Ref<VectorXd>;

/**
 * @brief Density functional approximation (DFA) information.
 *
 * It includes the exchange types and their weights in the DFA.
 */
struct DFAInfo {
    const string &name; /**< DFA name. */
    double gga_x;       /**< GGA type exchange weight. */
    double hf_x;        /**< Hatree-Fock exchange weight. */

    DFAInfo(double gga_x_wt, double hf_x_wt, const string &name = "")
        : name{name}, gga_x{gga_x_wt}, hf_x{hf_x_wt}
    {
    }
};

/**
 * @brief Base class for LOSC curvature.
 */
class CurvatureBase {
  protected:
    const DFAInfo dfa_info_; /**< DFA information. */
    size_t nlo_;             /**< number of LOs */
    size_t nbasis_;          /**< number of AOs basis. */
    size_t nfitbasis_; /**< number of fitbasis for density fitting in curvature
                          matrix construction. */
    size_t npts_;      /**< number of grid for numerical integration. */

    /**
     * @brief LO coefficient matrix under AO.
     * @details Dimension: [nbasis, nlo]
     * \f[
     * \psi_p = C_{\mu p} \phi_{\mu}
     * \f]
     * in which $\psi_p$ is the \f$p\f$-th LO, \f$C_{\mu p}\f$ is the LO
     * coefficient matrix and \f$\phi_{\mu}\f$ is the \f$\mu\f$-th AO.
     *
     * This matrix can be obtained from the losc::LocalizerBase::compute() or
     * any functions overwrite it from the derived localizer class.
     *
     * @see losc::LoscLocalizerV2::compute()
     */
    ConstRefMat &C_lo_;

    /**
     * @brief The three-body integral \f$\langle p|ii \rangle\f$ matrix used in
     * density fitting.
     * @details Dimension: [nfitbasis, nlo].
     *
     * \f[
     * \langle \phi_p | \phi_i \phi_i \rangle
     * = \int \frac{\phi_p(\rm{r}) \phi_i(\rm{r'}) \phi_i(\rm{r'})}{|\rm{r}
     *   - \rm{r'}|} \mbox{d} \rm{r} \mbox{d} \rm{r'},
     * \f]
     * in which index \f$p\f$ is for fitbasis and index \f$i\f$ is for LOs.
     *
     * @par Availability
     * You have to construct the matrix by yourself. See
     * losc::utils::convert_df_pmn2pii_blockwise() for help.
     */
    ConstRefMat &df_pii_;

    /**
     * @brief Inverse of \f$V_{pq}\f$ matrix used in density fitting.
     * @details Dimension: [nfitbasis, nfitbasis].
     *
     * \f[
     * V_{pq} =
     * \langle \phi_p| \frac{1}{\rm{r}_{12}} |\phi_q \rangle,
     * \f] in which indices \f$p, q\f$ are for fitting basis.
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
    ConstRefMat &df_Vpq_inverse_;

    /**
     * @brief AO basis value on grid.
     * @details Dimension: [npts, nbasis]. \f$V_{\mu \nu}\f$ is the \f$\mu\f$-th
     * grid point value for the \f$\nu\f$-th AO.
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
    ConstRefMat &grid_basis_value_;

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
    ConstRefVec &grid_weight_;

  public:
    /**
     * @brief Curvature base class constructor.
     *
     * @param [in] dfa: Type of DFA.
     * @param [in] C_lo: LO coefficient matrix under AO with dimension
     * [nbasis, nlo]. See CurvatureBase::C_lo_.
     * @param [in] df_pii: Three-body integral <p|ii> used in density fitting
     * with dimension [nfitbasis, nlo]. See CurvatureBase::df_pii_.
     * @param [in] df_Vpq_inverse: Inverse of Vpq matrix used in density
     * fitting with dimension [nfitbasis, nfitbasis]. See
     * CurvatureBase::df_Vpq_inverse_.
     * @param [in] grid_basis_value: Value of AO basis on grid with dimension
     * [npts, nbasis]. See CurvatureBase::grid_basis_value_.
     * @param [in] grid_weight: Coefficient vector for numerical integral on
     * grid with size [npts]. See CurvatureBase::grid_weight_.
     */
    CurvatureBase(const DFAInfo &dfa_info, ConstRefMat &C_lo,
                  ConstRefMat &df_pii, ConstRefMat &df_Vpq_inverse,
                  ConstRefMat &grid_basis_value, ConstRefVec &grid_weight);

    /**
     * @brief Compute the LOSC curvature matrix.
     * @return MatrixXd: the LOSC curvature matrix with dimension
     * [nlo, nlo].
     */
    virtual MatrixXd kappa() = 0;
};

/**
 * @brief LOSC curvature class for version 1.
 * @details This class take the responsibility to generate curvature
 * version 1 matrix. Curvature version 1 is used in the original
 * LOSC paper (https://doi.org/10.1093/nsr/nwx11).
 */
class CurvatureV1 : public CurvatureBase {
  private:
    /**
     * @brief Paramerter \f$C_x\f$ in curvature.
     * @details See Eq. (10) in the original LOSC paper
     * (https://doi.org/10.1093/nsr/nwx111).
     */
    double cx_ = 0.930526;

    /**
     * @brief Paramerter \f$\tau\f$ in curvature.
     * @details See Eq. (10) in the original LOSC paper
     * (https://doi.org/10.1093/nsr/nwx111).
     */
    double tau_ = 1.2378;

    MatrixXd compute_kappa_J();
    MatrixXd compute_kappa_xc();

  public:
    /**
     * @brief Class constructor for curvature version 1.
     * @details See CurvatureV1::compute() for more details of curvature version
     * 1 matrix.
     * @param [in] dfa: Type of DFA.
     * @param [in] C_lo: LO coefficient matrix under AO with dimension
     * [nbasis, nlo]. See CurvatureBase::C_lo_.
     * @param [in] df_pii: Three-body integral <p|ii> used in density fitting
     * with dimension [nfitbasis, nlo]. See CurvatureBase::df_pii_.
     * @param [in] df_Vpq_inverse: Inverse of Vpq matrix used in density
     * fitting with dimension [nfitbasis, nfitbasis]. See
     * CurvatureBase::df_Vpq_inverse_.
     * @param [in] grid_basis_value: Value of AO basis on grid with dimension
     * [npts, nbasis]. See CurvatureBase::grid_basis_value_.
     * @param [in] grid_weight: Coefficient vector for numerical integral on
     * grid with size [npts]. See CurvatureBase::grid_weight_.
     */
    CurvatureV1(const DFAInfo &dfa_info, ConstRefMat &C_lo, ConstRefMat &df_pii,
                ConstRefMat &df_Vpq_inverse, ConstRefMat &grid_basis_value,
                ConstRefVec &grid_weight)
        : CurvatureBase(dfa_info, C_lo, df_pii, df_Vpq_inverse,
                        grid_basis_value, grid_weight)
    {
    }

    /**
     * @brief Compute the LOSC curvature version 1 matrix.
     * @return MatrixXd: The LOSC curvature version 1 matrix with dimension of
     * [nlo, nlo]. The curvature matrix is defined as \f$\kappa\f$ in Eq. 3
     * of the original LOSC paper (https://doi.org/10.1093/nsr/nwx11).
     */
    virtual MatrixXd kappa() override;
};

/**
 * @brief LOSC curvature class for version 2.
 * @details This class take the responsibility to generate curvature version 2
 * matrix. Curvature version 2 is used in the LOSC2 paper
 * (see J. Phys. Chem. Lett. 2020, 11, 4, 1528–1535).
 */
class CurvatureV2 : public CurvatureBase {
  private:
    /**
     * @brief Paramerter \f$\zeta\f$ in curvature.
     * @details See Eq. 8 in the LOSC2 paper.
     */
    double zeta_ = 8.0;

    /**
     * @brief Paramerter \f$C_x\f$ in curvature.
     * @details See Eq. 2 in the LOSC2 paper.
     */
    double cx_ = 0.930526;

    /**
     * @brief Paramerter \f$\tau\f$ in curvature.
     * @details See Eq. 2 in the LOSC2 paper.
     */
    double tau_ = 1.2378;

  public:
    /**
     * @brief Class constructor for curvature version 2.
     * @param [in] dfa: Type of DFA.
     * @param [in] C_lo: LO coefficient matrix under AO with dimension
     * [nbasis, nlo]. See CurvatureBase::C_lo_.
     * @param [in] df_pii: Three-body integral <p|ii> used in density fitting
     * with dimension [nfitbasis, nlo]. See CurvatureBase::df_pii_.
     * @param [in] df_Vpq_inverse: Inverse of Vpq matrix used in density
     * fitting with dimension [nfitbasis, nfitbasis]. See
     * CurvatureBase::df_Vpq_inverse_.
     * @param [in] grid_basis_value: Value of AO basis on grid with dimension
     * [npts, nbasis]. See CurvatureBase::grid_basis_value_.
     * @param [in] grid_weight: Coefficient vector for numerical integral on
     * grid with size [npts]. See CurvatureBase::grid_weight_.
     */
    CurvatureV2(const DFAInfo &dfa, ConstRefMat &C_lo, ConstRefMat &df_pii,
                ConstRefMat &df_Vpq_inverse, ConstRefMat &grid_basis_value,
                ConstRefVec &grid_weight)
        : CurvatureBase(dfa, C_lo, df_pii, df_Vpq_inverse, grid_basis_value,
                        grid_weight)
    {
    }

    /**
     * @brief Compute the LOSC curvature version 2 matrix.
     * @return MatrixXd: The LOSC curvature version 2 matrix with dimension of
     * [nlo, nlo]. The curvature matrix is defined as \f$\tilde{\kappa}\f$ in
     * Eq. 8 of the LOSC2 paper (J. Phys. Chem. Lett. 2020, 11, 4, 1528–1535).
     */
    virtual MatrixXd kappa() override;
};

/**
 * @brief LOSC library utils namespace.
 * @details Collection of helper functions and so on.
 */
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
 * The <p|mn> matrix has dimension of [nfitbasis, nbasis * (nbasis + 1) / 2].
 * The <p|mn> integral is associated with one fitbasis and two AO basis. It
 * is defined as
 * \f[
 * \langle \phi_p | \phi_m \phi_n \rangle
 * = \int \frac{\phi_p(\rm{r}) \phi_m(\rm{r'}) \phi_n(\rm{r'})}{|\rm{r}
 *   - \rm{r'}|} \mbox{d} \rm{r} \mbox{d} \rm{r'},
 * \f]
 *
 * The <p|ii> matrix has dimension of [nfitbasis, nlo].
 * The <p|ii> integral is associated with one fitbasis and one LO, and it
 * is defined as
 * \f[
 * \langle \phi_p | \phi_i \phi_i \rangle
 * = \int \frac{\phi_p(\rm{r}) \phi_i(\rm{r'}) \phi_i(\rm{r'})}{|\rm{r}
 *   - \rm{r'}|} \mbox{d} \rm{r} \mbox{d} \rm{r'},
 * \f]
 *
 * This function convert a block (block of fitbasis indices) of <p|mn> matrix
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
 * @param [in] p_index: the corresponding fitbasis index in `df_pmn_block`.
 * For example, the i-th row in `df_pmn_block` matrix coresponding to
 * `p_index[i]`-th fitbasis.
 * @param [in] df_pmn_block: density fitting <p|mn> matrix block.
 * @param [in] C_lo: LO coefficient matrix under AO. See losc::LocalizerBase().
 * @param [in, out] df_pii: density fitting <p|ii> matrix. On exit, the rows
 * whose index store in `p_index` vector are updated.
 *
 * @note To build the whole <p|ii> matrix, you need to traverse all the <p|mn>
 * blocks.
 */
void convert_df_pmn2pii_blockwise(const vector<size_t> &p_index,
                                  ConstRefMat &df_pmn_block, ConstRefMat &C_lo,
                                  RefMat df_pii);

} // namespace utils

} // namespace losc

#endif // _LOSC_SRC_CURVATURE_H_
