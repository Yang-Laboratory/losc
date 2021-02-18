/**
 * @file curvature.hpp
 * @brief C++ interface for the curvature matrix of LOSC.
 */

#ifndef _LOSC_INCLUDE_LOSC_CURVATURE_HPP_
#define _LOSC_INCLUDE_LOSC_CURVATURE_HPP_

#include <losc/eigen_def.hpp>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace losc {

using std::string;

/**
 * @brief Density functional approximation (DFA) information.
 *
 * It includes the exchange types and their weights in the DFA.
 */
class DFAInfo {
    const string name_; /**< DFA name. */
    double gga_x_;      /**< Total weights of all GGA and LDA type exchanges. */
    double hf_x_;       /**< Total weights of Hatree-Fock exchange. */

  public:
    DFAInfo(double gga_x, double hf_x, const string &name = "")
        : name_{name}, gga_x_{gga_x}, hf_x_{hf_x}
    {
    }
    const string &name() const { return name_; }
    const double &gga_x() const { return gga_x_; }
    const double &hf_x() const { return hf_x_; }
};

/**
 * @brief Base class for LOSC curvature.
 */
class CurvatureBase {
  protected:
    const DFAInfo dfa_info_; /**< DFA information. */
    size_t nlo_;             /**< number of LOs */
    size_t nfitbasis_; /**< number of fitbasis for density fitting in curvature
                          matrix construction. */
    size_t npts_;      /**< number of grid for numerical integration. */

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
     */
    ConstRefMat df_pii_;

    /**
     * @brief Inverse of \f$V_{pq}\f$ matrix used in density fitting.
     * @details Dimension: [nfitbasis, nfitbasis].
     *
     * \f[
     * V_{pq} =
     * \langle \phi_p| \frac{1}{\rm{r}_{12}} |\phi_q \rangle,
     * \f] in which indices \f$p, q\f$ are for fitting basis.
     */
    ConstRefMat df_Vpq_inverse_;

    /**
     * @brief LOs' value on grid points.
     * @details Dimension: [npts, nlo].
     */
    ConstRefMat grid_lo_;

    /**
     * @brief Coefficient for grid points used in numerical integral.
     * @details Dimension: npts
     */
    ConstRefVec grid_weight_;

  public:
    /**
     * @brief Curvature base class constructor.
     *
     * @param [in] dfa: Type of DFA.
     * @param [in] df_pii: Three-body integral <p|ii> used in density fitting
     * with dimension of [nfitbasis, nlo]. See CurvatureBase::df_pii_.
     * @param [in] df_Vpq_inverse: Inverse of Vpq matrix used in density
     * fitting with dimension of [nfitbasis, nfitbasis]. See
     * CurvatureBase::df_Vpq_inverse_.
     * @param [in] grid_lo: LOs' value on grid points with dimension of
     * [npts, nlo]. See CurvatureBase::grid_lo_.
     * @param [in] grid_weight: Coefficient vector for numerical integral on
     * grid with size of [npts]. See CurvatureBase::grid_weight_.
     *
     * @note
     * This constructor requires all the pre-allocation of the whole grid matrix
     * of LOs (grid_lo). This grid_lo matrix could be very large.
     */
    CurvatureBase(const DFAInfo &dfa_info, ConstRefMat &df_pii,
                  ConstRefMat &df_Vpq_inverse, ConstRefMat &grid_lo,
                  ConstRefVec &grid_weight);

    /**
     * TODO:
     * Add a constructor that can support block-wise construction of grid_lo to
     * save memory.
     */

    /**
     * @brief The C interface to compute the LOSC curvature matrix in place.
     *
     * @param K [in, out]: a matrix with dimension of [nlo, nlo] at the input.
     * At exit, it is updated and stores the LOSC curvature matrix.
     */
    virtual void C_API_kappa(RefMat K) const = 0;

    /**
     * @brief Compute the LOSC curvature matrix.
     * @return LOSCMatrix: the LOSC curvature matrix with dimension of
     * [nlo, nlo].
     */
    virtual LOSCMatrix kappa() const
    {
        LOSCMatrix K(nlo_, nlo_);
        C_API_kappa(K);
        return std::move(K);
    }

    /**
     * Return number of LOs.
     */
    size_t nlo() const { return nlo_; }

    /**
     * Return number of fit basis.
     */
    size_t nfitbasis() const { return nfitbasis_; }

    /**
     * Return number of grid points.
     */
    size_t npts() const { return npts_; }
};

/**
 * @brief LOSC curvature class for version 1.
 * @details This class take the responsibility to generate curvature
 * version 1 matrix. Curvature version 1 is defined as \f$ \kappa \f$ in
 * Eq. 3 in the original LOSC paper (https://doi.org/10.1093/nsr/nwx11).
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

    LOSCMatrix compute_kappa_J() const;
    LOSCMatrix compute_kappa_xc() const;

  public:
    /**
     * @brief Class constructor for curvature version 1.
     * @details See CurvatureV1::compute() for more details of curvature version
     * 1 matrix.
     * @param [in] dfa: Type of DFA.
     * @param [in] df_pii: Three-body integral <p|ii> used in density fitting
     * with dimension of [nfitbasis, nlo]. See CurvatureBase::df_pii_.
     * @param [in] df_Vpq_inverse: Inverse of Vpq matrix used in density
     * fitting with dimension of [nfitbasis, nfitbasis]. See
     * CurvatureBase::df_Vpq_inverse_.
     * @param [in] grid_lo: LOs' value on grid points with dimension of
     * [npts, nlo]. See CurvatureBase::grid_basis_value_.
     * @param [in] grid_weight: Coefficient vector for numerical integral on
     * grid with size of [npts]. See CurvatureBase::grid_weight_.
     */
    CurvatureV1(const DFAInfo &dfa_info, ConstRefMat &df_pii,
                ConstRefMat &df_Vpq_inverse, ConstRefMat &grid_lo,
                ConstRefVec &grid_weight)
        : CurvatureBase(dfa_info, df_pii, df_Vpq_inverse, grid_lo, grid_weight)
    {
    }

    using CurvatureBase::kappa;

    /**
     * @see
     * losc::CurvatureBase::C_API_kappa(): calculate curvature matrix.
     */
    virtual void C_API_kappa(RefMat K) const override;

    /**
     * @brief Set parameter tau.
     */
    void set_tau(double tau) { tau_ = tau; }
};

/**
 * @brief LOSC curvature class for version 2.
 * @details This class take the responsibility to generate curvature version 2
 * matrix. The curvature matrix is defined as \f$\tilde{\kappa}\f$ in
 * Eq. 8 of the LOSC2 paper (J. Phys. Chem. Lett. 2020, 11, 4, 1528-1535).
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
     * @param [in] df_pii: Three-body integral <p|ii> used in density fitting
     * with dimension of [nfitbasis, nlo]. See CurvatureBase::df_pii_.
     * @param [in] df_Vpq_inverse: Inverse of Vpq matrix used in density
     * fitting with dimension of [nfitbasis, nfitbasis]. See
     * CurvatureBase::df_Vpq_inverse_.
     * @param [in] grid_lo: LOs' value on grid points with dimension of
     * [npts, nlo]. See CurvatureBase::grid_lo_.
     * @param [in] grid_weight: Coefficient vector for numerical integral on
     * grid with size of [npts]. See CurvatureBase::grid_weight_.
     */
    CurvatureV2(const DFAInfo &dfa, ConstRefMat &df_pii,
                ConstRefMat &df_Vpq_inverse, ConstRefMat &grid_lo,
                ConstRefVec &grid_weight)
        : CurvatureBase(dfa, df_pii, df_Vpq_inverse, grid_lo, grid_weight)
    {
    }

    using CurvatureBase::kappa;
    /**
     * @see
     * losc::CurvatureBase::C_API_kappa(): calculate curvature matrix.
     */
    virtual void C_API_kappa(RefMat K) const override;

    /**
     * @brief Set parameter tau.
     */
    void set_tau(double tau) { tau_ = tau; }

    /**
     * @brief Set parameter zeta.
     */
    void set_zeta(double zeta) { zeta_ = zeta; }
};

/**
 * @brief LOSC library utils namespace.
 * @details Collection of helper functions and so on.
 */
namespace utils {

using std::vector;

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
