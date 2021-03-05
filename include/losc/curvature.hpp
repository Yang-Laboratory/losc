/**
 * @file
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
 * @class __param__Curvature
 * @param [in] Curvature LOSC curvature matrix with dimension of `[nlo, nlo]`.
 */

/**
 * @brief Density functional approximation (DFA) information.
 */
class DFAInfo {
    const string name_; /**< DFA name. */
    double gga_x_;      /**< Total weights of all GGA and LDA type exchanges. */
    double hf_x_;       /**< Total weights of Hatree-Fock exchange. */

  public:
    /**
     * @brief Constructor of DFAInfo
     * @param [in] gga_x: The total weights of all GGA and LDA type exchanges.
     * @param [in] hf_x: The total weights of HF exchanges.
     * @param [in] name: The name of the DFA. Default to an empty string.
     * @par Example
     * Taking B3LYP functional as an example. The B3LYP functional is
     * \f[
     *      E^{\rm{B3LYP}}_{\rm{xc}} = E^{\rm{LDA}}_{\rm{x}} +
     *      a_0 (E^{\rm{HF}}_{\rm{x}} - E^{\rm{LDA}}_{\rm{x}}) +
     *      a_x (E^{\rm{GGA}}_{\rm{x}} - E^{\rm{LDA}}_{\rm{x}}) +
     *      E^{\rm{LDA}}_{\rm{c}} + a_c (E^{\rm{GGA}}_{\rm{c}} -
     *      E^{\rm{LDA}}_{\rm{c}}),
     * \f]
     * in which exchanges end with suffix ``_x`` and correlations end with
     * suffix ``_c``, \f$ a_0=0.20 \f$, \f$ a_x=0.72 \f$ and \f$ a_c=0.81 \f$.
     * Only the exchanges are considered and the correlations are ignored.
     * The GGA and LDA exchanges are viewed as the same type. Therefore,
     * the total weights of GGA and LDA exchanges are
     * \f$ 1 - a_0 + a_x \times (1 - 1) = 1 - a_0 = 0.80 \f$,
     * and the total weights of HF exchanges is
     * \f$ a_0 = 0.20 \f$. To construct a ``DFAInfo`` object for B3LYP
     * functional, one should do the following
     * @code
     * b3lyp = DFAInfo(0.80, 0.20, "B3LYP")
     * @endcode
     */
    DFAInfo(double gga_x, double hf_x, const string &name = "")
        : name_{name}, gga_x_{gga_x}, hf_x_{hf_x}
    {
    }

    /**
     * @brief Return the name of the DFA.
     */
    const string &name() const { return name_; }

    /**
     * @brief Return the total weights of GGA-like exchanges of the DFA.
     */
    const double &gga_x() const { return gga_x_; }

    /**
     * @brief Return the total weights of HF exchanges of the DFA.
     */
    const double &hf_x() const { return hf_x_; }
};

/**
 * @brief Base class for LOSC curvature.
 */
class CurvatureBase {
  protected:
    /// \cond
    const DFAInfo dfa_info_; /**< DFA information. */
    size_t nlo_;             /**< number of LOs */
    size_t nfitbasis_; /**< number of fitbasis for density fitting in curvature
                          matrix construction. */
    size_t npts_;      /**< number of grid for numerical integration. */

    // three-body integral <p|ii> for density fitting.
    // Dimension: [nfitbasis, nlo]
    ConstRefMat df_pii_;

    // inverse of <p|1/r|q> matrix for density fitting.
    // Dimension: [nfitbasis, nfitbasis]
    ConstRefMat df_Vpq_inverse_;

    // LOs' value on grid points.
    // Dimension: [npts, nlo].
    ConstRefMat grid_lo_;

    // Coefficient for grid points used in numerical integral.
    // Dimension: npts
    ConstRefVec grid_weight_;

    /// \endcond

  public:
    /**
     * @class __param__CurvatureBase
     * @param [in] dfa Type of DFA.
     * @param [in] df_pii Three-body integral \f$ \langle p |ii \rangle \f$
     * used in density fitting. Index `p` is for fitbasis and index `i`
     * is for LOs. The dimension of `df_pii` is `[nfitbasis, nlo]`.
     * @param [in] df_Vpq_inv Inverse of \f$ \langle p | 1/\mathbf{r}
     * | q \rangle \f$ matrix used in density fitting. Index `p` and `q`
     * are for LOs. The dimension of `df_Vpq_inv` is `[nfitbasis, nfitbasis]`.
     * @param [in] grid_lo LOs' value on grid points with dimension of
     * `[npts, nlo]`.
     * @param [in] grid_weight Coefficient vector for numerical integral on
     * grid with size of `npts`.
     * @note This constructor requires all the matrices are allocated in
     * advance. This could be memory consuming. For example, the `grid_lo`
     * matrix can be very large.
     */

    /**
     * @brief Constructor of CurvatureBase.
     * @copydoc __param__CurvatureBase
     */
    CurvatureBase(const DFAInfo &dfa, ConstRefMat &df_pii,
                  ConstRefMat &df_Vpq_inv, ConstRefMat &grid_lo,
                  ConstRefVec &grid_weight);

    /**
     * @brief Deconstructor of CurvatureBase.
     */
    ~CurvatureBase() {}

    /* TODO:
     * Add a constructor that can support block-wise construction of grid_lo to
     * save memory.
     */

    /**
     * @brief Compute the LOSC curvature matrix.
     * @return the LOSC curvature matrix with dimension of `[nlo, nlo]`.
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

    /**
     * @class __param__K
     * @param [in, out] K a matrix with dimension of `[nlo, nlo]` at the input.
     * At exit, it is updated and stores the LOSC curvature matrix.
     */

    /**
     * @brief The C interface to compute the LOSC curvature matrix in place.
     * @copydoc __param__K
     */
    virtual void C_API_kappa(RefMat K) const = 0;
};

/**
 * @brief LOSC curvature class for version 1.
 * @details This class take the responsibility to generate curvature
 * version 1 matrix. Curvature version 1 is defined as \f$ \kappa \f$ in
 * Eq. 3 in the original LOSC paper (https://doi.org/10.1093/nsr/nwx11).
 * @see losc::CurvatureV2
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
     * @copydoc __param__CurvatureBase
     */
    CurvatureV1(const DFAInfo &dfa, ConstRefMat &df_pii,
                ConstRefMat &df_Vpq_inv, ConstRefMat &grid_lo,
                ConstRefVec &grid_weight)
        : CurvatureBase(dfa, df_pii, df_Vpq_inv, grid_lo, grid_weight)
    {
    }

    /**
     * @brief Deconstructor of CurvatureV1.
     */
    ~CurvatureV1() {}

    /**
     * @brief The C interface to calculate LOSC curvature version 1 in place.
     * @copydoc __param__K
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
 * matrix. The curvature matrix is defined as \f$ \tilde{\kappa} \f$ in
 * Eq. 8 of the LOSC2 paper (J. Phys. Chem. Lett. 2020, 11, 4, 1528-1535).
 */
class CurvatureV2 : public CurvatureBase {
  private:
    /**
     * @brief Paramerter \f$ \zeta \f$ in curvature.
     * @details See Eq. 8 in the LOSC2 paper.
     */
    double zeta_ = 8.0;

    /**
     * @brief Paramerter \f$ C_x \f$ in curvature.
     * @details See Eq. 2 in the LOSC2 paper.
     */
    double cx_ = 0.930526;

    /**
     * @brief Paramerter \f$ \tau \f$ in curvature.
     * @details See Eq. 2 in the LOSC2 paper.
     */
    double tau_ = 1.2378;

  public:
    /**
     * @brief Class constructor for curvature version 2.
     * @copydoc __param__CurvatureBase
     */
    CurvatureV2(const DFAInfo &dfa, ConstRefMat &df_pii,
                ConstRefMat &df_Vpq_inv, ConstRefMat &grid_lo,
                ConstRefVec &grid_weight)
        : CurvatureBase(dfa, df_pii, df_Vpq_inv, grid_lo, grid_weight)
    {
    }

    /**
     * @brief Deconstructor of CurvatureV2.
     */
    ~CurvatureV2() {}

    /**
     * @brief The C interface to calculate LOSC curvature version 2 in place.
     * @copydoc __param__K
     */
    virtual void C_API_kappa(RefMat K) const override;

    /**
     * @brief Set parameter \f$ \tau \f$.
     */
    void set_tau(double tau) { tau_ = tau; }

    /**
     * @brief Set parameter \f$ \zeta \f$.
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
 * @param [in] p_index the corresponding fitbasis index in `df_pmn_block`.
 * For example, the i-th row in `df_pmn_block` matrix coresponding to
 * `p_index[i]`-th fitbasis.
 * @param [in] df_pmn_block density fitting <p|mn> matrix block.
 * @copydoc __param__C_lo
 * @param [in, out] df_pii density fitting <p|ii> matrix. On exit, the rows
 * whose index store in `p_index` vector are updated.
 * @note To build the whole <p|ii> matrix, you need to traverse all the <p|mn>
 * blocks.
 */
void convert_df_pmn2pii_blockwise(const vector<size_t> &p_index,
                                  ConstRefMat &df_pmn_block, ConstRefMat &C_lo,
                                  RefMat df_pii);

} // namespace utils

} // namespace losc

#endif // _LOSC_SRC_CURVATURE_H_
