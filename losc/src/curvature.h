#ifndef _LOSC_SRC_CURVATURE_H_
#define _LOSC_SRC_CURVATURE_H_

#include <memory>
#include <vector>
#include <matrix/matrix.h>

namespace losc {

using std::vector;
using matrix::Matrix;
using SharedMatrix = std::shared_ptr<Matrix>;
using SharedDoubleVector = std::shared_ptr<vector<double>>;

enum DFAType {
    LDA,
    GGA,
    B3LYP,
};

class CurvatureBase {
    protected:
    DFAType dfa_type_;

    size_t nlo_;        /* number of LOs */
    size_t nbasis_;     /* number of AOs */
    size_t nfitbasis_;  /* number of fitbasis */
    size_t npts_;       /* number of grid */

    double para_alpha_ = 0.0;
    double para_beta_ = 0.0;

    /**
     * LO coefficient matrix under AO.
     * Dimension: [nlo, nbasis].
     * It can be obtained from `losc::LocalizerBase::compute()` function.
     */
    SharedMatrix C_lo_;

    /**
     * Three-body integral of <p|mn> matrix used in density fitting.
     * Dimension: [nfitbasis, nbasis * (nbasis + 1) / 2].
     * Integral <p|mn> is defined as <p|1/r_{12}|mn>, where
     * m, n is the AO basis index and p is the fitbasis index.
     * The location of <p|mn> element is [p, m*(m+1)/2 + n].
     */
    SharedMatrix df_pmn_;

    /**
     * Inverse of Vpq matrix used in density fitting.
     * Dimension: [nfitbasis, nfitbasis].
     * Element Vpq is defined as <p|1/r_{12}|q> where p and q are the fitbasis index.
     */
    SharedMatrix df_Vpq_inverse_;

    /**
     * AO basis value on grid.
     * Dimension: [npts, nbasis].
     */
    SharedMatrix grid_basis_value_;

    /**
     * coefficient vector for numerical integral on grid.
     * Dimension: npts.
     */
    SharedDoubleVector grid_weight_;

    public:
    /**
     * @ param [in] dfa: type of DFA.
     * @ param [in] C_lo: LO coefficient matrix under AO.
     * @ param [in] df_pmn: three-body integral <p|mn> used in density fitting.
     * @ param [in] df_Vpq_inverse: inverse of Vpq matrix used in density fitting.
     * @ param [in] grid_basis_value: value of AO basis on grid.
     * @ param [in] grid_weight: coefficient vector for numerical integral on grid.
     */
    CurvatureBase(enum DFAType dfa, const SharedMatrix& C_lo, const SharedMatrix& df_pmn,
                  const SharedMatrix& df_Vpq_inverse, const SharedMatrix& grid_basis_value,
                  const SharedDoubleVector& grid_weight);

    /**
     * Compute the Losc curvature matrix.
     *
     * @ return: the Losc curvature matrix.
     */
    virtual SharedMatrix compute() = 0;
};

/**
 * Losc curvature version 1. This is the version used in the original
 * Losc paper (https://doi.org/10.1093/nsr/nwx11).
 */
class CurvatureV1 : public CurvatureBase
{
    private:
    /**
     * Paramerter $C_x$ in curvature.
     * See Eq. (10) in the original Losc paper (https://doi.org/10.1093/nsr/nwx111).
     */
    double para_cx_ = 0.930526;

    /**
     * Paramerter $\tau$ in curvature.
     * See Eq. (10) in the original Losc paper (https://doi.org/10.1093/nsr/nwx111).
     */
    double para_tau_ = 1.2378;

    SharedMatrix compute_kappa_J();
    SharedMatrix compute_kappa_xc();

    public:
    /**
     * @ param [in] dfa: type of DFA.
     * @ param [in] C_lo: LO coefficient matrix under AO.
     * @ param [in] df_pmn: three-body integral <p|mn> used in density fitting.
     * @ param [in] df_Vpq_inverse: inverse of Vpq matrix used in density fitting.
     * @ param [in] grid_basis_value: value of AO basis on grid.
     * @ param [in] grid_weight: coefficient vector for numerical integral on grid.
     */
    CurvatureV1(enum DFAType dfa, const SharedMatrix& C_lo, const SharedMatrix& df_pmn,
                const SharedMatrix& df_Vpq_inverse, const SharedMatrix& grid_basis_value,
                const SharedDoubleVector& grid_weight)
        : CurvatureBase(dfa, C_lo, df_pmn, df_Vpq_inverse, grid_basis_value, grid_weight) {}

    /**
     * Compute the Losc curvature version 1 matrix.
     *
     * @ return: the Losc curvature matrix.
     */
    virtual SharedMatrix compute() override;
};

/**
 * Losc curvature version 2. This is the version used in the Losc2 paper (xxx).
 */
class CurvatureV2 : public CurvatureBase
{
    private:
    /**
     * Paramerter $\zeta$ in curvature.
     * See Eq. (xxx) in the Losc2 paper.
     */
    double para_zeta_ = 8.0;

    /**
     * Paramerter $C_x$ in curvature.
     * See Eq. (xxx) in the Losc2 paper.
     */
    double para_cx_ = 0.930526;

    /**
     * Paramerter $\tau$ in curvature.
     * See Eq. (xxx) in the Losc2 paper.
     */
    double para_tau_ = 1.2378;

    public:

    /**
     * @ param [in] dfa: type of DFA.
     * @ param [in] C_lo: LO coefficient matrix under AO.
     * @ param [in] df_pmn: three-body integral <p|mn> used in density fitting.
     * @ param [in] df_Vpq_inverse: inverse of Vpq matrix used in density fitting.
     * @ param [in] grid_basis_value: value of AO basis on grid.
     * @ param [in] grid_weight: coefficient vector for numerical integral on grid.
     */
    CurvatureV2(enum DFAType dfa, const SharedMatrix& C_lo, const SharedMatrix& df_pmn,
                const SharedMatrix& df_Vpq_inverse, const SharedMatrix& grid_basis_value,
                const SharedDoubleVector& grid_weight)
        : CurvatureBase(dfa, C_lo, df_pmn, df_Vpq_inverse, grid_basis_value, grid_weight) {}

    /**
     * Compute the Losc curvature version 2 matrix.
     *
     * @ return: the Losc curvature matrix.
     */
    virtual SharedMatrix compute() override;
};

}   // namespace losc

#endif // _LOSC_SRC_CURVATURE_H_
