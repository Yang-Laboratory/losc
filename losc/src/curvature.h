#ifndef _LOSC_CURVATURE_H_
#define _LOSC_CURVATURE_H_

#include <matrix/matrix.h>
#include <memory>
#include <vector>

namespace losc {

using std::vector;
using matrix::Matrix;
using SharedMatrix = std::shared_ptr<Matrix>;
using SharedDoubleVector = std::shared_ptr<vector<double>>;

class CurvatureBase {
    protected:
    SharedMatrix kappa_;

    public:
    virtual void compute() = 0;
    SharedMatrix get_curvature() {return kappa_;}
};

class CurvatureV1 : public CurvatureBase {
    private:
    size_t nlo_;
    size_t nbasis_;
    size_t nfitbasis_;
    size_t npts_;

    double para_cx_ = 0.930526;
    double para_exf_ = 1.2378;
    double para_alpha_ = 0.0;
    double para_beta_ = 0.0;

    /**
     * LO coefficient matrix under AO.
     * Dimension: nlo x nbasis.
     * It can be provided by `losc::LocalizerBase::compute()` function.
     */
    SharedMatrix C_lo_;

    /**
     * three-body integral of <p|mn>.
     * m, n is the AO basis index, and p is the fitbasis index.
     * Dimension: [nfitbasis, nbasis * (nbasis + 1)/2 ].
     */
    SharedMatrix df_pmn_;

    /**
     * Inverse of Vpq matrix.
     * Dimension: [nfitbasis, nfitbasis].
     */
    SharedMatrix df_Vpq_inverse_;

    /**
     * AO basis value on grid.
     * dimension: [npts, nbasis].
     */
    SharedMatrix grid_basis_value_;

    /**
     * grid weight.
     * dimension: npts.
     */
    SharedDoubleVector grid_weight_;

    SharedMatrix kappa_J_;
    SharedMatrix kappa_xc_;

    void compute_kappa_J();

    void compute_kappa_xc();

    public:
    CurvatureV1(SharedMatrix C_lo, SharedMatrix df_pmn, SharedMatrix df_Vpq_inverse,
                SharedMatrix grid_basis_value, SharedDoubleVector grid_weight);
    virtual void compute() override;
};

class CurvatureV2 : public CurvatureBase {
    private:
    size_t nlo_;
    size_t nbasis_;
    size_t nfitbasis_;
    size_t npts_;

    double para_tau_ = 8.0;
    double para_cx_ = 0.930526;
    double para_exf_ = 1.2378;
    double para_alpha_ = 0.0;
    double para_beta_ = 0.0;


    /**
     * LO coefficient matrix under AO.
     * Dimension: nlo x nbasis.
     * It can be provided by `losc::LocalizerBase::compute()` function.
     */
    SharedMatrix C_lo_;

    /**
     * three-body integral of <p|mn>.
     * m, n is the AO basis index, and p is the fitbasis index.
     * Dimension: [nfitbasis, nbasis * (nbasis + 1)/2 ].
     */
    SharedMatrix df_pmn_;

    /**
     * Inverse of Vpq matrix.
     * Dimension: [nfitbasis, nfitbasis].
     */
    SharedMatrix df_Vpq_inverse_;

    /**
     * AO basis value on grid.
     * dimension: [npts, nbasis].
     */
    SharedMatrix grid_basis_value_;

    /**
     * grid weight.
     * dimension: npts.
     */
    SharedDoubleVector grid_weight_;


    public:
    CurvatureV2(SharedMatrix C_lo, SharedMatrix df_pmn, SharedMatrix df_Vpq_inverse,
                SharedMatrix grid_basis_value, SharedDoubleVector grid_weight);

    virtual void compute() override;
};

}

#endif // _LOSC_CURVATURE_H_
