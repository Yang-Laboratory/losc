#ifndef _LOSC_CURVATURE_H_
#define _LOSC_CURVATURE_H_

#include <matrix/matrix.h>
#include <memory>
#include <vector>

namespace losc {

using matrix::Matrix;
using SharedMatrix = std::shared_ptr<Matrix>;
using std::vector;

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

    /**
     * LO coefficient matrix under AO.
     * Dimension: nlo x nbasis.
     * It can be provided by `localization::LocalizerBase::compute()` function.
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

    SharedMatrix kappa_J_;
    SharedMatrix kappa_xc_;

    void compute_kappa_J();

    public:
    CurvatureV1(SharedMatrix C_lo, SharedMatrix df_pmn, SharedMatrix df_Vpq_inverse);
    virtual void compute() override;
};

}

#endif // _LOSC_CURVATURE_H_
