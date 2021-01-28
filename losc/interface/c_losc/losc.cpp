#include "losc.h"
#include "../../src/curvature.h"
#include "matrix_internal.h"
#include <Eigen/Dense>

namespace losc {

//***********************
// ==> Class DFAInfo
//***********************

DFAInfo *losc_DFAInfo_create(double gga_x_wt, double hf_x_wt, const char *name)
{
    return new DFAInfo(gga_x_wt, hf_x_wt, name);
}

void _losc_DFAInfo_free(DFAInfo **pptr_m)
{
    if (*pptr_m)
        delete (*pptr_m);
}

//***********************
// ==> Class CurvatureV1
//***********************

CurvatureV1 *losc_CurvatureV1_create(const DFAInfo *dfa_info,
                                     const losc_matrix *C_lo,
                                     const losc_matrix *df_pii,
                                     const losc_matrix *df_Vpq_inv,
                                     const losc_matrix *grid_basis_value,
                                     const double *grid_weight)
{
    const size_t npts = grid_basis_value->row_;
    return new CurvatureV1(*dfa_info, losc_matrix_to_eigen_const(C_lo),
                           losc_matrix_to_eigen_const(df_pii),
                           losc_matrix_to_eigen_const(df_Vpq_inv),
                           losc_matrix_to_eigen_const(grid_basis_value),
                           Eigen::Map<const VectorXd>(grid_weight, npts));
}

losc_matrix *losc_CurvatureV1_kappa(const CurvatureV1 *cur) {
    auto kappa = losc_matrix_create(cur->nlo(), cur->nlo());
    cur->kappa(losc_matrix_to_eigen(kappa));
    return kappa;
}

//***********************
// ==> Class CurvatureV2
//***********************

CurvatureV2 *losc_CurvatureV2_create(const DFAInfo *dfa_info,
                                     const losc_matrix *C_lo,
                                     const losc_matrix *df_pii,
                                     const losc_matrix *df_Vpq_inv,
                                     const losc_matrix *grid_basis_value,
                                     const double *grid_weight)
{
    const size_t npts = grid_basis_value->row_;
    return new CurvatureV2(*dfa_info, losc_matrix_to_eigen_const(C_lo),
                           losc_matrix_to_eigen_const(df_pii),
                           losc_matrix_to_eigen_const(df_Vpq_inv),
                           losc_matrix_to_eigen_const(grid_basis_value),
                           Eigen::Map<const VectorXd>(grid_weight, npts));
}

losc_matrix *losc_CurvatureV2_kappa(const CurvatureV2 *cur) {
    auto kappa = losc_matrix_create(cur->nlo(), cur->nlo());
    cur->kappa(losc_matrix_to_eigen(kappa));
    return kappa;
}

} // namespace losc
