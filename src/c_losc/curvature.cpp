#include "matrix_impl.hpp"
#include <c_losc/curvature.h>
#include <losc/curvature.hpp>

// To use C++ class defined in a namespace in C code, we need to do typecast
// between the C struct and C++ class. You can use C-style cast of pointers or
// use C++ `reinterpret_cast` template function. See some examples here
// https://stackoverflow.com/questions/18945419/using-c-with-namespace-in-c

extern "C" {
//**********************************************
// ==> Binding `losc::DFAInfo` methods.
//**********************************************
LoscDFAInfo *losc_dfa_info_create(double gga_x, double hf_x, const char *name)
{
    return reinterpret_cast<LoscDFAInfo *>(
        new losc::DFAInfo(gga_x, hf_x, name));
}

void _losc_dfa_info_free(LoscDFAInfo **pptr_self)
{
    if (*pptr_self)
        delete reinterpret_cast<losc::DFAInfo *>(*pptr_self);
}

//**********************************************
// ==> Binding `losc::CurvatureBase` methods.
//**********************************************

static size_t losc_curvature_base_nlo(const LoscCurvatureBase *self)
{
    auto base = reinterpret_cast<losc::CurvatureBase *>(self->_p_base);
    return base->nlo();
}

static size_t losc_curvature_base_nfitbasis(const LoscCurvatureBase *self)
{
    auto base = reinterpret_cast<losc::CurvatureBase *>(self->_p_base);
    return base->nfitbasis();
}

static size_t losc_curvature_base_npts(const LoscCurvatureBase *self)
{
    auto base = reinterpret_cast<losc::CurvatureBase *>(self->_p_base);
    return base->npts();
}

static void losc_curvature_base_kappa(const LoscCurvatureBase *self,
                                      losc_matrix *K)
{
    auto base = reinterpret_cast<losc::CurvatureBase *>(self->_p_base);
    base->C_API_kappa(losc_matrix_to_eigen(K));
}

static LoscCurvatureBase *
create_losc_curvature_base(losc::CurvatureBase *_p_base)
{
    auto p_base = new LoscCurvatureBase;
    p_base->_p_base = reinterpret_cast<_LoscCurvatureBase *>(_p_base);
    p_base->nlo = losc_curvature_base_nlo;
    p_base->nfitbasis = losc_curvature_base_nfitbasis;
    p_base->npts = losc_curvature_base_npts;
    p_base->kappa = losc_curvature_base_kappa;
    return p_base;
}

//**********************************************
// ==> Binding `losc::CurvatureV1` methods.
//**********************************************
static void losc_curvature_v1_set_tau(const LoscCurvatureV1 *self, double tau)
{
    auto v1 = reinterpret_cast<losc::CurvatureV1*>(self->_p_v1);
    v1->set_tau(tau);
}

LoscCurvatureV1 *losc_curvature_v1_create(const LoscDFAInfo *dfa_info,
                                          const losc_matrix *df_pii,
                                          const losc_matrix *df_Vpq_inv,
                                          const losc_matrix *grid_lo,
                                          const double *grid_weight)
{
    auto cpp_dfa = reinterpret_cast<const losc::DFAInfo *>(dfa_info);
    const size_t npts = grid_lo->row_;
    // create a real `losc::CurvatureV1` object.
    losc::CurvatureV1 *p_v1 = new losc::CurvatureV1(
        *cpp_dfa, losc_matrix_to_eigen_const(df_pii),
        losc_matrix_to_eigen_const(df_Vpq_inv),
        losc_matrix_to_eigen_const(grid_lo),
        Eigen::Map<const Eigen::VectorXd>(grid_weight, npts));
    // create a `LoscCurvatureV1` struct that holds a bunch of function
    // pointers.
    auto rst = new LoscCurvatureV1();
    // Assign member variables.
    // type cast C++ class into C interface.
    rst->_p_v1 = reinterpret_cast<_LoscCurvatureV1 *>(p_v1);
    rst->p_base = create_losc_curvature_base(p_v1);
    // Assign new function pointers here in the future, if we want to export
    // new functions in `losc::CurvatureV1`.
    rst->set_tau = losc_curvature_v1_set_tau;
    return rst;
}

void *_losc_curvature_v1_free(LoscCurvatureV1 **pptr_self)
{
    if (*pptr_self) {
        auto self = (*pptr_self);
        // free `LoscCurvatureBase` struct.
        delete self->p_base;
        // free `losc::CurvatureV1` class.
        delete reinterpret_cast<losc::CurvatureV1 *>(self->_p_v1);
        // free `LoscCurvatureV1` struct.
        delete self;
        // set self pointer to null.
        *pptr_self = nullptr;
    }
}

//**********************************************
// ==> Binding `losc::CurvatureV2` methods.
//**********************************************
static void losc_curvature_v2_set_tau(const LoscCurvatureV2 *self, double tau)
{
    auto v2 = reinterpret_cast<losc::CurvatureV2*>(self->_p_v2);
    v2->set_tau(tau);
}

static void losc_curvature_v2_set_zeta(const LoscCurvatureV2 *self, double zeta)
{
    auto v2 = reinterpret_cast<losc::CurvatureV2*>(self->_p_v2);
    v2->set_zeta(zeta);
}

LoscCurvatureV2 *losc_curvature_v2_create(const LoscDFAInfo *dfa_info,
                                          const losc_matrix *df_pii,
                                          const losc_matrix *df_Vpq_inv,
                                          const losc_matrix *grid_lo,
                                          const double *grid_weight)
{
    auto cpp_dfa = reinterpret_cast<const losc::DFAInfo *>(dfa_info);
    const size_t npts = grid_lo->row_;
    // create a real `losc::CurvatureV2` object.
    losc::CurvatureV2 *p_v2 = new losc::CurvatureV2(
        *cpp_dfa, losc_matrix_to_eigen_const(df_pii),
        losc_matrix_to_eigen_const(df_Vpq_inv),
        losc_matrix_to_eigen_const(grid_lo),
        Eigen::Map<const Eigen::VectorXd>(grid_weight, npts));
    // create a `LoscCurvatureV2` struct that holds a bunch of function
    // pointers.
    auto rst = new LoscCurvatureV2();
    // Assign member variables.
    // type cast C++ class into C interface.
    rst->_p_v2 = reinterpret_cast<_LoscCurvatureV2 *>(p_v2);
    rst->p_base = create_losc_curvature_base(p_v2);
    // Assign new function pointers here in the future, if we want to export
    // new functions in `losc::CurvatureV2`.
    rst->set_tau = losc_curvature_v2_set_tau;
    rst->set_zeta = losc_curvature_v2_set_zeta;
    return rst;
}

void *_losc_curvature_v2_free(LoscCurvatureV2 **pptr_self)
{
    if (*pptr_self) {
        auto self = (*pptr_self);
        // free `LoscCurvatureBase` struct.
        delete self->p_base;
        // free `losc::CurvatureV2` class.
        delete reinterpret_cast<losc::CurvatureV2 *>(self->_p_v2);
        // free `LoscCurvatureV1` struct.
        delete self;
        // set self pointer to null.
        *pptr_self = nullptr;
    }
}
}
