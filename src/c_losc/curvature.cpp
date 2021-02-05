#include "matrix_impl.hpp"
#include <c_losc/curvature.h>
#include <losc/curvature.hpp>

extern "C" {

// Don't understand why. If I put the following line behind
// typedef losc::CurvatureBase losc_CurvatureBase, it will be a error.
// typedef struct losc_CurvatureBase losc_CurvatureBase;

// Here is the trick to deal with the namespace in C++. Using typedef
// to create an alias for the class with its full name (including namespace).
// Make sure the alias is the same as the name in C code. Therefore,
// The C code can access the C++ class within a namespace.
typedef losc::CurvatureBase _CurvatureBase;
typedef losc::CurvatureV1 _CurvatureV1;
typedef losc::CurvatureV2 _CurvatureV2;
typedef losc::DFAInfo LoscDFAInfo;

//**********************************************
// ==> Binding `losc::DFAInfo` methods.
//**********************************************
LoscDFAInfo *losc_dfa_info_create(double gga_x, double hf_x, const char *name)
{
    return new LoscDFAInfo(gga_x, hf_x, name);
}

void _losc_dfa_info_free(LoscDFAInfo **pptr_self)
{
    if (*pptr_self)
        delete (*pptr_self);
}

//**********************************************
// ==> Binding `losc::CurvatureBase` methods.
//**********************************************

static size_t losc_curvature_base_nlo(const LoscCurvatureBase *self)
{
    return self->_p_base->nlo();
}

static size_t losc_curvature_base_nfitbasis(const LoscCurvatureBase *self)
{
    return self->_p_base->nfitbasis();
}

static size_t losc_curvature_base_npts(const LoscCurvatureBase *self)
{
    return self->_p_base->npts();
}

static void losc_curvature_base_kappa(const LoscCurvatureBase *self,
                                      losc_matrix *K)
{
    self->_p_base->C_API_kappa(losc_matrix_to_eigen(K));
}

static LoscCurvatureBase *create_losc_curvature_base(_CurvatureBase *_p_base)
{
    auto p_base = new LoscCurvatureBase;
    p_base->_p_base = _p_base;
    p_base->nlo = losc_curvature_base_nlo;
    p_base->nfitbasis = losc_curvature_base_nfitbasis;
    p_base->npts = losc_curvature_base_npts;
    p_base->kappa = losc_curvature_base_kappa;
    return p_base;
}

//**********************************************
// ==> Binding `losc::CurvatureV1` methods.
//**********************************************
LoscCurvatureV1 *losc_curvature_v1_create(const LoscDFAInfo *dfa_info,
                                          const losc_matrix *df_pii,
                                          const losc_matrix *df_Vpq_inv,
                                          const losc_matrix *grid_lo,
                                          const double *grid_weight)
{
    const size_t npts = grid_lo->row_;
    // create a real `losc::CurvatureV1` object.
    _CurvatureV1 *_p_v1 =
        new _CurvatureV1(*dfa_info, losc_matrix_to_eigen_const(df_pii),
                         losc_matrix_to_eigen_const(df_Vpq_inv),
                         losc_matrix_to_eigen_const(grid_lo),
                         Eigen::Map<const Eigen::VectorXd>(grid_weight, npts));

    // create a `LoscCurvatureV1` struct that holds a bunch of function
    // pointers.
    auto rst = new LoscCurvatureV1();
    // Assign member variables.
    rst->_p_v1 = _p_v1;
    rst->p_base = create_losc_curvature_base(_p_v1);
    // Assign new function pointers here in the future, if we want to export
    // new functions in `losc::CurvatureV1`.
    return rst;
}

void *_losc_curvature_v1_free(LoscCurvatureV1 **pptr_self)
{
    if (*pptr_self) {
        auto self = (*pptr_self);
        // free `LoscCurvatureBase` struct.
        delete self->p_base;
        // free `losc::CurvatureV1` class.
        delete self->_p_v1;
        // free `LoscCurvatureV1` struct.
        delete self;
        // set self pointer to null.
        *pptr_self = nullptr;
    }
}

//**********************************************
// ==> Binding `losc::CurvatureV2` methods.
//**********************************************
LoscCurvatureV2 *losc_curvature_v2_create(const LoscDFAInfo *dfa_info,
                                          const losc_matrix *df_pii,
                                          const losc_matrix *df_Vpq_inv,
                                          const losc_matrix *grid_lo,
                                          const double *grid_weight)
{
    const size_t npts = grid_lo->row_;
    // create a real `losc::CurvatureV2` object.
    _CurvatureV2 *_p_v2 =
        new _CurvatureV2(*dfa_info, losc_matrix_to_eigen_const(df_pii),
                         losc_matrix_to_eigen_const(df_Vpq_inv),
                         losc_matrix_to_eigen_const(grid_lo),
                         Eigen::Map<const Eigen::VectorXd>(grid_weight, npts));

    // create a `LoscCurvatureV2` struct that holds a bunch of function
    // pointers.
    auto rst = new LoscCurvatureV2();
    // Assign member variables.
    rst->_p_v2 = _p_v2;
    rst->p_base = create_losc_curvature_base(_p_v2);
    // Assign new function pointers here in the future, if we want to export
    // new functions in `losc::CurvatureV2`.
    return rst;
}

void *_losc_curvature_v2_free(LoscCurvatureV2 **pptr_self)
{
    if (*pptr_self) {
        auto self = (*pptr_self);
        // free `LoscCurvatureBase` struct.
        delete self->p_base;
        // free `losc::CurvatureV2` class.
        delete self->_p_v2;
        // free `LoscCurvatureV1` struct.
        delete self;
        // set self pointer to null.
        *pptr_self = nullptr;
    }
}
}
