#include "localization.h"
#include "../../src/localization.h"
#include "internal/matrix_impl.h"

extern "C" {

typedef losc::LocalizerBase _LocalizerBase;
typedef losc::LoscLocalizerV2 _LoscLocalizerV2;

//**********************************************
// ==> Binding `losc::LocalizerBase` methods.
//**********************************************

static void losc_localizer_base_set_max_iter(const LoscLocalizerBase *self,
                                             size_t max_iter)
{
    self->_p_base->set_max_iter(max_iter);
}

static void losc_localizer_base_set_convergence(const LoscLocalizerBase *self,
                                                double tol)
{
    self->_p_base->set_convergence(tol);
}

static void
losc_localizer_base_set_random_permutation(const LoscLocalizerBase *self,
                                           bool flag)
{
    self->_p_base->set_random_permutation(flag);
}

static void losc_localizer_base_lo_U(const LoscLocalizerBase *self,
                                     losc_matrix *L, losc_matrix *U)
{
    self->_p_base->lo_U(losc_matrix_to_eigen(L), losc_matrix_to_eigen(U));
}

/**
 * Create a `struct LoscLocalizerBase` and assign all its function pointers
 */
static LoscLocalizerBase *create_losc_localizer_base(_LocalizerBase *_p_base)
{
    auto p_base = new LoscLocalizerBase;
    p_base->_p_base = _p_base;
    p_base->set_max_iter = losc_localizer_base_set_max_iter;
    p_base->set_convergence = losc_localizer_base_set_convergence;
    p_base->set_random_permutation = losc_localizer_base_set_random_permutation;
    p_base->lo_U = losc_localizer_base_lo_U;
    return p_base;
}

//**********************************************
// ==> Binding `losc::LoscLocalizerV2` methods.
//**********************************************
LoscLocalizerV2 *losc_localizer_v2_create(const losc_matrix *C_lo_basis,
                                          const losc_matrix *H_ao,
                                          const losc_matrix *D_ao[3])
{
    using losc::vector;
    using losc::RefConstMat;
    // create a real `losc::LoscLocalizerV2` object.
    vector<RefConstMat> D_ao_ref;
    for (size_t i = 0; i < 3; ++i) {
        D_ao_ref.push_back(losc_matrix_to_eigen_const(D_ao[i]));
    }
    _LoscLocalizerV2 *_p_v2 =
        new _LoscLocalizerV2(losc_matrix_to_eigen_const(C_lo_basis),
                             losc_matrix_to_eigen_const(H_ao), D_ao_ref);

    // create a `LoscLocalizerV2` struct that holds a bunch of function
    // pointers.
    auto rst = new LoscLocalizerV2;
    // Assign member variables.
    rst->_p_v2 = _p_v2;
    rst->p_base = create_losc_localizer_base(_p_v2);

    // Assign new function pointers here in the future, if we want to export
    // new functions in `losc::CurvatureV1`.
    return rst;
}

void *_losc_localizer_v2_free(LoscLocalizerV2 **pptr_self)
{
    if (*pptr_self) {
        auto self = (*pptr_self);
        // free `LoscLocalizerBase` struct.
        delete self->p_base;
        // free `losc::LocalizerV2` class.
        delete self->_p_v2;
        // free `LoscLocalizerV2` struct.
        delete self;
        // set self pointer to null.
        *pptr_self = nullptr;
    }
}
}
