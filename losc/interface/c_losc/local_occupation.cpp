#include "local_occupation.h"
#include "../../src/local_occupation.h"
#include "internal/matrix_impl.h"
#include <stdlib.h>

extern "C" {

losc_matrix *losc_local_occupation(const losc_matrix *C_lo,
                                   const losc_matrix *S, const losc_matrix *D)
{
    const size_t nlo = C_lo->col_;
    losc_matrix *LocalOcc = losc_matrix_create(nlo, nlo);
    losc::_local_occupation(losc_matrix_to_eigen_const(C_lo),
                            losc_matrix_to_eigen_const(S),
                            losc_matrix_to_eigen_const(D));
    return LocalOcc;
}
}
