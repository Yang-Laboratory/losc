/**
 * @file
 * @brief declearation for blas library.
 */
#ifndef _LOSC_SRC_BLAS_BASE_H_
#define _LOSC_SRC_BLAS_BASE_H_

namespace losc {

namespace blas {

static int izero[] = {0};
static int ione[] = {1};
static double dzero[] = {0.0};
static double done[] = {1.0};

extern "C" void drot_(const int *N, double *x, const int *incx, double *y,
                      const int *incy, const double *c, const double *s);

} // namespace blas
} // namespace losc

#endif // _LOSC_SRC_BLAS_BASE_H_
