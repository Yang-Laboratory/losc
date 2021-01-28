#include "matrix.h"
#include "matrix_internal.h"

// ==> export to C functions.
extern "C" {
losc_matrix *losc_matrix_create(size_t row, size_t col)
{
    return new losc_matirx(row, col);
}

losc_matrix *losc_matrix_create_from_data(size_t row, size_t col, double *data)
{
    return new losc_matirx(row, col, data);
}

void _losc_matrix_free(losc_matrix **pptr_m)
{
    if (pptr_m) {
        delete (*pptr_m);
    }
}
}
