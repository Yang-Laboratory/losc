/**
 * @file matrix.cpp
 * @brief The implementation of matrix interface with C.
 *
 * This is the implementation of the matrix interface, not the matrix itself.
 * The implementation of matrix is hid internally and not exposed to the
 * C users intentionally.
 */
#include "matrix.h"
#include "internal/matrix_impl.h"

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

double *losc_matrix_data(losc_matrix *m) { return m->data_; }

const double *losc_matrix_data_const(const losc_matrix *m) { return m->data_; }

size_t losc_matrix_rows(const losc_matrix *m) { return m->row_; }

size_t losc_matrix_cols(const losc_matrix *m) { return m->col_; }
}
