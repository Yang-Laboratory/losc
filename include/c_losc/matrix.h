/**
 * @file matrix.h
 * @brief The C interface of matrix.
 *
 * The `struct losc_matrix` is used to represent matrix object and shared
 * between the C code and C++ LOSC library. The C users only need interface
 * as provided in this header file to use `struct losc_matrix`.
 * The implementation of `struct losc_matrix` is hid.
 */

#ifndef _LOSC_INCLUDE_C_LOSC_MATRIX_H_
#define _LOSC_INCLUDE_C_LOSC_MATRIX_H_
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct losc_matrix losc_matirx;

/**
 * Create a LOSC matrix object with data allocated.
 *
 * @param row: number of rows.
 * @param col: number of columns.
 * @return a losc_matrix pointer.
 */
losc_matrix *losc_matrix_create(size_t row, size_t col);

/**
 * Create a LOSC matrix object with data owned externally.
 *
 * @param row: number of rows.
 * @param col: number of columns.
 * @param data: points to an array which represents a [row, col] matrix.
 * The matrix is assumed to be stored in column-major. The size of the
 * array is assumed to be at least row x col.
 * @return a losc_matrix pointer.
 */
losc_matrix *losc_matrix_create_from_data(size_t row, size_t col, double *data);

/**
 * Release the resource associated with a losc_matrix pointer and
 * assign it NULL.
 *
 * @param pptr_m: points to a losc_matrix pointer.
 */
void _losc_matrix_free(losc_matrix **pptr_m);
#define losc_matrix_free(ptr_m) _losc_matrix_free(&(ptr_m))

/**
 * Return the pointer that points to the head of matrix data.
 */
double *losc_matrix_data(losc_matrix *m);

/**
 * Return the const pointer that points to the head of matrix data.
 */
const double *losc_matrix_data_const(const losc_matrix *m);

/**
 * Return the number of rows.
 */
size_t losc_matrix_rows(const losc_matrix *m);

/**
 * Return the number of columns.
 */
size_t losc_matrix_cols(const losc_matrix *m);

#ifdef __cplusplus
}
#endif

#endif
