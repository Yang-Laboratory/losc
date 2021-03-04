/**
 * @file
 * @brief The C interface of matrix used in LOSC.
 * @details The struct losc_matrix is used to represent matrix object and
 * shared between the C code and C++ LOSC library. The C users only need
 * the interface as provided in this header file to use struct losc_matrix.
 * The implementation of struct losc_matrix is hid.
 */

#ifndef _LOSC_INCLUDE_C_LOSC_MATRIX_H_
#define _LOSC_INCLUDE_C_LOSC_MATRIX_H_
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief losc_matrix
 * @note The complete definition of losc_matrix is hid.
 */
typedef struct losc_matrix losc_matirx;

/**
 * @brief Create a LOSC matrix object with data allocated.
 * @param [in] row number of rows.
 * @param [in] col number of columns.
 * @return Create a matrix and return the pointer.
 */
losc_matrix *losc_matrix_create(size_t row, size_t col);

/**
 * @brief Create a LOSC matrix object with data owned externally.
 * @param [in] row: number of rows.
 * @param [in] col: number of columns.
 * @param [in] data: `data` points to an array which represents a matrix with
 * dimension of `[row, col]`. The matrix is assumed to be stored in row-major.
 * The size of the array is assumed to be at least `row x col`.
 * @return Create a matrix and return the pointer.
 * @note This created matrix does not own the data.
 */
losc_matrix *losc_matrix_create_from_data(size_t row, size_t col, double *data);

/**
 * @brief Free the losc_matrix.
 * @param [in, out] ptr a pointer with type `losc_matrix *`. At exit,
 * `ptr` is set to null.
 */
#define losc_matrix_free(ptr) _losc_matrix_free(&(ptr))
void _losc_matrix_free(losc_matrix **pptr_m);

/**
 * @brief Return the pointer that points to the head of matrix data.
 * @param [in] m a losc_matrix pointer.
 */
double *losc_matrix_data(losc_matrix *m);

/**
 * @brief Return the const pointer that points to the head of matrix data.
 * @param [in] m a const losc_matrix pointer.
 */
const double *losc_matrix_data_const(const losc_matrix *m);

/**
 * @brief Return the number of rows.
 */
size_t losc_matrix_rows(const losc_matrix *m);

/**
 * @brief Return the number of columns.
 */
size_t losc_matrix_cols(const losc_matrix *m);

#ifdef __cplusplus
}
#endif

#endif
