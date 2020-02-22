/**
 * @file
 * @brief declarations relate to matrix I/O.
 */
#ifndef _LOSC_TEST_MATRIX_IO_H_
#define _LOSC_TEST_MATRIX_IO_H_

#include <losc/losc.h>
#include <memory>
#include <vector>
#include <string>

namespace test {

using losc::Matrix;
using std::vector;
using std::string;

/**
 * @brief Write a vector of matrices into a txt file.
 * @details The matrices are write in order. For each matrix, the format is the
 * following. The first line constaints the number of row and columns of the
 * matrix. This line is garranted to be the format
 * `printf("Dimension,%zu,%zu\n", row, col) and it is used to indicate the
 * separation of different matrices. The following multiples lines contains all
 * the matrix elements. Each matrix elements are comma-separated. The matrix
 * elements are written in column-wise order.
 *
 * @param [in] Mat: The set of input matrices to be written.
 * @param [in] fname: The name of the output txt file.
 * @param [in] num_per_line: number of matrix elements per line.
 */
void write_matrices_to_txt(vector<std::shared_ptr<Matrix>> &Mat,
                           const string &fname, size_t num_per_line = 5);

/**
 * @brief Read a vector of matrices from a txt file generated from
 * matrix::write_matrices_to_txt().
 *
 * @param [in] fname: the name (relative or absolute path) of the matrix txt
 * file.
 * @return std::vector<std::shared_ptr<Matrix>>: a vector of matrices.
 * @see matrix::write_matrices_to_txt()
 */
std::vector<std::shared_ptr<Matrix>> read_matrices_from_txt(const string &fname);

} // namespace matrix

#endif // _LOSC_TEST_MATRIX_IO_H_
