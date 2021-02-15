#ifndef _LOSC_TESTS_MATRIX_HELPER_HPP_
#define _LOSC_TESTS_MATRIX_HELPER_HPP_

#include <Eigen/Dense>
#include <iostream>
#include <losc/eigen_def.hpp>
#include <ostream>
#include <string>

namespace test {

using losc::LOSCMatrix;
using std::string;

/**
 * @brief Print out the full matrix to a `std::ostream`.
 * @param [in] os: The output stream object to print the matrix. Default
 * `os = std::cout`.
 * @param [in] elements_per_line: number of elements being printed per line.
 *  Default `element_per_line = 5`.
 *
 * @note This function will not change the settings of the `os`.
 * @note The format of printing matrix in this matrix is fixed. If you want
 * print the matrix in a more flexible way, use the Eige::LOSCMatrix
 * interface. Eigen library provides the overloaded `operator <<` to insert
 * an Eigen matrix into an `std::ostream`. It also provide `Eigen::IOFormat`
 * class to configure the output matrix format.
 */
void mtx_show_full(const LOSCMatrix &m, std::ostream &os = std::cout,
                   size_t elements_per_line = 5);

/**
 * @brief Print out the lower triangular matrix (including diagonal
 * elements).
 * @param [in] os: The output stream object to print the matrix. Default
 * `os = std::cout`.
 * @param [in] elements_per_line: number of elements being printed per line.
 *  Default `elements_per_line = 5`.
 *
 * @note This function will not change the settings of the `os`.
 * @see losc::show_full() for more format information.
 */
void mtx_show_lower(const LOSCMatrix &m, std::ostream &os = std::cout,
                    size_t elements_per_line = 5);

/**
 * @brief Check if the matrix is square or not.
 * @return bool
 */
inline bool mtx_is_square(const LOSCMatrix &m)
{
    return (m.rows() == m.cols());
}

/**
 * @brief Check if or not the matrix is numerically symmetric based on
 * the input threshold. Default `threshold = 1e-10`.
 *
 * @param [in] threshold: the threshold of testing numerical equality
 * between two float-point numbers.
 * @return bool.
 */
bool mtx_is_symmetric(const LOSCMatrix &m, double threshold = 1e-10);

/**
 * @brief Check if or not the matrix is numerically equal to the other
 * matrix by element-wise (coefficient-wise) comparision.
 * Default `threshold = 1e-10`.
 *
 * @param [in] other: The other matrix to be compared with.
 * @return bool
 *
 * @see Einen library also provides a similar function `Eigen::isApprox`.
 */
bool mtx_is_cwise_equal(const LOSCMatrix &m, const LOSCMatrix &other,
                        double threshold = 1e-10);

/**
 * @brief Check if or not the matrix has the same dimension to the other
 * matrix.
 *
 * @param [in] other: The other matrix to be compared with.
 * @return bool.
 */
bool mtx_is_same_dimension_to(const LOSCMatrix &m, const LOSCMatrix &other);

/**
 * @brief Make the matrix to be symmetric.
 *
 * @param [in] uplo: when \p uplo equals to "U", the upper triangular part
 * is used. when \p uplo equals to "L", the lower triangular part is used.
 */
void mtx_to_symmetric(LOSCMatrix &m, const string &uplo);

/**
 * @brief Make the matrix to be random with elements uniformly distributed
 * in range [a, b).
 *
 * @param [in] a: left range bound.
 * @param [in] b: right range bound.
 *
 * @note The random number generator is initialized with a NON-FIXED seed.
 * So the randomness behavior is not repeatable at running time.
 */
void mtx_randomize(LOSCMatrix &m, double a, double b);

/**
 * @brief Make the matrix to be random with elements uniformly distributed
 * in range [a, b).
 *
 * @param [in] a: left range bound.
 * @param [in] b: right range bound.
 * @return Matrix&: the randomized matrix itself.
 *
 * @note The random number generator is initialized with a FIXED seed. So
 * the randomness behavior is repeatable at running time.
 */
void mtx_randomize_seed_fixed(LOSCMatrix &m, double a, double b);

} // namespace test

#endif // _LOSC_SRC_MATRIX_H_
