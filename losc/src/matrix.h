#ifndef _LOSC_SRC_MATRIX_H_
#define _LOSC_SRC_MATRIX_H_

#include <Eigen/Dense>
#include <iostream>
#include <ostream>
#include <string>

namespace losc {

using Eigen::MatrixXd;
using std::string;

/**
 * @brief Matrix class used in the Losc libary.
 * @details The Matrix class is publicly inherited from Eigen::MatrixXd,
 * includeing all the public member functions. So you can
 * construct and operate losc::Matrix the same way as the Eigen::MatrixXd. In
 * addition to the interface provided by Eigen::MatrixXd, losc::Matrix
 * also provides some simple functions for convenience.
 *
 * @note Since it is inherited from Eigen::MatrixXd, there is not choice to
 * set the matrix storage-order. Eigen::MatrixXd use the default column-major
 * order, so the losc::Matrix. Keep this in mind!
 */
class Matrix : public MatrixXd {
  public:
    using MatrixXd::MatrixXd;

    /**
     * @brief Print out the full matrix to a `std::ostream`.
     * @param [in] os: The output stream object to print the matrix. Default
     * `os = std::cout`.
     * @param [in] elements_per_line: number of elements being printed per line.
     *  Default `element_per_line = 5`.
     *
     * @note This function will not change the settings of the `os`.
     * @note The format of printing matrix in this matrix is fixed. If you want
     * print the matrix in a more flexible way, use the Eige::MatrixXd
     * interface. Eigen library provides the overloaded `operator <<` to insert
     * an Eigen matrix into an `std::ostream`. It also provide `Eigen::IOFormat`
     * class to configure the output matrix format.
     */
    void show_full(std::ostream &os = std::cout,
                   size_t elements_per_line = 5) const;

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
    void show_lower(std::ostream &os = std::cout,
                    size_t elements_per_line = 5) const;

    /**
     * @brief Check if the matrix is square or not.
     * @return bool
     */
    bool is_square() const { return (this->rows() == this->cols()); }

    /**
     * @brief Check if or not the matrix is numerically symmetric based on
     * the input threshold. Default `threshold = 1e-10`.
     *
     * @param [in] threshold: the threshold of testing numerical equality
     * between two float-point numbers.
     * @return bool.
     */
    bool is_symmetric(double threshold = 1e-10) const;

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
    bool is_cwise_equal(const Matrix &other, double threshold = 1e-10) const;

    /**
     * @brief Check if or not the matrix has the same dimension to the other
     * matrix.
     *
     * @param [in] other: The other matrix to be compared with.
     * @return bool.
     */
    bool is_same_dimension_to(const Matrix &other) const;

    /**
     * @brief Make the matrix to be symmetric.
     *
     * @param [in] uplo: when \p uplo equals to "U", the upper triangular part
     * is used. when \p uplo equals to "L", the lower triangular part is used.
     */
    void to_symmetric(const string &uplo);

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
    void randomize(double a, double b);

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
    void randomize_seed_fixed(double a, double b);
};

} // namespace losc

#endif // _LOSC_SRC_MATRIX_H_
