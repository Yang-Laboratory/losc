/**
 * @file
 * @brief Definition matrix class.
 */

#include "matrix_helper.hpp"

#include <cmath>
#include <iomanip> // std::setprecision
#include <ios>
#include <iostream>
#include <losc/eigen_def.hpp>
#include <losc/exception.hpp>
#include <random>

using namespace losc;

namespace test {

static std::mt19937 g_rand_generator_mt19937_seed_fixed(1);

/**
 * @details The matrix elements are printed in the format `%16.8e`.
 */
void mtx_show_full(const LOSCMatrix &m, std::ostream &os,
                   size_t elements_per_line)
{
    size_t row = m.rows();
    size_t col = m.cols();
    std::ios_base::fmtflags os_flags(os.flags());
    os << std::fixed;
    os << "dimension: " << row << " x " << col << ", showing the full matrix."
       << std::endl;
    for (size_t i = 0; i < row; ++i) {
        os << std::setw(5) << i + 1 << ":\n";
        os << std::scientific;
        for (int j = 1; j <= col; ++j) {
            os << std::setprecision(8) << std::setw(16) << m(i, j - 1);
            if (j != col)
                os << ",";
            if (j % elements_per_line == 0)
                os << std::endl;
        }
        os << std::endl;
    }
    os << std::endl;
    os.flush();
    os.flags(os_flags);
}

/**
 * @details The matrix elements are printed in the format `%16.8e`.
 */
void mtx_show_lower(const LOSCMatrix &m, std::ostream &os,
                    size_t elements_per_line)
{
    std::ios_base::fmtflags os_flags(os.flags());
    os << std::fixed;
    os << "dimension: " << m.rows() << " x " << m.cols()
       << ", showing the lower triangular parts." << std::endl;
    for (size_t i = 0; i < m.rows(); i++) {
        os << std::setw(5) << i + 1 << ":\n";
        os << std::scientific;
        for (int j = 0; j <= i; j++) {
            os << std::setprecision(8) << std::setw(15) << m(i, j);
            if (j != i)
                os << ",";
            if ((j + 1) % elements_per_line == 0)
                os << std::endl;
        }
        os << std::endl;
    }
    os << std::endl;
    os.flush();
    os.flags(os_flags);
}

bool mtx_is_symmetric(const LOSCMatrix &m, double threshold)
{
    threshold = std::abs(threshold);
    if (!mtx_is_square(m)) {
        return false;
    }
    for (size_t i = 0; i < m.rows(); i++) {
        for (size_t j = 0; j < i; j++) {
            if (std::abs(m(i, j) - m(j, i)) > threshold) {
                return false;
            }
        }
    }
    return true;
}

/**
 * @note The dimension of two matrices will be compared as well. If not matched,
 * it will return false.
 */
bool mtx_is_cwise_equal(const LOSCMatrix &m, const LOSCMatrix &other,
                        double threshold)
{
    if (!mtx_is_same_dimension_to(m, other)) {
        return false;
    }
    LOSCMatrix diff = m - other;
    return diff.isZero(threshold);
}

bool mtx_is_same_dimension_to(const LOSCMatrix &m, const LOSCMatrix &other)
{
    return ((m.rows() == other.rows()) && (m.cols() == other.cols()));
}

void mtx_to_symmetric(LOSCMatrix &m, const string &uplo)
{
    const size_t row = m.rows();
    const size_t col = m.cols();
    if (!mtx_is_square(m)) {
        throw losc::exception::DimensionError(
            "Cannot symmetrize a matrix that is not squared.");
    }
    if (uplo == "U") {
        for (size_t i = 0; i < row; i++) {
            for (size_t j = 0; j < i; j++) {
                m(i, j) = m(j, i);
            }
        }
    } else if (uplo == "L") {
        for (size_t i = 0; i < row; i++) {
            for (size_t j = 0; j < i; j++) {
                m(j, i) = m(i, j);
            }
        }
    }
}

/**
 * @details This function is implemented with `std::random_device` to generate
 * random seed to `std::mt1993` radom number generator, and use
 * `std::uniform_real_distribution` to keep uniform distributed.
 *
 * @see Eigen::setRandom() can also set a matrix element normally distributed in
 * range [-1 : 1].
 */
void mtx_randomize(LOSCMatrix &m, double a, double b)
{
    std::random_device
        rd; // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(a, b);
    for (size_t i = 0; i < m.rows(); i++) {
        for (size_t j = 0; j < m.cols(); j++) {
            m(i, j) = dis(gen);
        }
    }
}

/**
 * @details `std::mt1993` seeded with constant value 1 is used as the radom
 * number generator, and use `std::uniform_real_distribution` to keep uniform
 * distributed.
 *
 * @see Eigen::setRandom() can also set a matrix element normally distributed in
 * range [-1 : 1].
 */
void mtx_randomize_seed_fixed(LOSCMatrix &m, double a, double b)
{
    std::uniform_real_distribution<> dis(a, b);
    for (size_t i = 0; i < m.rows(); i++) {
        for (size_t j = 0; j < m.cols(); j++) {
            m(i, j) = dis(g_rand_generator_mt19937_seed_fixed);
        }
    }
}

} // namespace test
