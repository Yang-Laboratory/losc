/**
 * @file
 * @brief Definition matrix class.
 */

#include "matrix.h"

#include "exception.h"
#include <cmath>
#include <iomanip> // std::setprecision
#include <ios>
#include <iostream>
#include <random>

namespace losc {

static std::mt19937 g_rand_generator_mt19937_seed_fixed(1);

/**
 * @details The matrix elements are printed in the format `%16.8e`.
 */
void Matrix::show_full(std::ostream &os, size_t elements_per_line) const
{
    size_t row = this->rows();
    size_t col = this->cols();
    const Matrix &self = *this;
    std::ios_base::fmtflags os_flags(os.flags());
    os << std::fixed;
    os << "dimension: " << row << " x " << col << ", showing the full matrix."
       << std::endl;
    for (size_t i = 0; i < row; ++i) {
        os << std::setw(5) << i + 1 << ":\n";
        os << std::scientific;
        for (int j = 1; j <= col; ++j) {
            os << std::setprecision(8) << std::setw(16) << self(i, j - 1);
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
void Matrix::show_lower(std::ostream &os, size_t elements_per_line) const
{
    size_t row = this->rows();
    size_t col = this->cols();
    const Matrix &self = *this;
    std::ios_base::fmtflags os_flags(os.flags());
    os << std::fixed;
    os << "dimension: " << row << " x " << col
       << ", showing the lower triangular parts." << std::endl;
    for (size_t i = 0; i < row; i++) {
        os << std::setw(5) << i + 1 << ":\n";
        os << std::scientific;
        for (int j = 0; j <= i; j++) {
            os << std::setprecision(8) << std::setw(15) << self(i, j);
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

bool Matrix::is_symmetric(double threshold) const
{
    const Matrix &self = *this;
    threshold = std::abs(threshold);
    if (!self.is_square()) {
        return false;
    }
    for (size_t i = 0; i < self.rows(); i++) {
        for (size_t j = 0; j < i; j++) {
            if (std::abs(self(i, j) - self(j, i)) > threshold) {
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
bool Matrix::is_cwise_equal(const Matrix &other, double threshold) const
{
    const Matrix &self = *this;
    if (!self.is_same_dimension_to(other)) {
        return false;
    }
    Matrix diff = *this - other;
    return diff.isZero(threshold);
}

bool Matrix::is_same_dimension_to(const Matrix &other) const
{
    return ((this->rows() == other.rows()) && (this->cols() == other.cols()));
}

void Matrix::to_symmetric(const string &uplo)
{
    Matrix &self = *this;
    const size_t row = self.rows();
    const size_t col = self.cols();
    if (!self.is_square()) {
        throw exception::DimensionError(
            "Cannot symmetrize a matrix that is not squared.");
    }
    if (uplo == "U") {
        for (size_t i = 0; i < row; i++) {
            for (size_t j = 0; j < i; j++) {
                self(i, j) = self(j, i);
            }
        }
    } else if (uplo == "L") {
        for (size_t i = 0; i < row; i++) {
            for (size_t j = 0; j < i; j++) {
                self(j, i) = self(i, j);
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
void Matrix::randomize(double a, double b)
{
    std::random_device
        rd; // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(a, b);
    for (size_t i = 0; i < this->rows(); i++) {
        for (size_t j = 0; j < this->cols(); j++) {
            (*this)(i, j) = dis(gen);
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
void Matrix::randomize_seed_fixed(double a, double b)
{
    std::uniform_real_distribution<> dis(a, b);
    for (size_t i = 0; i < this->rows(); i++) {
        for (size_t j = 0; j < this->cols(); j++) {
            (*this)(i, j) = dis(g_rand_generator_mt19937_seed_fixed);
        }
    }
}

} // namespace losc
