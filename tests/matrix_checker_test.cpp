#include "matrix_helper.hpp"
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <vector>

using Eigen::MatrixXd;
using std::vector;

using namespace test;

/**
 * Test is_symmetric() method.
 */
TEST(MatrixCheckerTest, is_symmetric_test)
{
    MatrixXd A = MatrixXd::Zero(10, 10);
    // zero matrix should be regarded as symmetric.
    EXPECT_TRUE(mtx_is_symmetric(A));
    EXPECT_TRUE(mtx_is_symmetric(A, 1e-16));
    EXPECT_TRUE(mtx_is_symmetric(A, 1e-300));

    mtx_randomize(A, 0, 1);
    EXPECT_FALSE(
        mtx_is_symmetric(A)); // random matrix usually is not symmetric.
    EXPECT_TRUE(mtx_is_symmetric(A, 2)); // lower threshold to be symmetric.
    mtx_to_symmetric(A, "L");
    EXPECT_TRUE(mtx_is_symmetric(A)); // real symmetric matrix.
    EXPECT_TRUE(mtx_is_symmetric(
        A, -1e-16)); // happended to set it to negative threshold.

    // if customized threshold working?
    double thred = 1e-8;
    A(2, 3) = A(3, 2) - thred;
    EXPECT_FALSE(mtx_is_symmetric(A, thred / 10));
    EXPECT_TRUE(mtx_is_symmetric(A, thred * 10));

    // non-square matrix.
    MatrixXd B(2, 3);
    EXPECT_FALSE(mtx_is_symmetric(B));
    EXPECT_FALSE(mtx_is_symmetric(B, 100));
    EXPECT_FALSE(mtx_is_symmetric(B, 1e-300));
}

/**
 * Test is_cwise_equal() method.
 */
TEST(MatrixCheckerTest, is_cwise_equal_test)
{
    MatrixXd A = MatrixXd::Zero(10, 8);
    MatrixXd B = MatrixXd::Zero(10, 8);
    // real zero matrix should be equal.
    EXPECT_TRUE(mtx_is_cwise_equal(A, B));
    EXPECT_TRUE(mtx_is_cwise_equal(A, B, 1e-16));
    EXPECT_TRUE(mtx_is_cwise_equal(A, B, 1e-300));

    // customized threshold
    double thred = 1e-8;
    A(1, 1) += thred;
    EXPECT_FALSE(mtx_is_cwise_equal(A, B, thred / 10));
    EXPECT_TRUE(mtx_is_cwise_equal(A, B, thred * 10));

    // random matrix usually is not equal to each other.
    mtx_randomize(A, 0, 1);
    mtx_randomize(B, 0, 1);
    EXPECT_FALSE(mtx_is_cwise_equal(A, B));
}
