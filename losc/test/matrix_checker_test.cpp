#include <gtest/gtest.h>
#include <vector>

#include <losc/losc.h>

using std::vector;
using losc::Matrix;

/**
 * Test is_symmetric() method.
 */
TEST(MatrixCheckerTest, is_symmetric_test)
{
    Matrix A = Matrix::Zero(10, 10);
    // zero matrix should be regarded as symmetric.
    EXPECT_TRUE(A.is_symmetric());
    EXPECT_TRUE(A.is_symmetric(1e-16));
    EXPECT_TRUE(A.is_symmetric(1e-300));

    A.randomize(0, 1);
    EXPECT_FALSE(A.is_symmetric()); // random matrix usually is not symmetric.
    EXPECT_TRUE(A.is_symmetric(2)); // lower threshold to be symmetric.
    A.to_symmetric("L");
    EXPECT_TRUE(A.is_symmetric());  // real symmetric matrix.
    EXPECT_TRUE(A.is_symmetric(-1e-16)); // happended to set it to negative threshold.

    // if customized threshold working?
    double thred = 1e-8;
    A(2,3) = A(3, 2) - thred;
    EXPECT_FALSE(A.is_symmetric(thred/10));
    EXPECT_TRUE(A.is_symmetric(thred * 10));

    // non-square matrix.
    Matrix B(2, 3);
    EXPECT_FALSE(B.is_symmetric());
    EXPECT_FALSE(B.is_symmetric(100));
    EXPECT_FALSE(B.is_symmetric(1e-300));
}

/**
 * Test is_cwise_equal() method.
 */
TEST(MatrixCheckerTest, is_cwise_equal_test)
{
    Matrix A = Matrix::Zero(10, 8);
    Matrix B = Matrix::Zero(10, 8);
    // real zero matrix should be equal.
    EXPECT_TRUE(A.is_cwise_equal(B));
    EXPECT_TRUE(A.is_cwise_equal(B, 1e-16));
    EXPECT_TRUE(A.is_cwise_equal(B, 1e-300));

    // customized threshold
    double thred = 1e-8;
    A(1,1) += thred;
    EXPECT_FALSE(A.is_cwise_equal(B, thred/10));
    EXPECT_TRUE(A.is_cwise_equal(B, thred * 10));

    // random matrix usually is not equal to each other.
    A.randomize(0, 1);
    B.randomize(0, 1);
    EXPECT_FALSE(A.is_cwise_equal(B));
}
