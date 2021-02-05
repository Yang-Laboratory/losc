#include <cmath>
#include <gtest/gtest.h>
#include <stdio.h>
#include <string>
#include <vector>

#include "matrix_io.hpp"
#include "matrix_helper.hpp"
#include <losc/localization.hpp>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::string;
using std::vector;

bool is_permutation_matrix(MatrixXd &A, double threshold = 1e-10)
{
    if (!test::mtx_is_square(A)) {
        return false;
    }

    for (size_t i = 0; i < A.rows(); ++i) {
        size_t n_one = 0;
        for (size_t j = 0; j < A.cols(); ++j) {
            if (!(A(i, j) < threshold) &&
                !(std::abs((A(i, j) - 1.0)) < threshold)) {
                return false;
            } else if (std::abs(A(i, j) - 1.0) < threshold) {
                n_one++;
                if (n_one > 1) {
                    return false;
                }
            }
        }
    }

    for (size_t j = 0; j < A.cols(); ++j) {
        size_t n_one = 0;
        for (size_t i = 0; i < A.rows(); ++i) {
            if (!(A(i, j) < threshold) &&
                !(std::abs((A(i, j) - 1.0)) < threshold)) {
                return false;
            } else if (std::abs(A(i, j) - 1.0) < threshold) {
                n_one++;
                if (n_one > 1) {
                    return false;
                }
            }
        }
    }
    return true;
}

struct LocalizationTest : public ::testing::TestWithParam<string> {
    string file_path;
    string dir_path;

    virtual void SetUp() override
    {
        file_path = __FILE__;
        dir_path = file_path.substr(0, file_path.rfind("/"));
    }

    virtual void TearDown() override {}
};

TEST_P(LocalizationTest, H2)
{
    // load data
    string mol = GetParam();
    const char *mol_str = mol.c_str();
    string lo_basis_path = dir_path + "/./data/" + mol + "/lo_basis.txt";
    string lo_ref_path = dir_path + "/./data/" + mol + "/lo.txt";
    string H_ao_path = dir_path + "/./data/" + mol + "/dfa_h.txt";
    string D_ao_path = dir_path + "/./data/" + mol + "/dipole_ao.txt";
    string S_ao_path = dir_path + "/./data/" + mol + "/ao_overlap.txt";
    auto C_lo_basis = test::read_matrices_from_txt(lo_basis_path);
    auto C_lo_ref = test::read_matrices_from_txt(lo_ref_path);
    auto H_ao = test::read_matrices_from_txt(H_ao_path);
    auto D_ao_tmp = test::read_matrices_from_txt(D_ao_path);
    auto S_ao = test::read_matrices_from_txt(S_ao_path);
    vector<losc::RefConstMat> D_ao(D_ao_tmp.begin(), D_ao_tmp.end());

    for (int is = 0; is < 2; is++) {
        C_lo_ref[is].transposeInPlace();
        C_lo_basis[is].transposeInPlace();
    }

    for (int is = 0; is < 2; is++) {
        // Do localization.
        losc::LocalizerV2 localizer(C_lo_basis[is], H_ao[is], D_ao);
        localizer.set_random_permutation(false);
        // localizer.set_print(losc::kPrintLevelNormal);
        MatrixXd C_lo_calc = localizer.lo();

        // Test.
        // True condition:
        // 1. Directly compare the calculated LO coefficient matrix with the
        // reference, to see if they can match up to 8-th digit.
        // 2. Try to check if the calculated one is related to reference with
        // a permutation transformation (negative sign is also okay, since it
        // will change the physical observables).
        bool is_equal = test::mtx_is_cwise_equal(C_lo_calc, C_lo_ref[is], 1e-8);
        size_t nlo = C_lo_basis[is].cols();
        MatrixXd P(nlo, nlo);
        bool is_permutation = false;
        if (!is_equal) {
            size_t nlo = C_lo_basis[is].cols();
            size_t nbasis = C_lo_basis[is].rows();
            MatrixXd &C1 = C_lo_calc;
            MatrixXd &C2 = C_lo_ref[is];
            MatrixXd &S = S_ao[0];
            MatrixXd C1SC1 = C1.transpose() * S * C1;
            EXPECT_TRUE(C1SC1.isIdentity(1e-8));

            P.noalias() = C1.transpose() * S * C2;

            // P = abs(P)
            P = P.array().abs().matrix();
            is_permutation = is_permutation_matrix(P, 1e-8);
        }
        EXPECT_TRUE(is_equal || is_permutation);

        // If testing failed, print out matrices to check.
        if (!(is_equal || is_permutation)) {
            printf("%s, spin=%d: equal to the reference: %s\n", mol_str, is,
                   is_equal ? "yes" : "no");
            printf("%s, spin=%d: relate to the reference with a permutation "
                   "matrix: %s\n",
                   mol_str, is, is_permutation ? "yes" : "no");
            std::cout << "Mol: " << mol << std::endl;
            std::cout << "data dir path: " << dir_path << std::endl;
            std::cout << "LO basis coefficient: " << lo_basis_path << std::endl;
            std::cout << "LO coefficient reference : " << lo_ref_path
                      << std::endl;
            std::cout << "DFA H under AO: " << H_ao_path << std::endl;
            std::cout << "Dipole under AO: " << D_ao_path << std::endl;
            std::cout << "AO overlap: " << S_ao_path << std::endl;
            printf("%s: AO overlap matrix:\n", mol_str);
            test::mtx_show_full(S_ao[0]);
            for (int xyz = 0; xyz < 3; xyz++) {
                printf("%s: Dipole matrix under AO: xyz=[%d]\n", mol_str, xyz);
                test::mtx_show_full(D_ao[xyz]);
            }
            printf("%s: LO basis coefficient matrix: spin=%d\n", mol_str, is);
            test::mtx_show_full(C_lo_basis[is]);
            printf("%s: Hamiltonian under AO: spin=%d\n", mol_str, is);
            test::mtx_show_full(H_ao[is]);
            printf("%s: LO coefficient matrix: reference. spin=%d\n", mol_str,
                   is);
            test::mtx_show_full(C_lo_ref[is]);
            printf("%s: LO coefficient matrix: calculated. spin=%d\n", mol_str,
                   is);
            test::mtx_show_full(C_lo_calc);
            if (!is_permutation) {
                printf("%s: Permutation matrix: calculated. spin=%d\n", mol_str,
                       is);
                test::mtx_show_full(P);
            }
        }
    }
}

INSTANTIATE_TEST_SUITE_P(Localizationtest, LocalizationTest,
                         ::testing::Values("H2_1.0", "H2_10", "H2+_10", "H2O"));
