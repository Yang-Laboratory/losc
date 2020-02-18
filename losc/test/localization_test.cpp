#include <matrix/matrix.h>
#include <gtest/gtest.h>
#include <vector>
#include <memory>
#include <stdio.h>
#include <string>
#include <cmath>

#include <losc/localization.h>

using std::vector;
using std::shared_ptr;
using matrix::Matrix;
using std::string;

using SharedMatrix = std::shared_ptr<Matrix>;


bool is_permutation_matrix(Matrix &A, double thredshold=1e-10)
{
    if (!A.is_square()) {
        return false;
    }

    for (size_t i = 0; i < A.row(); ++i) {
        size_t n_one = 0;
        for (size_t j = 0; j < A.col(); ++j) {
            if (!(A(i, j) < thredshold) &&
                !(std::fabs((A(i, j) - 1.0)) < thredshold)
               ) {
                   return false;
            } else if (std::fabs(A(i, j) - 1.0) < thredshold) {
                n_one++;
                if (n_one > 1) {
                    return false;
                }
            }
        }
    }

    for (size_t j = 0; j < A.col(); ++j) {
        size_t n_one = 0;
        for (size_t i = 0; i < A.row(); ++i) {
            if (!(A(i, j) < thredshold) &&
                !(std::fabs((A(i, j) - 1.0)) < thredshold)
               ) {
                   return false;
            } else if (std::fabs(A(i, j) - 1.0) < thredshold) {
                n_one++;
                if (n_one > 1) {
                    return false;
                }
            }
        }
    }
    return true;
}

struct LocalizationTest:public ::testing::TestWithParam<string> {
    size_t nlo = 4;
    size_t nbasis = 4;

    vector<shared_ptr<Matrix>> C_lo_basis;
    vector<shared_ptr<Matrix>> H_ao;
    vector<shared_ptr<Matrix>> D_ao;
    vector<shared_ptr<Matrix>> C_lo_ref;
    vector<shared_ptr<Matrix>> S_ao;

    string file_path;
    string dir_path;

    virtual void SetUp() override
    {
        file_path = __FILE__;
        dir_path = file_path.substr(0, file_path.rfind("/"));
    }

    virtual void TearDown() override
    {

    }
};

TEST_P(LocalizationTest, H2)
{
    string mol = GetParam();
    const char *mol_str = mol.c_str();
    string lo_basis_path = dir_path + "/./data/" + mol + "/lo_basis.txt";
    string lo_ref_path = dir_path + "/./data/" + mol + "/lo.txt";
    string H_ao_path = dir_path + "/./data/" + mol + "/dfa_h.txt";
    string D_ao_path = dir_path + "/./data/" + mol + "/dipole_ao.txt";
    string S_ao_path = dir_path + "/./data/" + mol + "/ao_overlap.txt";

    std::cout << "Mol: " << mol << std::endl;
    std::cout << "data dir path: " << dir_path << std::endl;
    std::cout << "LO basis coefficient: " << lo_basis_path << std::endl;
    std::cout << "LO coefficient reference : " << lo_ref_path << std::endl;
    std::cout << "DFA H under AO: " << H_ao_path << std::endl;
    std::cout << "Dipole under AO: " << D_ao_path << std::endl;
    std::cout << "AO overlap: " << S_ao_path << std::endl;

    auto C_lo_basis = matrix::read_matrices_from_txt(lo_basis_path);
    auto C_lo_ref = matrix::read_matrices_from_txt(lo_ref_path);
    auto H_ao = matrix::read_matrices_from_txt(H_ao_path);
    auto D_ao = matrix::read_matrices_from_txt(D_ao_path);
    auto S_ao = matrix::read_matrices_from_txt(S_ao_path);

    for (int is = 0; is < 2; is++) {
        // do localization.
        losc::LoscLocalizerV2 localizer(C_lo_basis[is], H_ao[is], D_ao);
        localizer.set_random_permutation(false);
        //localizer.set_print(losc::kPrintLevelNormal);
        auto C_lo_calc = localizer.compute();

        // verify the correctness of calculated LO coefficient matrix.
        bool is_equal = C_lo_calc->is_equal_to(*C_lo_ref[is], 1e-8);
        printf("%s, spin=%d: equal to the reference: %s\n", mol_str, is, is_equal ? "yes" : "no");
        SharedMatrix P;
        bool is_permutation = false;
        if (! is_equal) {
            SharedMatrix C1 = C_lo_calc;
            SharedMatrix C2 = C_lo_ref[is];
            SharedMatrix S = S_ao[0];
            SharedMatrix C1SC1 = std::make_shared<Matrix> (nlo, nlo);
            matrix::mult_dgemm_ABAT(*C1, *S, *C1SC1);
            EXPECT_TRUE(C1SC1->is_identity(1e-8));

            P = std::make_shared<Matrix> (nlo, nlo);
            SharedMatrix C2S = std::make_shared<Matrix> (nlo, nbasis);
            matrix::mult_dgemm(1.0, *C2, "N", *S, "N", 0.0, *C2S);
            matrix::mult_dgemm(1.0, *C2S, "N", *C1, "T", 0.0, *P);

            // P = abs(P)
            for (size_t i = 0; i < P->size(); ++i) {
                P->data()[i] = std::fabs(P->data()[i]);
            }
            is_permutation = is_permutation_matrix(*P, 1e-8);
            printf("%s, spin=%d: relate to the reference with a permutation matrix: %s\n",
                   mol_str, is, is_permutation ? "yes" : "no");
        }

        // Test failed. Print out matrices to check.
        if (!(is_equal || is_permutation)) {
            printf("%s: AO overlap matrix:\n", mol_str);
            S_ao[0]->show_full();
            for (int xyz = 0; xyz < 3; xyz++) {
                printf("%s: Dipole matrix under AO: xyz=[%d]\n", mol_str, xyz);
                D_ao[xyz]->show_full();
            }
            printf("%s: LO basis coefficient matrix: spin=%d\n", mol_str, is);
            C_lo_basis[is]->show_full();

            printf("%s: Hamiltonian under AO: spin=%d\n", mol_str, is);
            H_ao[is]->show_full();

            printf("%s: LO coefficient matrix: reference. spin=%d\n", mol_str, is);
            C_lo_ref[is]->show_full();

            printf("%s: LO coefficient matrix: calculated. spin=%d\n", mol_str, is);
            C_lo_calc->show_full();

            if (!is_equal) {
                printf("%s: calculated LO coefficient matrix may relates to reference with a permutation matrix:\n", mol_str);
                printf("The permutation matrix is:\n");
                P->show_full();
            }
        }

        EXPECT_TRUE(is_equal || is_permutation);
    }
}

INSTANTIATE_TEST_SUITE_P(
    Localizationtest,
    LocalizationTest,
    ::testing::Values(
            "H2_1.0", "H2_10", "H2+_10",
            "H2O")
);
