#include <cmath>
#include <gtest/gtest.h>
#include <matrix/matrix.h>
#include <memory>
#include <stdio.h>
#include <string>
#include <vector>

#include <losc/correction.h>

using matrix::Matrix;
using std::shared_ptr;
using std::string;
using std::vector;

using SharedMatrix = std::shared_ptr<Matrix>;

struct CorrectionTest : public ::testing::TestWithParam<string> {
    vector<shared_ptr<Matrix>> S;
    vector<shared_ptr<Matrix>> C_lo;
    vector<shared_ptr<Matrix>> C_co;
    vector<shared_ptr<Matrix>> K;
    vector<shared_ptr<Matrix>> L;
    vector<shared_ptr<Matrix>> H_dfa;
    vector<shared_ptr<Matrix>> H_losc_ref;

    string file_path;
    string dir_path;

    virtual void SetUp() override
    {
        file_path = __FILE__;
        dir_path = file_path.substr(0, file_path.rfind("/"));
    }

    virtual void TearDown() override {}
};

TEST_P(CorrectionTest, losc_H)
{
    // load data.
    string mol = GetParam();
    const char *mol_str = mol.c_str();
    string S_path = dir_path + "/./data/" + mol + "/ao_overlap.txt";
    string C_lo_path = dir_path + "/./data/" + mol + "/lo.txt";
    string C_co_path = dir_path + "/./data/" + mol + "/lo_basis.txt";
    string K_path = dir_path + "/./data/" + mol + "/kappa.txt";
    string L_path = dir_path + "/./data/" + mol + "/localocc.txt";
    string H_dfa_path = dir_path + "/./data/" + mol + "/dfa_h.txt";
    string H_losc_path = dir_path + "/./data/" + mol + "/losc_H_corr.txt";

    S = matrix::read_matrices_from_txt(S_path);
    C_lo = matrix::read_matrices_from_txt(C_lo_path);
    C_co = matrix::read_matrices_from_txt(C_co_path);
    K = matrix::read_matrices_from_txt(K_path);
    L = matrix::read_matrices_from_txt(L_path);
    H_dfa = matrix::read_matrices_from_txt(H_dfa_path);
    H_losc_ref = matrix::read_matrices_from_txt(H_losc_path);

    std::cout << "Mol: " << mol << std::endl;
    std::cout << "data dir path: " << dir_path << std::endl;

    for (int is = 0; is < 2; is++) {
        // compute curvature matrix.
        auto H_losc_calc =
            losc::losc_hamiltonian_correction(*S[0], *C_lo[is], *K[is], *L[is]);

        // verify the results.
        bool status = H_losc_calc->is_equal_to(*H_losc_ref[is], 1e-8);
        if (!status) {
            printf("%s, Reference H_losc: spin=%d\n", mol_str, is);
            H_losc_ref[is]->show_full();
            printf("%s, Calculated H_losc: spin=%d\n", mol_str, is);
            H_losc_calc->show_full();
        }
        EXPECT_TRUE(status);
    }
}

INSTANTIATE_TEST_SUITE_P(correction_test, CorrectionTest,
                         ::testing::Values("H2_1.0", "H2_10", "H2+_10", "H2O"));
