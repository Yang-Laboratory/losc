#include <matrix/matrix.h>
#include <gtest/gtest.h>
#include <vector>
#include <memory>
#include <stdio.h>
#include <string>
#include <cmath>

#include <losc/curvature.h>

using std::vector;
using std::shared_ptr;
using matrix::Matrix;
using std::string;

using SharedMatrix = std::shared_ptr<Matrix>;

struct CurvatureTest:public ::testing::TestWithParam<string> {
    vector<shared_ptr<Matrix>> C_lo;
    vector<shared_ptr<Matrix>> df_pmn;
    vector<shared_ptr<Matrix>> df_Vpq_inverse;
    vector<shared_ptr<Matrix>> grid_basis_value;
    vector<shared_ptr<Matrix>> kappa_ref;
    std::shared_ptr<vector<double>> grid_weight;

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

TEST_P(CurvatureTest, test)
{
    // load data.
    string mol = GetParam();
    const char *mol_str = mol.c_str();
    string lo_ref_path = dir_path + "/./data/" + mol + "/lo.bin";
    string df_pmn_path = dir_path + "/./data/" + mol + "/df_pmn.bin";
    string df_Vpq_inv_path = dir_path + "/./data/" + mol + "/df_Vpq_inverse.bin";
    string grid_basis_value_path = dir_path + "/./data/" + mol + "/grid_basis.bin";
    string grid_weight_path = dir_path + "/./data/" + mol + "/grid_weight.bin";
    string kappa_ref_path = dir_path + "/./data/" + mol + "/kappa.bin";

    C_lo = matrix::read_matrices_from_binary(lo_ref_path);
    df_pmn = matrix::read_matrices_from_binary(df_pmn_path);
    df_Vpq_inverse = matrix::read_matrices_from_binary(df_Vpq_inv_path);
    grid_basis_value = matrix::read_matrices_from_binary(grid_basis_value_path);
    kappa_ref = matrix::read_matrices_from_binary(kappa_ref_path);

    auto grid_wt = matrix::read_matrices_from_binary(grid_weight_path);
    const size_t npts = grid_wt[0]->size();
    grid_weight = std::make_shared<vector<double>>(npts);
    for (size_t i = 0; i < npts; ++i) {
        (*grid_weight)[i] = grid_wt[0]->data()[i];
    }

    std::cout << "Mol: " << mol << std::endl;
    std::cout << "data dir path: " << dir_path << std::endl;

    for (int is = 0; is < 2; is++) {
        // compute curvature matrix.
        losc::CurvatureV2 kappa_man(losc::B3LYP, C_lo[is], df_pmn[0], df_Vpq_inverse[0], grid_basis_value[0], grid_weight);
        kappa_man.compute();
        auto kappa_calc = kappa_man.get_curvature();

        // verify the results.
        bool status = kappa_calc->is_equal_to(*kappa_ref[is], 1e-8);
        if (!status) {
            printf("%s, Reference Kappa: spin=%d\n", mol_str, is);
            kappa_ref[is]->show_full();
            printf("%s, Calculated Kappa: spin=%d\n", mol_str, is);
            kappa_calc->show_full();
        }
        EXPECT_TRUE(status);
    }
}

INSTANTIATE_TEST_SUITE_P(
    curvature_test,
    CurvatureTest,
    ::testing::Values(
            "H2_1.0", "H2_10", "H2+_10",
            "H2O")
);
