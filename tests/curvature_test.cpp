#include <cmath>
#include <gtest/gtest.h>
#include <memory>
#include <stdio.h>
#include <string>
#include <vector>

#include "matrix_io.h"
#include <losc/curvature.h>
#include <losc/losc.h>

using losc::Matrix;
using std::shared_ptr;
using std::string;
using std::vector;

using shared_ptr<Matrix> = std::shared_ptr<Matrix>;

struct CurvatureTest : public ::testing::TestWithParam<string> {
    string file_path;
    string dir_path;
    virtual void SetUp() override
    {
        file_path = __FILE__;
        dir_path = file_path.substr(0, file_path.rfind("/"));
    }

    virtual void TearDown() override {}
};

TEST_P(CurvatureTest, non_block_Kappa_J)
{
    // load data.
    string mol = GetParam();
    const char *mol_str = mol.c_str();
    string lo_ref_path = dir_path + "/./data/" + mol + "/lo.txt";
    string df_pmn_path = dir_path + "/./data/" + mol + "/df_pmn.txt";
    string df_Vpq_inv_path =
        dir_path + "/./data/" + mol + "/df_Vpq_inverse.txt";
    string grid_basis_value_path =
        dir_path + "/./data/" + mol + "/grid_basis.txt";
    string grid_weight_path = dir_path + "/./data/" + mol + "/grid_weight.txt";
    string kappa_ref_path = dir_path + "/./data/" + mol + "/kappa.txt";
    auto C_lo = test::read_matrices_from_txt(lo_ref_path);
    auto df_pmn = test::read_matrices_from_txt(df_pmn_path);
    auto df_Vpq_inverse = test::read_matrices_from_txt(df_Vpq_inv_path);
    auto grid_basis_value = test::read_matrices_from_txt(grid_basis_value_path);
    auto kappa_ref = test::read_matrices_from_txt(kappa_ref_path);
    auto grid_wt = test::read_matrices_from_txt(grid_weight_path);
    const size_t npts = grid_wt[0]->size();
    auto grid_weight = std::make_shared<vector<double>>(npts);
    for (size_t i = 0; i < npts; ++i) {
        (*grid_weight)[i] = grid_wt[0]->data()[i];
    }

    for (int is = 0; is < 2; is++) {
        (*C_lo[is]).transposeInPlace();
    }

    for (int is = 0; is < 2; is++) {
        // calculate df_pii
        const size_t nfitbasis = df_Vpq_inverse[0]->rows();
        const size_t nlo = C_lo[is]->cols();
        vector<size_t> p_index(nfitbasis);
        for (size_t i = 0; i < nfitbasis; ++i) {
            p_index[i] = i;
        }
        auto df_pii = std::make_shared<Matrix>(nfitbasis, nlo);
        losc::utils::convert_df_pmn2pii_blockwise(p_index, *df_pmn[0],
                                                  *C_lo[is], *df_pii);

        // Do calculation
        losc::CurvatureV2 kappa_man(losc::B3LYP, C_lo[is], df_pii,
                                    df_Vpq_inverse[0], grid_basis_value[0],
                                    grid_weight);
        auto kappa_calc = kappa_man.compute();

        // Test.
        // True consition is that the calculated curvature matrix matches the
        // reference to the 8-th digit.
        bool status = kappa_calc->is_cwise_equal(*kappa_ref[is], 1e-8);
        if (!status) {
            std::cout << "Mol: " << mol << std::endl;
            std::cout << "data dir path: " << dir_path << std::endl;
            printf("%s, Reference Kappa: spin=%d\n", mol_str, is);
            kappa_ref[is]->show_full();
            printf("%s, Calculated Kappa: spin=%d\n", mol_str, is);
            kappa_calc->show_full();
        }
        EXPECT_TRUE(status);
    }
}

TEST_P(CurvatureTest, block_Kappa_J)
{
    // load data.
    string mol = GetParam();
    const char *mol_str = mol.c_str();
    string lo_ref_path = dir_path + "/./data/" + mol + "/lo.txt";
    string df_pmn_path = dir_path + "/./data/" + mol + "/df_pmn.txt";
    string df_Vpq_inv_path =
        dir_path + "/./data/" + mol + "/df_Vpq_inverse.txt";
    string grid_basis_value_path =
        dir_path + "/./data/" + mol + "/grid_basis.txt";
    string grid_weight_path = dir_path + "/./data/" + mol + "/grid_weight.txt";
    string kappa_ref_path = dir_path + "/./data/" + mol + "/kappa.txt";
    auto C_lo = test::read_matrices_from_txt(lo_ref_path);
    auto df_pmn = test::read_matrices_from_txt(df_pmn_path);
    auto df_Vpq_inverse = test::read_matrices_from_txt(df_Vpq_inv_path);
    auto grid_basis_value = test::read_matrices_from_txt(grid_basis_value_path);
    auto kappa_ref = test::read_matrices_from_txt(kappa_ref_path);
    auto grid_wt = test::read_matrices_from_txt(grid_weight_path);
    const size_t npts = grid_wt[0]->size();
    auto grid_weight = std::make_shared<vector<double>>(npts);
    for (size_t i = 0; i < npts; ++i) {
        (*grid_weight)[i] = grid_wt[0]->data()[i];
    }

    for (int is = 0; is < 2; is++) {
        (*C_lo[is]).transposeInPlace();
    }

    for (int is = 0; is < 2; is++) {
        // calculate df_pii
        const size_t nfitbasis = df_Vpq_inverse[0]->rows();
        const size_t nbasis = C_lo[is]->rows();
        const size_t nlo = C_lo[is]->cols();

        auto df_pii = std::make_shared<Matrix>(nfitbasis, nlo);
        const size_t block_size = 3;
        size_t nblock = nfitbasis / block_size;
        size_t block_tail_size = nfitbasis % block_size;
        if (block_tail_size != 0)
            nblock++;
        for (size_t block_i = 0; block_i < nblock; ++block_i) {
            size_t size = block_size;
            if (block_i == nblock - 1) {
                size = block_tail_size;
            }
            vector<size_t> index(size);
            for (size_t i = 0; i < size; ++i) {
                index[i] = block_i * block_size + i;
            }
            Matrix pmn_block =
                df_pmn[0]->block(index[0], 0, size, nbasis * (nbasis + 1) / 2);

            losc::utils::convert_df_pmn2pii_blockwise(index, pmn_block,
                                                      *C_lo[is], *df_pii);
        }

        // Do calculation
        losc::CurvatureV2 kappa_man(losc::B3LYP, C_lo[is], df_pii,
                                    df_Vpq_inverse[0], grid_basis_value[0],
                                    grid_weight);
        auto kappa_calc = kappa_man.compute();

        // Test.
        // True consition is that the calculated curvature matrix matches the
        // reference to the 8-th digit.
        bool status = kappa_calc->is_cwise_equal(*kappa_ref[is], 1e-8);
        if (!status) {
            std::cout << "Mol: " << mol << std::endl;
            std::cout << "data dir path: " << dir_path << std::endl;
            printf("%s, Reference Kappa: spin=%d\n", mol_str, is);
            kappa_ref[is]->show_full();
            printf("%s, Calculated Kappa: spin=%d\n", mol_str, is);
            kappa_calc->show_full();
        }
        EXPECT_TRUE(status);
    }
}

INSTANTIATE_TEST_SUITE_P(curvature_test, CurvatureTest,
                         ::testing::Values("H2_1.0", "H2_10", "H2+_10", "H2O"));
