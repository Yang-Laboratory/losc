#include <cmath>
#include <fstream>
#include <gtest/gtest.h>
#include <memory>
#include <stdio.h>
#include <string>
#include <utility>
#include <vector>

#include "matrix_helper.hpp"
#include "matrix_io.hpp"
#include <losc/curvature.hpp>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::string;
using std::vector;

// input matrices data.
typedef struct Data {
    string mol;
    string df_pmn_path;
    string df_Vpq_inv_path;
    string grid_basis_value_path;
    string grid_weight_path;
    string lo_path;
    string kappa_ref_path;

    MatrixXd df_pmn;
    MatrixXd df_Vpq_inv;
    MatrixXd grid_basis_value;
    VectorXd grid_weight;
    MatrixXd grid_lo;
    MatrixXd C_lo;
    MatrixXd kappa_ref;
} Data;

// load input matrices data from txt file.
static Data load_data(string mol)
{
    using std::move;
    string file_path = __FILE__;
    string dir_path = file_path.substr(0, file_path.rfind("/"));
    Data rst;
    rst.lo_path = dir_path + "/data/" + mol + "/lo.txt";
    rst.df_pmn_path = dir_path + "/data/" + mol + "/df_pmn.txt";
    rst.df_Vpq_inv_path = dir_path + "/data/" + mol + "/df_Vpq_inv.txt";
    rst.grid_basis_value_path = dir_path + "/data/" + mol + "/grid_basis.txt";
    rst.grid_weight_path = dir_path + "/data/" + mol + "/grid_wt.txt";
    rst.kappa_ref_path = dir_path + "/data/" + mol + "/kappa.txt";
    rst.C_lo = test::read_matrices_from_txt(rst.lo_path)[0].transpose();
    rst.kappa_ref = move(test::read_matrices_from_txt(rst.kappa_ref_path)[0]);
    rst.df_pmn = move(test::read_matrices_from_txt(rst.df_pmn_path)[0]);
    rst.df_Vpq_inv = move(test::read_matrices_from_txt(rst.df_Vpq_inv_path)[0]);
    rst.grid_basis_value =
        move(test::read_matrices_from_txt(rst.grid_basis_value_path)[0]);
    rst.grid_weight =
        move(test::read_matrices_from_txt(rst.grid_weight_path)[0].transpose());
    rst.grid_lo = rst.grid_basis_value * rst.C_lo;
    return move(rst);
}

struct CurvatureTest : public ::testing::TestWithParam<string> {
    virtual void SetUp() override {}
    virtual void TearDown() override {}
};

TEST_P(CurvatureTest, non_block_kappa)
{
    // load data.
    string mol = GetParam();
    Data data = load_data(mol);

    // calculate df_pii
    const size_t nfitbasis = data.df_Vpq_inv.rows();
    const size_t nlo = data.C_lo.cols();
    vector<size_t> p_index(nfitbasis);
    for (size_t i = 0; i < nfitbasis; ++i) {
        p_index[i] = i;
    }
    MatrixXd df_pii(nfitbasis, nlo);
    losc::utils::convert_df_pmn2pii_blockwise(p_index, data.df_pmn, data.C_lo,
                                              df_pii);

    // Do calculation
    losc::DFAInfo b3lyp(0.8, 0.2, "b3lyp");
    losc::CurvatureV2 kappa_man(b3lyp, df_pii, data.df_Vpq_inv, data.grid_lo,
                                data.grid_weight);
    MatrixXd kappa_calc = kappa_man.kappa();

    // Test.
    // True consition is that the calculated curvature matrix matches the
    // reference to the 8-th digit.
    bool status = test::mtx_is_cwise_equal(kappa_calc, data.kappa_ref, 1e-8);
    if (!status) {
        std::ofstream log;
        log.open("CurvatureTest.no_block_kappa." + mol + ".log");
        log << "Kappa matrix: reference\n";
        test::mtx_show_full(data.kappa_ref, log);
        log << "Kappa matrix: calculated\n";
        test::mtx_show_full(kappa_calc, log);
        log.close();
    }
    EXPECT_TRUE(status);
}

TEST_P(CurvatureTest, block_kappa)
{
    // load data.
    string mol = GetParam();
    Data data = load_data(mol);

    // calculate df_pii
    const size_t nfitbasis = data.df_Vpq_inv.rows();
    const size_t nbasis = data.C_lo.rows();
    const size_t nlo = data.C_lo.cols();
    MatrixXd df_pii(nfitbasis, nlo);
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
        MatrixXd pmn_block =
            data.df_pmn.block(index[0], 0, size, nbasis * (nbasis + 1) / 2);

        losc::utils::convert_df_pmn2pii_blockwise(index, pmn_block, data.C_lo,
                                                  df_pii);
    }

    // Do calculation
    losc::DFAInfo b3lyp(0.8, 0.2, "b3lyp");
    losc::CurvatureV2 kappa_man(b3lyp, df_pii, data.df_Vpq_inv, data.grid_lo,
                                data.grid_weight);
    MatrixXd kappa_calc = kappa_man.kappa();

    // Test.
    // True consition is that the calculated curvature matrix matches the
    // reference to the 8-th digit.
    bool status = test::mtx_is_cwise_equal(kappa_calc, data.kappa_ref, 1e-8);
    EXPECT_TRUE(status);
    if (!status) {
        std::ofstream log;
        string log_file = "CurvatureTest.block_kappa." + mol + ".log";
        log.open(log_file);
        log << "Kappa matrix: reference\n";
        test::mtx_show_full(data.kappa_ref, log);
        log << "Kappa matrix: calculated\n";
        test::mtx_show_full(kappa_calc, log);
        log.close();
        std::cout << "Test Failure: See log file for more information: "
                  << log_file << std::endl;
    }
}

INSTANTIATE_TEST_SUITE_P(curvature_test, CurvatureTest,
                         ::testing::Values("H2.1A", "H2.10A", "H2O"));
