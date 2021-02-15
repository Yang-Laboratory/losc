#include <cmath>
#include <fstream>
#include <gtest/gtest.h>
#include <stdio.h>
#include <string>
#include <utility>
#include <vector>

#include "matrix_helper.hpp"
#include "matrix_io.hpp"
#include <losc/local_occupation.hpp>

using std::move;
using std::string;
using std::vector;
using namespace losc;

struct LocalOccupationTest : public ::testing::TestWithParam<string> {
    string file_path;
    string dir_path;

    virtual void SetUp() override
    {
        file_path = __FILE__;
        dir_path = file_path.substr(0, file_path.rfind("/"));
    }

    virtual void TearDown() override {}
};

TEST_P(LocalOccupationTest, test)
{
    // load data.
    string mol = GetParam();
    const char *mol_str = mol.c_str();
    string S_path = dir_path + "/data/" + mol + "/ao_overlap.txt";
    string C_lo_path = dir_path + "/data/" + mol + "/lo.txt";
    string D_path = dir_path + "/data/" + mol + "/dfa_density.txt";
    string L_ref_path = dir_path + "/data/" + mol + "/localocc.txt";
    LOSCMatrix S = move(test::read_matrices_from_txt(S_path)[0]);
    LOSCMatrix C_lo = test::read_matrices_from_txt(C_lo_path)[0].transpose();
    LOSCMatrix D = move(test::read_matrices_from_txt(D_path)[0]);
    LOSCMatrix L_ref = move(test::read_matrices_from_txt(L_ref_path)[0]);

    // Do calculation.
    LOSCMatrix L_calc = losc::local_occupation(C_lo, S, D);

    // Testing.
    // True condition is the calculated Losc local occupation matrix matches
    // the reference to the 8-th digit.
    bool status = test::mtx_is_cwise_equal(L_calc, L_ref, 1e-8);
    EXPECT_TRUE(status);
    if (!status) {
        std::ofstream log;
        string log_file = "LocalOccTest.test." + mol + ".log";
        log.open(log_file);
        log << "data dir path: " << dir_path << std::endl;
        log << "LocalOcc matrix: reference\n";
        test::mtx_show_full(L_ref, log);
        log << "LocalOcc matrix: calculated\n";
        test::mtx_show_full(L_calc, log);
        log.close();
        std::cout << "Test Failure: See log file for more information: "
                  << log_file << std::endl;
    }
}

INSTANTIATE_TEST_SUITE_P(correction_test, LocalOccupationTest,
                         ::testing::Values("H2.1A", "H2.10A", "H2O"));
