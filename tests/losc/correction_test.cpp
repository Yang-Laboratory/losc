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
#include <losc/correction.hpp>

using Eigen::MatrixXd;
using std::move;
using std::string;
using std::vector;

/**
 * Get Losc total energy correction from a txt file.
 * The txt file only contains one line.
 */
static double get_losc_E_correction(string &fname)
{
    std::fstream fin;
    fin.open(fname);
    string energy_str;
    std::getline(fin, energy_str);
    fin.close();
    return std::stod(energy_str);
}

struct CorrectionTest : public ::testing::TestWithParam<string> {
    string file_path;
    string dir_path;

    virtual void SetUp() override
    {
        file_path = __FILE__;
        dir_path = file_path.substr(0, file_path.rfind("/"));
    }

    virtual void TearDown() override {}
};

TEST_P(CorrectionTest, losc_Hamiltonian_correction)
{
    // load data.
    string mol = GetParam();
    const char *mol_str = mol.c_str();
    string S_path = dir_path + "/data/" + mol + "/ao_overlap.txt";
    string C_lo_path = dir_path + "/data/" + mol + "/lo.txt";
    string K_path = dir_path + "/data/" + mol + "/kappa.txt";
    string L_path = dir_path + "/data/" + mol + "/localocc.txt";
    string H_losc_path = dir_path + "/data/" + mol + "/losc_H_corr.txt";
    MatrixXd S = move(test::read_matrices_from_txt(S_path)[0]);
    MatrixXd C_lo =
        move(test::read_matrices_from_txt(C_lo_path)[0].transpose());
    MatrixXd K = move(test::read_matrices_from_txt(K_path)[0]);
    MatrixXd L = move(test::read_matrices_from_txt(L_path)[0]);
    MatrixXd H_losc_ref = move(test::read_matrices_from_txt(H_losc_path)[0]);

    // Do calculation.
    MatrixXd H_losc_calc = losc::ao_hamiltonian_correction(S, C_lo, K, L);

    // Testing.
    // True condition is the calculated Losc H matrix matches the reference
    // to the 8-th digit.
    bool status = test::mtx_is_cwise_equal(H_losc_calc, H_losc_ref, 1e-8);
    EXPECT_TRUE(status);
    if (!status) {
        std::ofstream log;
        string log_file = "CorrectionTest.losc_H." + mol + ".log";
        log.open(log_file);
        log << "data dir path: " << dir_path << std::endl;
        log << "H_losc: reference" << std::endl;
        test::mtx_show_full(H_losc_ref, log);
        log << "H_losc: calculated" << std::endl;
        test::mtx_show_full(H_losc_calc, log);
        log.close();
        std::cout << "Test Failure: See log file for more information: "
                  << log_file << std::endl;
    }
}

TEST_P(CorrectionTest, losc_energy_correction)
{
    // load data.
    string mol = GetParam();
    const char *mol_str = mol.c_str();
    string K_path = dir_path + "/data/" + mol + "/kappa.txt";
    string L_path = dir_path + "/data/" + mol + "/localocc.txt";
    string E_path = dir_path + "/data/" + mol + "/energy.txt";
    MatrixXd K = move(test::read_matrices_from_txt(K_path)[0]);
    MatrixXd L = move(test::read_matrices_from_txt(L_path)[0]);
    double E_ref = get_losc_E_correction(E_path);

    // assume it is all RKS calculations.
    double E_calc = 2 * losc::energy_correction(K, L);

    // Testing.
    // True condition is the calculated total energy correction matches the
    // reference to the 8-th digit.
    bool status = std::abs(E_calc - E_ref) < 1e-8;
    EXPECT_TRUE(status);
    if (!status) {
        std::ofstream log;
        string log_file = "CorrectionTest.losc_energy." + mol + ".log";
        log.open(log_file);
        log << "Calculated total energy correction: " << E_calc << std::endl;
        log << "Referenced total energy correction: " << E_ref << std::endl;
        log.close();
        std::cout << "Test Failure: See log file for more information: "
                  << log_file << std::endl;
    }
}

TEST_P(CorrectionTest, losc_orbE_projection)
{
    // load data.
    string mol = GetParam();
    const char *mol_str = mol.c_str();
    string H_dfa_path = dir_path + "/data/" + mol + "/dfa_h.txt";
    string H_losc_path = dir_path + "/data/" + mol + "/losc_H_corr.txt";
    string C_co_path = dir_path + "/data/" + mol + "/lo_basis.txt";
    string eig_path = dir_path + "/data/" + mol + "/losc_eig.txt";
    MatrixXd H_dfa = move(test::read_matrices_from_txt(H_dfa_path)[0]);
    MatrixXd H_losc = move(test::read_matrices_from_txt(H_losc_path)[0]);
    MatrixXd C_co = test::read_matrices_from_txt(C_co_path)[0].transpose();
    MatrixXd eig_ref = move(test::read_matrices_from_txt(eig_path)[0]);

    // Do calculation.
    vector<double> eig_calc =
        losc::orbital_energy_post_scf(H_dfa, H_losc, C_co);

    // TESTING
    // Return is vector<double>, transfer it into matrix to compare.
    // True condition is that the calculated eigenvalue matches to 8-th
    // digit.
    MatrixXd eig_calc_M(1, eig_calc.size());
    eig_calc_M.setZero();
    for (size_t i = 0; i < eig_calc.size(); ++i) {
        eig_calc_M(0, i) = eig_calc[i];
    }
    bool status = test::mtx_is_cwise_equal(eig_calc_M, eig_ref, 1e-8);
    EXPECT_TRUE(status);
    if (!status) {
        std::ofstream log;
        string log_file = "CorrectionTest.losc_orbital_energy." + mol + ".log";
        log.open(log_file);
        log << "Calculated eig diag:\n";
        test::mtx_show_full(eig_calc_M, log);
        log << "Reference eig diag:\n";
        test::mtx_show_full(eig_ref, log);
        std::cout << "Test Failure: See log file for more information: "
                  << log_file << std::endl;
    }
}

INSTANTIATE_TEST_SUITE_P(correction_test, CorrectionTest,
                         ::testing::Values("H2.1A", "H2.10A", "H2O"));
