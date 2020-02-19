#include <cmath>
#include <fstream>
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
    string S_path = dir_path + "/./data/" + mol + "/ao_overlap.txt";
    string C_lo_path = dir_path + "/./data/" + mol + "/lo.txt";
    string K_path = dir_path + "/./data/" + mol + "/kappa.txt";
    string L_path = dir_path + "/./data/" + mol + "/localocc.txt";
    string H_losc_path = dir_path + "/./data/" + mol + "/losc_H_corr.txt";
    auto S = matrix::read_matrices_from_txt(S_path);
    auto C_lo = matrix::read_matrices_from_txt(C_lo_path);
    auto K = matrix::read_matrices_from_txt(K_path);
    auto L = matrix::read_matrices_from_txt(L_path);
    auto H_losc_ref = matrix::read_matrices_from_txt(H_losc_path);

    for (int is = 0; is < 2; is++) {
        // Do calculation.
        auto H_losc_calc =
            losc::losc_hamiltonian_correction(*S[0], *C_lo[is], *K[is], *L[is]);

        // Testing.
        // True condition is the calculated Losc H matrix matches the reference
        // to the 8-th digit.
        bool status = H_losc_calc->is_equal_to(*H_losc_ref[is], 1e-8);
        if (!status) {
            std::cout << "Mol: " << mol << std::endl;
            std::cout << "data dir path: " << dir_path << std::endl;
            printf("%s, Reference H_losc: spin=%d\n", mol_str, is);
            H_losc_ref[is]->show_full();
            printf("%s, Calculated H_losc: spin=%d\n", mol_str, is);
            H_losc_calc->show_full();
        }
        EXPECT_TRUE(status);
    }
}

TEST_P(CorrectionTest, losc_energy_correction)
{
    // load data.
    string mol = GetParam();
    const char *mol_str = mol.c_str();
    string K_path = dir_path + "/./data/" + mol + "/kappa.txt";
    string L_path = dir_path + "/./data/" + mol + "/localocc.txt";
    string E_path = dir_path + "/./data/" + mol + "/energy.txt";
    auto K = matrix::read_matrices_from_txt(K_path);
    auto L = matrix::read_matrices_from_txt(L_path);
    auto E_ref = get_losc_E_correction(E_path);

    // Do calculation.
    double E_calc = 0.0;
    for (int is = 0; is < 2; is++) {
        E_calc += losc::losc_total_energy_correction(*K[is], *L[is]);
    }

    // Testing.
    // True condition is the calculated total energy correction matches the
    // reference to the 8-th digit.
    bool status = std::fabs(E_calc - E_ref) < 1e-8;
    EXPECT_TRUE(status);
    if (!status) {
        std::cout << "Calculated total energy correction: " << E_calc
                  << std::endl;
        std::cout << "Referenced total energy correction: " << E_ref
                  << std::endl;
    }
}

TEST_P(CorrectionTest, losc_orbE_correction)
{
    // load data.
    string mol = GetParam();
    const char *mol_str = mol.c_str();
    string S_path = dir_path + "/./data/" + mol + "/ao_overlap.txt";
    string C_co_path = dir_path + "/./data/" + mol + "/lo_basis.txt";
    string C_lo_path = dir_path + "/./data/" + mol + "/lo.txt";
    string K_path = dir_path + "/./data/" + mol + "/kappa.txt";
    string L_path = dir_path + "/./data/" + mol + "/localocc.txt";
    string eig_path = dir_path + "/./data/" + mol + "/losc_eig_direct.txt";
    auto S = matrix::read_matrices_from_txt(S_path);
    auto C_co = matrix::read_matrices_from_txt(C_co_path);
    auto C_lo = matrix::read_matrices_from_txt(C_lo_path);
    auto K = matrix::read_matrices_from_txt(K_path);
    auto L = matrix::read_matrices_from_txt(L_path);
    auto eig_ref = matrix::read_matrices_from_txt(eig_path);

    for (int is = 0; is < 2; is++) {
        // Do calculation.
        auto eig_calc = losc::losc_orbital_energy_correction(
            *S[0], *C_co[is], *C_lo[is], *K[is], *L[is]);

        // Testing.
        // Returned is vector<double>, transfer it into matrix to compare.
        // True condition is that the calculated eigenvalue matches to 8-th
        // digit.
        Matrix eig_calc_M(1, eig_calc.size(), eig_calc);
        bool status = eig_calc_M.is_equal_to(*eig_ref[is], 1e-8);
        EXPECT_TRUE(status);
        if (!status) {
            std::cout << "Calculated eig diag:\n";
            eig_calc_M.show_full();
            std::cout << "Reference eig diag:\n";
            (*eig_ref[is]).show_full();
        }
    }
}

TEST_P(CorrectionTest, losc_orbE_projection)
{
    // load data.
    string mol = GetParam();
    const char *mol_str = mol.c_str();
    string H_dfa_path = dir_path + "/./data/" + mol + "/dfa_h.txt";
    string H_losc_path = dir_path + "/./data/" + mol + "/losc_H_corr.txt";
    string C_co_path = dir_path + "/./data/" + mol + "/lo_basis.txt";
    string eig_path = dir_path + "/./data/" + mol + "/losc_eig_proj.txt";
    auto H_dfa = matrix::read_matrices_from_txt(H_dfa_path);
    auto H_losc = matrix::read_matrices_from_txt(H_losc_path);
    auto C_co = matrix::read_matrices_from_txt(C_co_path);
    auto eig_ref = matrix::read_matrices_from_txt(eig_path);

    for (int is = 0; is < 2; is++) {
        // Do calculation.
        auto eig_calc = losc::losc_corrected_orbital_energy_by_projection(
            *H_dfa[is], *H_losc[is], *C_co[is]);

        // TESTING
        // Return is vector<double>, transfer it into matrix to compare.
        // True condition is that the calculated eigenvalue matches to 8-th
        // digit.
        Matrix eig_calc_M(1, eig_calc.size(), eig_calc);
        bool status = eig_calc_M.is_equal_to(*eig_ref[is], 1e-8);
        EXPECT_TRUE(status);
        if (!status) {
            std::cout << "Calculated eig diag:\n";
            eig_calc_M.show_full();
            std::cout << "Reference eig diag:\n";
            (*eig_ref[is]).show_full();
        }
    }
}

TEST_P(CorrectionTest, losc_orbE_diagonalization)
{
    // load data.
    string mol = GetParam();
    const char *mol_str = mol.c_str();
    string H_dfa_path = dir_path + "/./data/" + mol + "/dfa_h.txt";
    string H_losc_path = dir_path + "/./data/" + mol + "/losc_H_corr.txt";
    string S_path = dir_path + "/./data/" + mol + "/ao_overlap.txt";
    string eig_path = dir_path + "/./data/" + mol + "/losc_eig_diag.txt";
    auto H_dfa = matrix::read_matrices_from_txt(H_dfa_path);
    auto H_losc = matrix::read_matrices_from_txt(H_losc_path);
    auto S = matrix::read_matrices_from_txt(S_path);
    auto eig_ref = matrix::read_matrices_from_txt(eig_path);
    // calculate S^(-1/2);
    vector<double> S_eig(S[0]->row());
    Matrix Shalf(S_eig.size(), S_eig.size());
    Matrix S_eigV(S_eig.size(), S_eig.size());
    matrix::diagonalize_sym_matrix_dsyev("L", *S[0], S_eig);
    for (size_t i = 0; i < S_eig.size(); ++i) {
        S_eigV(i, i) = 1.0 / std::sqrt(S_eig[i]);
    }
    matrix::mult_dgemm_ATBA(*S[0], S_eigV, Shalf);

    for (int is = 0; is < 2; is++) {
        // Do calculation.
        auto eig_calc = losc::losc_corrected_orbital_energy_by_diagonalize(
            *H_dfa[is], *H_losc[is], Shalf);

        // TESTING:
        // Return is vector<double>, transfer it into matrix to compare.
        // True condition is that the calculated eigenvalue matches to 8-th
        // digit.
        Matrix eig_calc_M(1, eig_calc.size(), eig_calc);
        bool status = eig_calc_M.is_equal_to(*eig_ref[is], 1e-8);
        EXPECT_TRUE(status);
        if (!status) {
            std::cout << "Calculated eig diag:\n";
            eig_calc_M.show_full();
            std::cout << "Reference eig diag:\n";
            (*eig_ref[is]).show_full();
        }
    }
}

INSTANTIATE_TEST_SUITE_P(correction_test, CorrectionTest,
                         ::testing::Values("H2_1.0", "H2_10", "H2+_10", "H2O"));
