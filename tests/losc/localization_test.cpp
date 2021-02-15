#include <cmath>
#include <fstream>
#include <gtest/gtest.h>
#include <map>
#include <stdio.h>
#include <string>
#include <utility>
#include <vector>

#include "matrix_helper.hpp"
#include "matrix_io.hpp"
#include <losc/localization.hpp>

using std::move;
using std::string;
using std::vector;
using namespace losc;

bool is_unitary(const LOSCMatrix &U, double prec = 1.e-8)
{
    LOSCMatrix UUT = U * U.transpose();
    LOSCMatrix UTU = U.transpose() * U;
    return UUT.isIdentity(prec) && UTU.isIdentity(prec);
}

bool is_permutation_matrix(LOSCMatrix &A, double threshold = 1e-10)
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

struct CaseInfo {
    string dir;   // dir name of the test case.
    string label; // the label of localization setting.
    double gamma;
};

struct LocalizationTest : public ::testing::TestWithParam<CaseInfo> {
    virtual void SetUp() override {}

    virtual void TearDown() override {}
};

TEST_P(LocalizationTest, localizer)
{
    string file_path = __FILE__;
    string dir_path = file_path.substr(0, file_path.rfind("/"));
    CaseInfo test_case = GetParam();
    string dir = test_case.dir;     // source dir of the data
    string label = test_case.label; // label for different localizations.
    string test_name =
        test_case.dir + ":" + test_case.label; // test name of the localization.

    // common data for the parent DFA.
    string H_ao_path = dir_path + "/data/" + dir + "/dfa_h.txt";
    string D_ao_path = dir_path + "/data/" + dir + "/dipole_ao.txt";
    string S_ao_path = dir_path + "/data/" + dir + "/ao_overlap.txt";
    LOSCMatrix H_ao = move(test::read_matrices_from_txt(H_ao_path)[0]);
    LOSCMatrix S_ao = move(test::read_matrices_from_txt(S_ao_path)[0]);
    auto D_ao_tmp = test::read_matrices_from_txt(D_ao_path);
    vector<losc::RefConstMat> D_ao(D_ao_tmp.begin(), D_ao_tmp.end());

    // different LO data
    string lo_basis_file = "lo_basis.txt";
    string lo_ref_file = "lo.txt";
    string U_ref_file = "u.txt";
    if (!label.empty()) {
        lo_basis_file = "lo_basis." + label + ".txt";
        lo_ref_file = "lo." + label + ".txt";
        U_ref_file = "u." + label + ".txt";
    }
    string lo_basis_path = dir_path + "/data/" + dir + "/" + lo_basis_file;
    string lo_ref_path = dir_path + "/data/" + dir + "/" + lo_ref_file;
    string U_ref_path = dir_path + "/data/" + dir + "/" + U_ref_file;
    LOSCMatrix C_lo_basis =
        test::read_matrices_from_txt(lo_basis_path)[0].transpose();
    LOSCMatrix C_lo_ref =
        test::read_matrices_from_txt(lo_ref_path)[0].transpose();
    LOSCMatrix U_ref = test::read_matrices_from_txt(U_ref_path)[0];

    // Do localization.
    losc::LocalizerV2 localizer(C_lo_basis, H_ao, D_ao);
    localizer.set_random_permutation(false);
    localizer.set_gamma(test_case.gamma);
    // localizer.set_print(losc::kPrintLevelNormal);
    vector<LOSCMatrix> rst = localizer.lo_U();
    LOSCMatrix &C_lo_calc = rst[0];
    LOSCMatrix &U_calc = rst[1];

    // Test.
    // True condition:
    // 1. Directly compare the calculated LO coefficient matrix with the
    // reference, to see if they can match up to 8-th digit.
    // 2. Try to check if the calculated one is related to reference with
    // a permutation transformation (negative sign is also okay, since it
    // will change the physical observables).
    bool is_equal = test::mtx_is_cwise_equal(C_lo_calc, C_lo_ref, 1e-8);
    size_t nlo = C_lo_basis.cols();
    LOSCMatrix P(nlo, nlo);
    bool is_permutation = false;
    if (!is_equal) {
        size_t nlo = C_lo_basis.cols();
        size_t nbasis = C_lo_basis.rows();
        LOSCMatrix &C1 = C_lo_calc;
        LOSCMatrix &C2 = C_lo_ref;
        LOSCMatrix C1SC1 = C1.transpose() * S_ao * C1;
        EXPECT_TRUE(C1SC1.isIdentity(1e-8));
        P.noalias() = C1.transpose() * S_ao * C2;

        // P = abs(P)
        P = P.array().abs().matrix();
        is_permutation = is_permutation_matrix(P, 1e-8);
    }
    EXPECT_TRUE(is_equal || is_permutation);

    // If testing failed, write matrices into log file.
    if (!(is_equal || is_permutation)) {
        std::ofstream log;
        string log_file = "LocalizationTest.localizer." + test_name + ".log";
        log.open("LocalizationTest.localizer." + test_name + ".log");
        log << "Test name: " << test_name << std::endl;
        log << "Equal to the reference: " << (is_equal ? "yes" : "no")
            << std::endl;
        log << "Relate to the reference with a permutation matrix: "
            << (is_permutation ? "yes" : "no") << std::endl;
        log << std::endl;
        log << "data dir path: " << dir_path << std::endl;
        log << "LO basis coefficient: " << lo_basis_path << std::endl;
        log << "LO coefficient reference : " << lo_ref_path << std::endl;
        log << "DFA H under AO: " << H_ao_path << std::endl;
        log << "Dipole under AO: " << D_ao_path << std::endl;
        log << "AO overlap: " << S_ao_path << std::endl;
        log << std::endl;
        log << "AO overlap matrix:" << std::endl;
        test::mtx_show_full(S_ao, log);
        for (int xyz = 0; xyz < 3; xyz++) {
            log << "Dipole matrix under AO: "
                << "xyz=" << xyz << std::endl;
            test::mtx_show_full(D_ao[xyz], log);
        }
        log << "LO basis coefficient matrix:" << std::endl;
        test::mtx_show_full(C_lo_basis, log);
        log << "DFA Hamiltonian under AO:" << std::endl;
        test::mtx_show_full(H_ao, log);
        log << "LO coefficient matrix: reference" << std::endl;
        test::mtx_show_full(C_lo_ref, log);
        log << "LO coefficient matrix: calculated" << std::endl;
        test::mtx_show_full(C_lo_calc, log);
        log << "U calculated is unitary:" << is_unitary(U_calc) << std::endl;
        log << "U reference is unitary:" << is_unitary(U_ref) << std::endl;
        log << "U matrix: reference" << std::endl;
        test::mtx_show_full(U_ref, log);
        log << "U matrix: calculated" << std::endl;
        test::mtx_show_full(U_calc, log);
        if (!is_permutation) {
            log << "Permutation matrix: calculated" << std::endl;
            test::mtx_show_full(P, log);
        }
        log.close();
        std::cout << "Test Failure: See log file for more information: "
                  << log_file << std::endl;
    }
}

std::vector<CaseInfo> test_targets{
    // ==> LO = CO <==
    // To check if localization is right for many orbitals.
    {
        .dir = "H2.1A",
        .label = "",
        .gamma = 0.707,
    },
    {
        .dir = "H2O",
        .label = "gamma=default.win=0_7",
        .gamma = 0.707,
    },
    // ==> LO != CO <==
    // stretched H2. The LOs are unique.
    {
        .dir = "H2.10A",
        .label = "",
        .gamma = 0.707,
    },
    // Localization of just 2 orbitals. This checks if the one-step rotation
    // angle is right.
    {
        .dir = "H2O",
        .label = "gamma=0.win=5_7",
        .gamma = 0.0,
    },
    // Localization of occupied orbitals + a few unoccupied orbitals.
    // The solution of localization for mostly occupied orbitals should be
    // unique. This case checks the correctness of localization for multiple
    // orbitals.
    {
        .dir = "H2O",
        .label = "gamma=0.win=0_7",
        .gamma = 0.0,
    },
    // Localization of ALL the orbitals (including a lot of virtual orbitals).
    // This is a crazy case. Because of the large number of virtual orbitals,
    // the solution of localization may not be numerically stable.
    // After testing with complation by cc/debug, cc/release, icc/debug
    // and icc/release, it looks like the results for this case is stable.
    // So keep it here for testing purpose.
    {
        .dir = "H2O",
        .label = "gamma=0.win=all",
        .gamma = 0.0,
    },
};

INSTANTIATE_TEST_SUITE_P(Localizationtest, LocalizationTest,
                         ::testing::ValuesIn(test_targets));
