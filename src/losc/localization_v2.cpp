/**
 * @file
 * @brief Implementation of Losc localization v2.
 */

#include "eigen_helper.hpp"
#include <algorithm>
#include <cmath>
#include <losc/exception.hpp>
#include <losc/localization.hpp>
#include <random>
#include <stdarg.h>
#include <string>

namespace losc {

void LocalizerV2::js_optimize_one_pair(const size_t i, const size_t j,
                                       const vector<LOSCMatrix> &D_lo,
                                       const LOSCMatrix &H_lo,
                                       double &theta_val,
                                       double &delta_val) const
{
    // Spacial part: constant x1 and x2.
    double x1 = 0.0;
    double x2 = 0.0;
    double Dii_Djj = 0.0;
    double Dij = 0.0;
    for (size_t xyz = 0; xyz < 3; ++xyz) {
        ConstRefMat D = D_lo[xyz];
        Dii_Djj = D(i, i) - D(j, j);
        Dij = D(i, j);
        x1 += -0.5 * Dii_Djj * Dii_Djj + 2 * Dij * Dij;
        x2 -= 2 * Dij * Dii_Djj;
    }

    // Energy part: constant y1 and y2.
    const double Hii_Hjj = H_lo(i, i) - H_lo(j, j);
    const double Hij = H_lo(i, j);
    const double y1 = -0.5 * Hii_Hjj * Hii_Hjj + 2 * Hij * Hij;
    const double y2 = -2 * Hij * Hii_Hjj;

    /* constant values */
    const double a1 = (1 - gamma_) * x1 + gamma_ * c_ * y1;
    const double a2 = (1 - gamma_) * x2 + gamma_ * c_ * y2;

    /* theta */
    double theta = 0.0;
    // atan2 returns a value in range [-pi, pi], which is of a
    // larger range than what atan gives ([-pi/2, pi/2]).
    // This is actually a desired behavior, as we may want the
    // derivative of U to be stable.
    // Then U should not jump from -pi/4 to pi/4 (or pi/4 to -pi/4).
    theta = 0.5 * atan2(-a1 - sqrt(a1 * a1 + a2 * a2), a2);
    if (theta > M_PI / 4.)
        theta -= M_PI / 2.;
    if (theta < -M_PI / 4.)
        theta += M_PI / 2.;
    theta_val = theta;
    double delta = a1 * cos(4 * theta) + a2 * sin(4 * theta) - a1;
    delta_val = delta; // unit in bohr^2.
}

void LocalizerV2::js_rotate_one_pair(const size_t i, const size_t j,
                                     const double theta, RefMat U,
                                     vector<LOSCMatrix> &D_lo,
                                     LOSCMatrix &H_lo) const
{
    // rotate U
    rotate_two_vectors(U.col(i), U.col(j), theta);

    // rotate D_lo
    for (int xyz = 0; xyz < 3; ++xyz) {
        LOSCMatrix &D = D_lo[xyz];
        rotate_two_vectors(D.row(i), D.row(j), theta);
        rotate_two_vectors(D.col(i), D.col(j), theta);
    }

    // rotate H_lo
    rotate_two_vectors(H_lo.row(i), H_lo.row(j), theta);
    rotate_two_vectors(H_lo.col(i), H_lo.col(j), theta);
}

LocalizerV2::LocalizerV2(ConstRefMat &C_lo_basis, ConstRefMat &H_ao,
                         const vector<RefConstMat> &Dipole_ao)
    : LocalizerBase(C_lo_basis), H_ao_{H_ao}, Dipole_ao_{Dipole_ao}
{
    if (!mtx_match_dimension(H_ao_, nbasis_, nbasis_)) {
        throw exception::DimensionError(H_ao_, nbasis_, nbasis_,
                                        "LocalizerV2: wrong dimension for "
                                        "DFA Hamiltonian matrix under AO.");
    }
    if (Dipole_ao_.size() < 3) {
        throw exception::DimensionError(
            "LocalizerV2: No enough dipole matrices under AO is given: x, "
            "y and z "
            "components are needed.");
        vector<std::string> xyz_name = {"x", "y", "z"};
        for (size_t i = 0; i < 3; ++i) {
            if (!mtx_match_dimension(Dipole_ao_[i], nbasis_, nbasis_)) {
                std::string msg = "LocalizerV2: wrong dimension for dipole "
                                  "matrix under AO for " +
                                  xyz_name[i] + "component.";
                throw exception::DimensionError(Dipole_ao_[i], nbasis_, nbasis_,
                                                msg);
            }
        }
    }
}

void LocalizerV2::C_API_lo_U(RefMat L, RefMat U)
{
    // Sanity check of the input.
    if (!mtx_match_dimension(L, nbasis_, nlo_)) {
        throw exception::DimensionError(
            L, nbasis_, nlo_,
            "LocalizerV2::lo_U: dimension error of LO coefficient matrix.");
    }
    if (!mtx_match_dimension(U, nlo_, nlo_)) {
        throw exception::DimensionError(
            U, nlo_, nlo_,
            "LocalizerV2::lo_U: dimension error of the U matrix.");
    }

    LOSCMatrix L_init = C_lo_basis_ * U;
    // calculate dipole on LO initial guess.
    // D_lo = U^T * C_lo_basis^T * D_ao * C_lo_basis * U
    vector<LOSCMatrix> D_lo;
    for (size_t xyz = 0; xyz < 3; xyz++) {
        D_lo.push_back(L_init.transpose() * Dipole_ao_[xyz] * L_init);
    }

    // calculate Hamiltonian matrix on LO initial guess.
    // H_lo = U^T * C_lo_basis^T * H_ao * C_lo_basis * U
    LOSCMatrix H_lo = L_init.transpose() * H_ao_ * L_init;

    // using jacobi sweep for localization.
    std::mt19937 g(0);
    size_t iter = 0;
    double cycle_delta = 10000.0;
    while (iter < js_max_iter_ && std::abs(cycle_delta) > js_tol_) {
        cycle_delta = 0.0;
        vector<size_t> order;
        for (size_t i = 0; i < nlo_; ++i)
            order.push_back(i);

        /* random permutation */
        if (js_random_permutation_) {
            std::shuffle(order.begin(), order.end(), g);
        }

        /* orbital pair rotation */
        for (size_t ni = 0; ni < nlo_; ++ni) {
            size_t i = order[ni];
            for (size_t nj = 0; nj < ni; ++nj) {
                size_t j = order[nj];
                double delta = 0.0;
                double theta = 0.0;
                // get the optimized rotation angle.
                // JS_pair_operator(JS_pair_optimize, pArgs, is, i, j,
                // JS_constant, &theta, &delta);
                js_optimize_one_pair(i, j, D_lo, H_lo, theta, delta);
                // apply rotation to related matrix.
                js_rotate_one_pair(i, j, theta, U, D_lo, H_lo);
                cycle_delta += delta;
            }
        }
        ++iter;
    }

    // modify the localization flags.
    if (iter >= js_max_iter_ && abs(cycle_delta) > js_tol_) {
        converged_ = false;
    } else {
        converged_ = true;
    }
    nsteps_ = iter;

    // calculate the LO coefficient matrix.
    // C_lo = C_lo_basis * U
    L.noalias() = C_lo_basis_ * U;
}

vector<LOSCMatrix> LocalizerV2::lo_U(const string &guess)
{
    vector<LOSCMatrix> rst{LOSCMatrix(nbasis_, nlo_), LOSCMatrix(nlo_, nlo_)};
    RefMat L = rst[0];
    RefMat U = rst[1];
    set_u_guess(U, guess);
    C_API_lo_U(L, U);
    return std::move(rst);
}

vector<LOSCMatrix> LocalizerV2::lo_U(ConstRefMat &U_guess, double threshold)
{
    vector<LOSCMatrix> rst{LOSCMatrix(nbasis_, nlo_), LOSCMatrix(nlo_, nlo_)};
    RefMat L = rst[0];
    RefMat U = rst[1];
    set_u_guess(U, U_guess, threshold);
    C_API_lo_U(L, U);
    return std::move(rst);
}

double LocalizerV2::cost_func(ConstRefMat &lo) const
{
    vector<LOSCMatrix> D;
    for (size_t xyz = 0; xyz < 3; ++xyz) {
        D.push_back(lo.transpose() * Dipole_ao_[xyz] * lo);
    }
    LOSCMatrix H(lo.transpose() * H_ao_ * lo);

    double rst = 0;
    for (size_t xyz = 0; xyz < 3; ++xyz) {
        rst -= (1 - gamma_) * D[xyz].diagonal().array().square().sum();
    }
    rst -= gamma_ * c_ * H.diagonal().array().square().sum();
    return rst;
}

} // namespace losc
