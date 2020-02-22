/**
 * @file
 * @brief Implementation relates to Losc localization version 2.
 */
#include "localization.h"

#include "blas_base.h"
#include "exception.h"
#include "matrix.h"
#include <algorithm>
#include <cmath>
#include <random>
#include <stdarg.h>
#include <string>

namespace losc {

/**
 * @brief Get the optimized rotation angle value for one pair of orbitals.
 */
static void js_optimize_one_pair(const size_t i, const size_t j,
                                 const double para_c, const double para_gamma,
                                 const vector<shared_ptr<Matrix>> &dipole_lo,
                                 const shared_ptr<Matrix> &hamiltonian_lo,
                                 double &theta_val, double &delta_val)
{
    // Spacial part: constant x1 and x2.
    double x1 = 0.0;
    double x2 = 0.0;
    double Dii_Djj = 0.0;
    double Dij = 0.0;
    for (size_t xyz = 0; xyz < 3; ++xyz) {
        const Matrix &DMat = *dipole_lo[xyz];
        Dii_Djj = DMat(i, i) - DMat(j, j);
        Dij = DMat(i, j);
        x1 += -0.5 * Dii_Djj * Dii_Djj + 2 * Dij * Dij;
        x2 -= 2 * Dij * Dii_Djj;
    }

    // Energy part: constant y1 and y2.
    double y1 = 0.0;
    double y2 = 0.0;
    double Eii_Ejj = 0.0;
    double Eij = 0.0;
    const Matrix &EMat = *hamiltonian_lo;
    Eii_Ejj = EMat(i, i) - EMat(j, j);
    Eij = EMat(i, j);
    y1 = -0.5 * Eii_Ejj * Eii_Ejj + 2 * Eij * Eij;
    y2 = -2 * Eij * Eii_Ejj;

    /* constant values */
    const double a1 = (1 - para_gamma) * x1 + para_gamma * para_c * y1;
    const double a2 = (1 - para_gamma) * x2 + para_gamma * para_c * y2;

    /* theta */
    double theta = 0.0;
    // atan2 returns a value in range [-pi, pi], which is of a
    // larger range than what atan gives ([-pi/2, pi/2]).
    // This is actually a desired behavior, as we may want the
    // derivative of U to be stable.
    // Then U should not jump from -pi/4 to pi/4 (or pi/4 to -pi/4).
    theta = 0.5 * atan2(-a1 - sqrt(a1 * a1 + a2 * a2), a2);
    // try not to get too large rotation angles.
    // however we don't want restrict it in [-pi/4, pi/4],
    // as that will cause derivative discontinuity for numerical dU/dP.
    if (theta > 3. * M_PI / 8.)
        theta -= M_PI / 2.;
    if (theta < -3. * M_PI / 8.)
        theta += M_PI / 2.;
    theta_val = theta;
    double delta = a1 * cos(4 * theta) + a2 * sin(4 * theta) - a1;
    delta_val = delta; // unit in bohr^2.
}

/**
 * @brief Apply rotation to the pair of orbitals.
 */
static void js_rotate_one_pair(const size_t i, const size_t j,
                               const double theta,
                               const shared_ptr<Matrix> &u_matrix,
                               const vector<shared_ptr<Matrix>> &dipole_lo,
                               const shared_ptr<Matrix> &hamiltonian_lo)
{
    double costheta = cos(theta);
    double sintheta = sin(theta);
    const int nLMO_t = u_matrix->rows();

    double *U_Mat = u_matrix->data();
    losc::blas::drot_(&nLMO_t, U_Mat + i * nLMO_t, losc::blas::ione,
                      U_Mat + j * nLMO_t, losc::blas::ione, &costheta,
                      &sintheta);

    // rotate Dipole Matrix
    for (int xyz = 0; xyz < 3; ++xyz) {
        double *DMat = dipole_lo[xyz]->data();
        losc::blas::drot_(&nLMO_t, DMat + i * nLMO_t, losc::blas::ione,
                          DMat + j * nLMO_t, losc::blas::ione, &costheta,
                          &sintheta);
        losc::blas::drot_(&nLMO_t, DMat + i, &nLMO_t, DMat + j, &nLMO_t,
                          &costheta, &sintheta);
    }

    // rotate localization Hamiltonian matrix.
    double *HMat = hamiltonian_lo->data();
    losc::blas::drot_(&nLMO_t, HMat + i * nLMO_t, losc::blas::ione,
                      HMat + j * nLMO_t, losc::blas::ione, &costheta,
                      &sintheta);
    losc::blas::drot_(&nLMO_t, HMat + i, &nLMO_t, HMat + j, &nLMO_t, &costheta,
                      &sintheta);
}

LoscLocalizerV2::LoscLocalizerV2(const shared_ptr<Matrix> &C_lo_basis,
                                 const shared_ptr<Matrix> &H_ao,
                                 const vector<shared_ptr<Matrix>> &Dipole_ao)
    : LocalizerBase(C_lo_basis), H_ao_{H_ao}, Dipole_ao_{Dipole_ao}
{
    if (!H_ao_->is_square() && H_ao_->rows() != nbasis_) {
        throw exception::DimensionError(
            *H_ao, nbasis_, nbasis_,
            "wrong dimension for DFA Hamiltonian matrix under AO.");
    }

    if (Dipole_ao_.size() != 3) {
        throw exception::DimensionError(
            "No enough dipole matrices under AO is given: x, y and z "
            "components are needed.");
        vector<std::string> xyz_name = {"x", "y", "z"};
        for (size_t i = 0; i < 3; ++i) {
            if (!Dipole_ao_[i]->is_square() &&
                Dipole_ao_[i]->rows() != nbasis_) {
                std::string msg =
                    "wrong dimension for dipole matrix under AO for " +
                    xyz_name[i] + "component.";
                throw exception::DimensionError(*Dipole_ao[i], nbasis_, nbasis_,
                                                msg);
            }
        }
    }
}

void LoscLocalizerV2::message(std::string t, ...) const
{
    if (print_level_ >= kPrintLevelNormal) {
        va_list args;
        va_start(args, t.c_str());
        vfprintf(stdout, t.c_str(), args);
        va_end(args);
    }
}

shared_ptr<Matrix> LoscLocalizerV2::compute()
{
    Matrix &U = *U_;
    Matrix &C_lo_basis = *C_lo_basis_;

    // calculate dipole on LO initial guess.
    // D_lo = U^T * C_lo_basis^T * D_ao * C_lo_basis * U
    vector<shared_ptr<Matrix>> D_lo;
    for (size_t xyz = 0; xyz < 3; xyz++) {
        D_lo.push_back(std::make_shared<Matrix>(nlo_, nlo_));
        Matrix &D_ao = *Dipole_ao_[xyz];
        (*D_lo[xyz]).noalias() =
            U.transpose() * C_lo_basis.transpose() * D_ao * C_lo_basis * U;
    }

    // calculate Hamiltonian matrix on LO initial guess.
    // H_lo = U^T * C_lo_basis^T * H_ao * C_lo_basis * U
    shared_ptr<Matrix> H_lo = std::make_shared<Matrix>(nlo_, nlo_);
    Matrix &H_ao = *H_ao_;
    (*H_lo).noalias() =
        U.transpose() * C_lo_basis.transpose() * H_ao * C_lo_basis * U;

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
                js_optimize_one_pair(i, j, para_c_, para_gamma_, D_lo, H_lo,
                                     theta, delta);
                // apply rotation to related matrix.
                js_rotate_one_pair(i, j, theta, U_, D_lo, H_lo);
                cycle_delta += delta;
            }
        }
        ++iter;
    }

    message("Localization final iteration = %zu.\n", iter);
    if (iter >= js_max_iter_ && abs(cycle_delta) > js_tol_) {
        std::cout << "Warning: localization is not convergened.\n";
    }

    // calculate the LO coefficient matrix.
    // C_lo = C_lo_basis * U
    auto C_lo = std::make_shared<Matrix>(nlo_, nlo_);
    (*C_lo).noalias() = C_lo_basis * U;
    return C_lo;
}

} // namespace losc
