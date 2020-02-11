#include "localization.h"
#include <matrix/matrix.h>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <random>

namespace localization {

/**
 * Get the optimized rotation angle value for one pair of orbitals.
 */
static void js_optimize_one_pair(const size_t i, const size_t j,
                                 const double para_c, const double para_gamma,
                                 vector<SharedMatrix> dipole_lo,
                                 SharedMatrix hamiltonian_lo,
                                 double &theta_val, double &delta_val)
{
    // Spacial part: constant x1 and x2.
    double x1 = 0.0;
    double x2 = 0.0;
    double Dii_Djj = 0.0;
    double Dij =0.0;
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
    double Eij =0.0;
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
    if (theta > 3. * M_PI / 8.) theta -= M_PI / 2.;
    if (theta < -3. * M_PI / 8.) theta += M_PI / 2.;
    theta_val = theta;
    double delta = a1 * cos(4 * theta) + a2 * sin(4 * theta) - a1;
    delta_val = delta; // unit in bohr^2.
}

/**
 * Apply rotation to the pair of orbitals.
 */
static void js_rotate_one_pair(const size_t i, const size_t j, const double theta,
                               SharedMatrix u_matrix,
                               vector<SharedMatrix> dipole_lo,
                               SharedMatrix hamiltonian_lo)
{
    double costheta = cos(theta);
    double sintheta = sin(theta);
    const int nLMO_t = u_matrix->row();

    // rotate U matrix
    double *U_Mat = u_matrix->data();
    matrix::blas::drot_(&nLMO_t, U_Mat + i * nLMO_t, matrix::blas::ione,
                        U_Mat + j * nLMO_t, matrix::blas::ione,
                        &costheta, &sintheta);

    // rotate Dipole Matrix
    for (int xyz = 0; xyz < 3; ++xyz) {
        double *DMat = dipole_lo[xyz]->data();
        matrix::blas::drot_(&nLMO_t, DMat + i * nLMO_t, matrix::blas::ione,
                            DMat + j * nLMO_t, matrix::blas::ione,
                            &costheta, &sintheta);
        matrix::blas::drot_(&nLMO_t, DMat + i, &nLMO_t, DMat + j, &nLMO_t,
                            &costheta, &sintheta);
    }

    // rotate localization Hamiltonian matrix.
    double *HMat = hamiltonian_lo->data();
    matrix::blas::drot_(&nLMO_t, HMat + i * nLMO_t, matrix::blas::ione,
                        HMat + j * nLMO_t, matrix::blas::ione,
                        &costheta, &sintheta);
    matrix::blas::drot_(&nLMO_t, HMat + i, &nLMO_t, HMat + j, &nLMO_t,
                        &costheta, &sintheta);
}

Losc2Localizer::Losc2Localizer(shared_ptr<const LOBasisCoefficientMatrix> C_lo_basis,
                               shared_ptr<const HamiltonianAOMatrix> H_ao,
                               vector<shared_ptr<const DipoleAOMatrix>> Dipole_ao)
    : LocalizerBase(C_lo_basis), H_ao_{H_ao}, Dipole_ao_{Dipole_ao}
{
    if (! H_ao_->is_square() && H_ao_->row() != nbasis_) {
        std::cout << "Dimension error: Hamiltonian under AO.\n";
        std::exit(EXIT_FAILURE);
    }

    if (Dipole_ao_.size() != 3) {
        std::cout << "Dipole size error!\n";
        std::exit(EXIT_FAILURE);
        for (size_t i = 0; i < 3; ++i) {
            if (! Dipole_ao_[i]->is_square() && Dipole_ao_[i]->row() != nbasis_) {
                std::cout << "Dimension error: Dipole under AO in xyz = " << i << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
    }
}

/**
 * Do localization and compute the LO coefficient matrix under AO.
 */
void Losc2Localizer::compute()
{
    // calculate dipole on LO initial guess.
    // D_lo = U * C_lo_basis * D_ao * C_lo_basis^T * U^T
    vector<SharedMatrix> D_lo;
    for (size_t xyz = 0; xyz < 3; xyz++) {
        D_lo.push_back(std::make_shared<Matrix> (nlo_, nlo_));
        SharedMatrix CDC = std::make_shared<Matrix> (nlo_, nlo_);
        matrix::mult_dgemm_ABAT(*C_lo_basis_, *Dipole_ao_[xyz], *CDC);
        matrix::mult_dgemm_ABAT(*U_, *CDC, *D_lo[xyz]);
    }

    // calculate Hamiltonian matrix on LO initial guess.
    // H_lo = U * C_lo_basis * H_ao * C_lo_basis^T * U^T
    SharedMatrix H_lo = std::make_shared<Matrix> (nlo_, nlo_);
    SharedMatrix CHC = std::make_shared<Matrix> (nlo_, nlo_);
    matrix::mult_dgemm_ABAT(*C_lo_basis_, *H_ao_, *CHC);
    matrix::mult_dgemm_ABAT(*U_, *CHC, *H_lo);
    CHC.reset();

    // using jacobi sweep for localization.
    std::mt19937 g(0);
    size_t iter = 0;
    double cycle_delta = 10000.0;
    while (iter < js_max_iter_ && fabs(cycle_delta) > js_tol_) {
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
                //JS_pair_operator(JS_pair_optimize, pArgs, is, i, j, JS_constant, &theta, &delta);
                js_optimize_one_pair(i, j, para_c_, para_gamma_,
                                     D_lo, H_lo, theta, delta);
                // apply rotation to related matrix.
                js_rotate_one_pair(i, j, theta, U_, D_lo, H_lo);
                cycle_delta += delta;
            }
        }
        iter++;
    }
    if (iter >= js_max_iter_ && fabs(cycle_delta) > js_tol_) {
        std::cout << "Warning: localization is not convergened.\n";
    }

    // calculate the LO coefficient matrix.
    C_lo_ = std::make_shared<LOCoefficientMatrix> (nlo_, nlo_);
    matrix::mult_dgemm(1.0, *U_, "N", *C_lo_basis_, "N", 0.0, *C_lo_);
}

}
