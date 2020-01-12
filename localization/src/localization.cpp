#include "localization.h"
#include <matrix/matrix.h>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <random>


namespace localization {

static void js_one_pair_opt(size_t i, size_t j, double para_c, double para_gamma,
                            vector<Matrix *> dipole_lo, Matrix *hamiltonian_lo,
                            double *theta_val, double *delta_val)
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
    *theta_val = theta;
    double delta = a1 * cos(4 * theta) + a2 * sin(4 * theta) - a1;
    *delta_val = delta; // unit in bohr^2.
}

static void js_rotate_one_pair(size_t i, size_t j, double theta,
                               Matrix *u_matrix,
                               vector<Matrix *> &dipole_lo,
                               Matrix *hamiltonian_lo)
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

/**
 * Create a localized orbital coefficient matrix that expanded on atomic basis set
 * based on the input localization data.
 *
 * @ param[in] input_data: the input data for localization. On exit, the input_data.u_matrix
 *  is updated with final U matrix.
 * @ return Matrix: the LO coefficient matrix on atomic basis.
 */
Matrix *LocalizerLosc2(LocalizerLosc2InputData &input_data)
{
    const size_t nlmo = input_data.nlmo;
    const size_t nbasis = input_data.nbasis;
    const Matrix *lo_basis_coef = input_data.lo_basis_coef;
    const Matrix *hamiltonian_ao = input_data.hamiltonian_ao;
    vector<const Matrix *> &dipole_ao = input_data.dipole_ao;
    Matrix *u_matrix = input_data.u_matrix;

    // calculate dipole on LO initial guess.
    // D_lo = U * C_lo_basis * D_ao * C_lo_basis^T * U^T
    vector<Matrix *> dipole_lo = {nullptr, nullptr, nullptr};
    for (size_t xyz = 0; xyz < 3; xyz++) {
        dipole_lo[xyz] = new Matrix(nlmo, nlmo);
        Matrix *CDC = new Matrix(nlmo, nlmo);
        matrix::mult_dgemm_ABAT(*lo_basis_coef, *dipole_ao[xyz], *CDC);
        matrix::mult_dgemm_ABAT(*u_matrix, *CDC, *dipole_lo[xyz]);
        delete CDC;
    }

    // calculate Hamiltonian matrix on LO initial guess.
    // H_lo = U * C_lo_basis * H_ao * C_lo_basis^T * U^T
    Matrix *hamiltonian_lo = new Matrix(nlmo, nlmo);
    hamiltonian_lo = new Matrix(nlmo, nlmo);
    Matrix *CHC = new Matrix(nlmo, nlmo);
    matrix::mult_dgemm_ABAT(*lo_basis_coef, *hamiltonian_ao, *CHC);
    matrix::mult_dgemm_ABAT(*u_matrix, *CHC, *hamiltonian_lo);
    delete CHC;

    // using jacobi sweep for localization.
    std::mt19937 g(0);
    size_t iter = 0;
    double cycle_delta = 10000.0;
    while (iter < input_data.js_max_iter && fabs(cycle_delta) > input_data.js_tol) {
        cycle_delta = 0.0;
        vector<size_t> order;
        for (size_t i = 0; i < nlmo; i++)
            order.push_back(i);

        /* random permutation */
        if (input_data.use_js_random_permutation) {
            std::shuffle(order.begin(), order.end(), g);
        }

        /* orbital pair rotation */
        for (size_t ni = 0; ni < nlmo; ni++) {
            size_t i = order[ni];
            for (size_t nj = 0; nj < ni; nj++) {
                size_t j = order[nj];
                double delta = 0.0;
                double theta = 0.0;
                // get the optimized rotation angle.
                //JS_pair_operator(JS_pair_optimize, pArgs, is, i, j, JS_constant, &theta, &delta);
                js_one_pair_opt(i, j, input_data.para_c, input_data.para_gamma,
                                dipole_lo, hamiltonian_lo,
                                &theta, &delta);
                // apply rotation to related matrix.
                js_rotate_one_pair(i, j, theta, u_matrix, dipole_lo, hamiltonian_lo);
                cycle_delta += delta;
            }
        }
        iter++;
    }
    if (iter >= input_data.js_max_iter && fabs(cycle_delta) > input_data.js_tol) {
        std::cout << "Warning: localization is not convergened.\n";
    }

    for (size_t xyz = 0; xyz < 3; xyz++) {
        delete dipole_lo[xyz];
    }
    delete hamiltonian_lo;

    // create LO coefficient matrix.
    Matrix *lo_coef = new Matrix(nlmo, nbasis);
    matrix::mult_dgemm(1.0, *u_matrix, "N", *lo_basis_coef, "N",
                       0.0, *lo_coef);
    return lo_coef;
}

}
