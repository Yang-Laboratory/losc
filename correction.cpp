#include "correction.h"
#include <matrix/matrix.h>
#include <memory>
#include <vector>

namespace losc {

/**
 * calculate losc correction to Hamiltonian.
 *
 * H_eff = S * C^T  * A * C * S:
 * S: AO overlap matrix.
 * C: LO coef matrix.
 * A: A_ij = \delta_{ij} 1/2 K_ii - K_ij \lambda_ij, K is kappa, \lambda is local occ.
 */
SharedMatrix losc_hamiltonian_matrix(const Matrix &S, const Matrix &C_lo, const Matrix &Curvature,
                                     const Matrix &LocalOcc)
{
    const size_t nlo = C_lo.row();
    const size_t nbasis = C_lo.col();
    Matrix A(nlo, nbasis);
    for (size_t i = 0; i < nlo; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            const double K_ij = Curvature(i, j);
            const double L_ij = LocalOcc(i, j);
            if (i != j) {
                A(i, j) = - K_ij * L_ij;
            } else {
                A(i, j) = 0.5 * K_ij - K_ij * L_ij;
            }
        }
    }

    SharedMatrix H = std::make_shared<Matrix>(nbasis, nbasis);
    Matrix CS(nlo, nbasis);
    Matrix ACS(nlo, nbasis);
    matrix::mult_dgemm(1.0, C_lo, "N", S, "N", 0.0, CS);
    matrix::mult_dgemm(1.0, A, "N", CS, "N", 0.0, ACS);
    matrix::mult_dgemm(1.0, CS, "T", ACS, "N", 0.0, *H);
    return H;
}

/**
 * calculate losc correction to total energy.
 */
double losc_energy(const Matrix &Curvature, const Matrix &LocalOcc)
{
    const size_t nlo = Curvature.row();

    double energy = 0.0;
    for (size_t i = 0; i < nlo; ++i) {
        energy += 0.5 * Curvature(i, i) * LocalOcc(i, i) * (1.0 - LocalOcc(i, i));
        for (size_t j = 0; j < i; ++j) {
            energy -= Curvature(i, j) * LocalOcc(i, j) * LocalOcc(i, j);
        }
    }
    return energy;
}

/**
 * calculate losc correction to orbital energy with direct correction.
 */
vector<double> losc_orbital_energy_direct_correction(const Matrix &S, const Matrix &C_co,
                                                     const Matrix &C_lo, const Matrix &Curvature,
                                                     const Matrix &LocalOcc)
{
    // V_ij = <LO_i|CO_j> = C_LO x S x C_CO. If the LOs are expanded under COs,
    // the V matrix is the localization U matrix. Otherwise, V matrix has to be
    // calculated and used to construct the direct orbital energy correction.
    // Here, to make code more general, V matrix is always calculated and used.

    // !!!!!!
    // make sure dimension match between C_lo and C_co.

    const size_t nlo = C_lo.row();
    const size_t nbasis = C_lo.col();

    vector<double> eig_correction(nlo);
    SharedMatrix Clo_S = std::make_shared<Matrix>(nlo, nbasis);
    SharedMatrix V = std::make_shared<Matrix>(nlo, nlo);
    matrix::mult_dgemm(1.0, C_lo, "N", S, "N", 0.0, *Clo_S);
    matrix::mult_dgemm(1.0, *Clo_S, "N", C_co, "T", 0.0, *V);

    // get correction.
    for (size_t m = 0; m < nlo; ++m) {
        for (size_t i = 0; i < nlo; ++i) {
            const double K_ii = Curvature(i, i);
            const double L_ii = LocalOcc(i, i);
            const double V_im = (*V)(i, m);
            eig_correction[m] += K_ii * (0.5 - L_ii) * V_im * V_im;
            for (size_t j = 0; j < i; j++) {
                const double V_jm = (*V)(j, m);
                const double K_ij = Curvature(i, j);
                const double L_ij = LocalOcc(i, j);
                eig_correction[m] -= 2.0 * K_ij * L_ij * V_im * V_jm;
            }
        }
    }
    return eig_correction;
}

/**
 * calculate orbital energy by projection.
 */
vector<double> losc_orbital_energy_by_projection(const Matrix &H_dfa, const Matrix &H_losc,
                                                 const Matrix &C_co)
{
    const size_t nlo = C_co.row();
    const size_t nbasis = C_co.col();
    vector<double> eig(nlo);

    SharedMatrix H_tot = std::make_shared<Matrix>(nbasis, nbasis);
    for (size_t i = 0; i < H_dfa.size(); ++i) {
        H_tot->data()[i] = H_dfa.data()[i] + H_losc.data()[i];
    }

    SharedMatrix Cco_H = std::make_shared<Matrix>(nlo, nbasis);
    matrix::mult_dgemm(1.0, C_co, "N", *H_tot, "N", 0.0, *Cco_H);
    int nbasis_t = nbasis;
    for (size_t i = 0; i < nlo; ++i) {
        eig[i] = matrix::blas::ddot_(&nbasis_t, C_co.data() + i * nbasis,
                                     matrix::blas::ione,
                                     Cco_H->data() + i * nbasis,
                                     matrix::blas::ione);
    }
    return eig;
}

/**
 * calculate orbital energy by diagnoalization.
 */
vector<double> losc_orbital_energy_by_diagonalize(const Matrix &Shalf, const Matrix &H_dfa,
                                                  const Matrix &H_losc)
{
    const size_t nbasis = H_dfa.row();
    vector<double> eig(nbasis);

    SharedMatrix H_tot = std::make_shared<Matrix>(nbasis, nbasis);
    for (size_t i = 0; i < H_dfa.size(); ++i) {
        H_tot->data()[i] = H_dfa.data()[i] + H_losc.data()[i];
    }

    // diagonal S^(-1/2) * H * S^(-1/2).
    SharedMatrix SH = std::make_shared<Matrix>(nbasis, nbasis);
    SharedMatrix SHS = std::make_shared<Matrix>(nbasis, nbasis);
    matrix::mult_dgemm(1.0, Shalf, "N", *H_tot, "N", 0.0, *SH);
    matrix::mult_dgemm(1.0, *SH, "N", Shalf, "N", 0.0, *SHS);
    SH.reset();

    matrix::diagonalize_sym_matrix_dsyev("L", *SHS, eig);
    return eig;
}

}
