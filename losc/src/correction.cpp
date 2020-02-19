/**
 * @file
 * @brief definition relates to Losc correction.
 */
#include "correction.h"

#include "blas_base.h"
#include "exception.h"
#include <matrix/matrix.h>
#include <memory>
#include <sstream>
#include <vector>

namespace losc {

/**
 * @par Formula
 * In matrix express, the Losc correcting Hamiltonian matrix under AO
 * is defined as \f[ H_{losc} = S * C^T  * A * C * S, \f]
 * where matrix S is the input `S` matrix, matrix C is the input
 * `C_lo` matrix and matrix A is defined as
 * \f[ A_{ij} = \delta_{ij} 1/2 K_{ii} - K_{ij} \lambda_{ij}, \f]
 * in which  `K` is the input `Curvature` matrix, \f$\lambda\f$ is
 * the input `LocalOcc` matrix.
 */
SharedMatrix losc_hamiltonian_correction(const Matrix &S, const Matrix &C_lo,
                                         const Matrix &Curvature,
                                         const Matrix &LocalOcc)
{
    const size_t nlo = C_lo.row();
    const size_t nbasis = C_lo.col();

    if (nlo > nbasis) {
        throw exception::DimensionError(
            "wrong dimension for LO coefficient matrix: number of LO is larger "
            "than the number of AO.");
    }
    if (!S.is_square() || nbasis != S.row()) {
        throw exception::DimensionError(
            S, nbasis, nbasis, "wrong dimension for AO overlap matrix.");
    }
    if (!Curvature.is_square() || nlo != Curvature.row()) {
        throw exception::DimensionError(
            Curvature, nlo, nlo, "wrong dimension for curvature matrix.");
    }
    if (!LocalOcc.is_square() || nlo != LocalOcc.row()) {
        throw exception::DimensionError(
            LocalOcc, nlo, nlo, "wrong dimension for local occupation matrix.");
    }

    Matrix A(nlo, nlo);
    for (size_t i = 0; i < nlo; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            const double K_ij = Curvature(i, j);
            const double L_ij = LocalOcc(i, j);
            if (i != j) {
                A(i, j) = -K_ij * L_ij;
            } else {
                A(i, j) = 0.5 * K_ij - K_ij * L_ij;
            }
        }
    }
    A.to_symmetric("L");

    SharedMatrix H = std::make_shared<Matrix>(nbasis, nbasis);
    Matrix CAC(nbasis, nbasis);
    matrix::mult_dgemm_ATBA(C_lo, A, CAC);
    matrix::mult_dgemm_ATBA(S, CAC, *H);
    return H;
}

double losc_total_energy_correction(const Matrix &Curvature,
                                    const Matrix &LocalOcc)
{
    const size_t nlo = Curvature.row();

    if (!Curvature.is_square() || nlo != Curvature.row()) {
        throw exception::DimensionError(
            Curvature, nlo, nlo, "wrong dimension for curvature matrix.");
    }
    if (!LocalOcc.is_square() || nlo != LocalOcc.row()) {
        throw exception::DimensionError(
            LocalOcc, nlo, nlo, "wrong dimension for local occupation matrix.");
    }

    double energy = 0.0;
    for (size_t i = 0; i < nlo; ++i) {
        energy +=
            0.5 * Curvature(i, i) * LocalOcc(i, i) * (1.0 - LocalOcc(i, i));
        for (size_t j = 0; j < i; ++j) {
            energy -= Curvature(i, j) * LocalOcc(i, j) * LocalOcc(i, j);
        }
    }
    return energy;
}

vector<double> losc_orbital_energy_correction(const Matrix &S,
                                              const Matrix &C_co,
                                              const Matrix &C_lo,
                                              const Matrix &Curvature,
                                              const Matrix &LocalOcc)
{
    // V_ij = <LO_i|CO_j> = C_LO x S x C_CO. If the LOs are expanded under COs,
    // the V matrix is the localization U matrix. Otherwise, V matrix has to be
    // calculated and used to construct the direct orbital energy correction.
    // Here, to make code more general, V matrix is always calculated and used.
    const size_t nlo = C_lo.row();
    const size_t nbasis = C_lo.col();

    if (nlo > nbasis) {
        throw exception::DimensionError(
            "wrong dimension for LO coefficient matrix: number of LO is larger "
            "than the number of AO.");
    }
    if (!S.is_square() || nbasis != S.row()) {
        throw exception::DimensionError(
            S, nbasis, nbasis, "wrong dimension for AO overlap matrix.");
    }
    if (nlo != C_co.row() || nbasis != C_co.col()) {
        throw exception::DimensionError(
            C_co, nbasis, nbasis, "wrong dimension for CO coefficient matrix.");
    }
    if (!Curvature.is_square() || nlo != Curvature.row()) {
        throw exception::DimensionError(
            Curvature, nlo, nlo, "wrong dimension for curvature matrix.");
    }
    if (!LocalOcc.is_square() || nlo != LocalOcc.row()) {
        throw exception::DimensionError(
            LocalOcc, nlo, nlo, "wrong dimension for local occupation matrix.");
    }

    vector<double> eig_correction(nlo);
    Matrix Clo_S(nlo, nbasis);
    Matrix V(nlo, nlo);
    matrix::mult_dgemm(1.0, C_lo, "N", S, "N", 0.0, Clo_S);
    matrix::mult_dgemm(1.0, Clo_S, "N", C_co, "T", 0.0, V);

    // get correction.
    for (size_t m = 0; m < nlo; ++m) {
        for (size_t i = 0; i < nlo; ++i) {
            const double K_ii = Curvature(i, i);
            const double L_ii = LocalOcc(i, i);
            const double V_im = V(i, m);
            eig_correction[m] += K_ii * (0.5 - L_ii) * V_im * V_im;
            for (size_t j = 0; j < i; j++) {
                const double V_jm = V(j, m);
                const double K_ij = Curvature(i, j);
                const double L_ij = LocalOcc(i, j);
                eig_correction[m] -= 2.0 * K_ij * L_ij * V_im * V_jm;
            }
        }
    }
    return eig_correction;
}

vector<double> losc_corrected_orbital_energy_by_projection(const Matrix &H_dfa,
                                                           const Matrix &H_losc,
                                                           const Matrix &C_co)
{
    const size_t nlo = C_co.row();
    const size_t nbasis = C_co.col();

    if (nlo != C_co.row() || nbasis != C_co.col()) {
        throw exception::DimensionError(
            C_co, nbasis, nbasis, "wrong dimension for CO coefficient matrix.");
    }
    if (!H_dfa.is_square() || nbasis != H_dfa.row()) {
        throw exception::DimensionError(
            H_dfa, nbasis, nbasis,
            "wrong dimension for DFA Hamiltonian matrix.");
    }
    if (!H_losc.is_square() || nbasis != H_losc.row()) {
        throw exception::DimensionError(
            H_losc, nbasis, nbasis,
            "wrong dimension for losc correcting Hamiltonian matrix.");
    }

    vector<double> eig(nlo);
    Matrix H_tot(nbasis, nbasis);
    for (size_t i = 0; i < H_dfa.size(); ++i) {
        H_tot.data()[i] = H_dfa.data()[i] + H_losc.data()[i];
    }

    Matrix Cco_H(nlo, nbasis);
    matrix::mult_dgemm(1.0, C_co, "N", H_tot, "N", 0.0, Cco_H);
    int nbasis_t = nbasis;
    for (size_t i = 0; i < nlo; ++i) {
        eig[i] = losc::blas::ddot_(&nbasis_t, C_co.data() + i * nbasis,
                                   losc::blas::ione, Cco_H.data() + i * nbasis,
                                   losc::blas::ione);
    }
    return eig;
}

vector<double> losc_corrected_orbital_energy_by_diagonalize(
    const Matrix &H_dfa, const Matrix &H_losc, const Matrix &S_neg_half)
{
    const size_t nbasis = H_dfa.row();

    if (!S_neg_half.is_square()) {
        throw exception::DimensionError(S_neg_half, nbasis, nbasis,
                                        "wrong dimension for S^(-1/2) matrix");
    }
    if (!H_dfa.is_square() || nbasis != H_dfa.row()) {
        throw exception::DimensionError(
            H_dfa, nbasis, nbasis,
            "wrong dimension for DFA Hamiltonian matrix.");
    }
    if (!H_losc.is_square() || nbasis != H_losc.row()) {
        throw exception::DimensionError(
            H_losc, nbasis, nbasis,
            "wrong dimension for losc correcting Hamiltonian matrix.");
    }

    vector<double> eig(nbasis);

    SharedMatrix H_tot = std::make_shared<Matrix>(nbasis, nbasis);
    for (size_t i = 0; i < H_dfa.size(); ++i) {
        H_tot->data()[i] = H_dfa.data()[i] + H_losc.data()[i];
    }

    // diagonal S^(-1/2) * H * S^(-1/2).
    try {
        SharedMatrix SH = std::make_shared<Matrix>(nbasis, nbasis);
        SharedMatrix SHS = std::make_shared<Matrix>(nbasis, nbasis);
        matrix::mult_dgemm(1.0, S_neg_half, "N", *H_tot, "N", 0.0, *SH);
        matrix::mult_dgemm(1.0, *SH, "N", S_neg_half, "N", 0.0, *SHS);
        SH.reset();
        matrix::diagonalize_sym_matrix_dsyev("L", *SHS, eig);
    } catch (matrix::exception::MatrixException &e) {
        std::stringstream msg;
        msg << "failed to diagonalize losc corrected total Hamiltonian to get "
               "eigenvalues.\n"
            << e.what();
        throw exception::LoscException(msg.str());
    }
    return eig;
}

} // namespace losc
