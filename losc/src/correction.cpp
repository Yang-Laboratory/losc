/**
 * @file
 * @brief definition relates to Losc correction.
 */
#include "correction.h"

#include "blas_base.h"
#include "exception.h"
#include "matrix.h"
#include <Eigen/Eigenvalues>
#include <sstream>

namespace losc {

/**
 * @par Formula
 * In matrix express, the Losc correcting Hamiltonian matrix under AO
 * is defined as \f[ H_{losc} = S * C  * A * C^T * S, \f]
 * where matrix S is the input `S` matrix, matrix C is the input
 * `C_lo` matrix and matrix A is defined as
 * \f[ A_{ij} = \delta_{ij} 1/2 K_{ii} - K_{ij} \lambda_{ij}, \f]
 * in which  `K` is the input `Curvature` matrix, \f$\lambda\f$ is
 * the input `LocalOcc` matrix.
 */
shared_ptr<Matrix> losc_hamiltonian_correction(const Matrix &S,
                                               const Matrix &C_lo,
                                               const Matrix &Curvature,
                                               const Matrix &LocalOcc)
{
    const size_t nlo = C_lo.cols();
    const size_t nbasis = C_lo.rows();

    if (nlo > nbasis) {
        throw exception::DimensionError(
            "wrong dimension for LO coefficient matrix: number of LO is larger "
            "than the number of AO.");
    }
    if (!S.is_square() || nbasis != S.rows()) {
        throw exception::DimensionError(
            S, nbasis, nbasis, "wrong dimension for AO overlap matrix.");
    }
    if (!Curvature.is_square() || nlo != Curvature.rows()) {
        throw exception::DimensionError(
            Curvature, nlo, nlo, "wrong dimension for curvature matrix.");
    }
    if (!LocalOcc.is_square() || nlo != LocalOcc.rows()) {
        throw exception::DimensionError(
            LocalOcc, nlo, nlo, "wrong dimension for local occupation matrix.");
    }

    // build A matrix.
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

    // calculate Losc correcting Hamiltonian matrix.
    shared_ptr<Matrix> H = std::make_shared<Matrix>(nbasis, nbasis);
    (*H).noalias() = S * C_lo * A * C_lo.transpose() * S;
    return H;
}

double losc_total_energy_correction(const Matrix &Curvature,
                                    const Matrix &LocalOcc)
{
    const size_t nlo = Curvature.rows();

    if (!Curvature.is_square() || nlo != Curvature.rows()) {
        throw exception::DimensionError(
            Curvature, nlo, nlo, "wrong dimension for curvature matrix.");
    }
    if (!LocalOcc.is_square() || nlo != LocalOcc.rows()) {
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
    // V_ij = <LO_i|CO_j> = C_LO^T x S x C_CO. If the LOs are expanded under
    // COs, the V matrix is the localization U matrix. Otherwise, V matrix has
    // to be calculated and used to construct the direct orbital energy
    // correction. Here, to make code more general, V matrix is always
    // calculated and used.
    const size_t nlo = C_lo.cols();
    const size_t nbasis = C_lo.rows();

    if (nlo > nbasis) {
        throw exception::DimensionError(
            "wrong dimension for LO coefficient matrix: number of LO is larger "
            "than the number of AO.");
    }
    if (!S.is_square() || nbasis != S.rows()) {
        throw exception::DimensionError(
            S, nbasis, nbasis, "wrong dimension for AO overlap matrix.");
    }
    if (nlo != C_co.cols() || nbasis != C_co.rows()) {
        throw exception::DimensionError(
            C_co, nbasis, nbasis, "wrong dimension for CO coefficient matrix.");
    }
    if (!Curvature.is_square() || nlo != Curvature.rows()) {
        throw exception::DimensionError(
            Curvature, nlo, nlo, "wrong dimension for curvature matrix.");
    }
    if (!LocalOcc.is_square() || nlo != LocalOcc.cols()) {
        throw exception::DimensionError(
            LocalOcc, nlo, nlo, "wrong dimension for local occupation matrix.");
    }

    Matrix V(nlo, nlo);
    V = C_lo.transpose() * S * C_co;

    vector<double> eig_correction(nlo, 0.0);
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
    const size_t nlo = C_co.cols();
    const size_t nbasis = C_co.rows();

    if (nlo != C_co.cols() || nbasis != C_co.rows()) {
        throw exception::DimensionError(
            C_co, nbasis, nbasis, "wrong dimension for CO coefficient matrix.");
    }
    if (!H_dfa.is_square() || nbasis != H_dfa.rows()) {
        throw exception::DimensionError(
            H_dfa, nbasis, nbasis,
            "wrong dimension for DFA Hamiltonian matrix.");
    }
    if (!H_losc.is_square() || nbasis != H_losc.rows()) {
        throw exception::DimensionError(
            H_losc, nbasis, nbasis,
            "wrong dimension for losc correcting Hamiltonian matrix.");
    }

    Matrix H_tot(nbasis, nbasis);
    H_tot.noalias() = H_dfa + H_losc;

    vector<double> eig(nlo, 0.0);
    for (size_t i = 0; i < nlo; ++i) {
        eig[i] = C_co.col(i).transpose() * H_tot * C_co.col(i);
    }
    return eig;
}

vector<double> losc_corrected_orbital_energy_by_diagonalize(
    const Matrix &H_dfa, const Matrix &H_losc, const Matrix &S_neg_half)
{
    const size_t nbasis = H_dfa.rows();

    if (!S_neg_half.is_square()) {
        throw exception::DimensionError(S_neg_half, nbasis, nbasis,
                                        "wrong dimension for S^(-1/2) matrix");
    }
    if (!H_dfa.is_square() || nbasis != H_dfa.rows()) {
        throw exception::DimensionError(
            H_dfa, nbasis, nbasis,
            "wrong dimension for DFA Hamiltonian matrix.");
    }
    if (!H_losc.is_square() || nbasis != H_losc.rows()) {
        throw exception::DimensionError(
            H_losc, nbasis, nbasis,
            "wrong dimension for losc correcting Hamiltonian matrix.");
    }

    // SHS = S^(-1/2) * H * S^(-1/2).
    Matrix SHS(nbasis, nbasis);
    SHS.noalias() = S_neg_half * (H_dfa + H_losc) * S_neg_half;
    Eigen::SelfAdjointEigenSolver<MatrixXd> eigensolver(SHS,
                                                        Eigen::EigenvaluesOnly);
    if (eigensolver.info() != Eigen::Success) {
        std::stringstream msg;
        msg << "failed to diagonalize losc corrected total Hamiltonian to get "
               "eigenvalues.\n";
        throw exception::LoscException(msg.str());
    }
    // copy eigenvalues to results.
    vector<double> eig(eigensolver.eigenvalues().data(),
                       eigensolver.eigenvalues().data() + nbasis);
    return eig;
}

} // namespace losc
