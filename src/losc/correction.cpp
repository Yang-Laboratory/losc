#include "eigen_helper.hpp"
#include <losc/correction.hpp>
#include <losc/exception.hpp>
#include <sstream>
#include <utility>

namespace losc {

void C_API_ao_hamiltonian_correction(ConstRefMat &S, ConstRefMat &C_lo,
                                     ConstRefMat &Curvature,
                                     ConstRefMat &LocalOcc, RefMat H_losc)
{
    const size_t nlo = C_lo.cols();
    const size_t nbasis = C_lo.rows();

    if (nlo > nbasis) {
        throw exception::DimensionError(
            "losc::ao_hamiltonian_correction(): the number of LOs is larger "
            "than the number of AOs.");
    }
    if (!mtx_match_dimension(S, nbasis, nbasis)) {
        throw exception::DimensionError(
            S, nbasis, nbasis,
            "losc::ao_hamiltonian_correction(): mismatch AO overlap matrix.");
    }
    if (!mtx_match_dimension(Curvature, nlo, nlo)) {
        throw exception::DimensionError(
            Curvature, nlo, nlo,
            "losc::ao_hamiltonian_correction(): mismatch curvature matrix.");
    }
    if (!mtx_match_dimension(LocalOcc, nlo, nlo)) {
        throw exception::DimensionError(LocalOcc, nlo, nlo,
                                        "losc::ao_hamiltonian_correction(): "
                                        "mismatch local occupation matrix.");
    }

    // build A matrix.
    MatrixXd A(nlo, nlo);
    for (size_t i = 0; i < nlo; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            const double K_ij = Curvature(i, j);
            const double L_ij = LocalOcc(i, j);
            if (i != j) {
                A(i, j) = -K_ij * L_ij;
                A(j, i) = A(i, j);
            } else {
                A(i, j) = 0.5 * K_ij - K_ij * L_ij;
            }
        }
    }

    H_losc = S * C_lo * A * C_lo.transpose() * S;
}

MatrixXd ao_hamiltonian_correction(ConstRefMat &S, ConstRefMat &C_lo,
                                   ConstRefMat &Curvature,
                                   ConstRefMat &LocalOcc)
{
    const size_t nbasis = S.cols();
    MatrixXd H_losc(nbasis, nbasis);
    C_API_ao_hamiltonian_correction(S, C_lo, Curvature, LocalOcc, H_losc);
    return std::move(H_losc);
}

double energy_correction(ConstRefMat &Curvature, ConstRefMat &LocalOcc)
{
    const size_t nlo = Curvature.rows();
    if (!mtx_match_dimension(Curvature, nlo, nlo)) {
        throw exception::DimensionError(
            Curvature, nlo, nlo,
            "losc::energy_correction(): mismatch curvature matrix.");
    }
    if (!mtx_match_dimension(LocalOcc, nlo, nlo)) {
        throw exception::DimensionError(LocalOcc, nlo, nlo,
                                        "losc::energy_correction(): "
                                        "mismatch local occupation matrix.");
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

void C_API_orbital_energy_post_scf(ConstRefMat &H_dfa, ConstRefMat &H_losc,
                                   ConstRefMat &C_co, double *eig)
{
    const size_t nbasis = C_co.rows();
    const size_t nlo = C_co.cols();

    if (!mtx_match_dimension(H_dfa, nbasis, nbasis)) {
        throw exception::DimensionError(H_dfa, nbasis, nbasis,
                                        "losc::post_scf_orbital_energy(): "
                                        "mismatch DFA Hamiltonian matrix.");
    }
    if (!mtx_match_dimension(H_losc, nbasis, nbasis)) {
        throw exception::DimensionError(
            H_losc, nbasis, nbasis,
            "losc::post_scf_orbital_energy(): "
            "mismatch LOSC effective Hamiltonian matrix.");
    }

    MatrixXd H_tot = H_dfa + H_losc;
    for (size_t i = 0; i < nlo; ++i) {
        eig[i] = C_co.col(i).transpose() * H_tot * C_co.col(i);
    }
}

vector<double> orbital_energy_post_scf(ConstRefMat &H_dfa, ConstRefMat &H_losc,
                                       ConstRefMat &C_co)
{
    const size_t n = C_co.cols();
    vector<double> eig(n, 0);
    C_API_orbital_energy_post_scf(H_dfa, H_losc, C_co, eig.data());
    return std::move(eig);
}

} // namespace losc
