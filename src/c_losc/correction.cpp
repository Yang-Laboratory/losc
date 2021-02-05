#include "matrix_impl.hpp"
#include <c_losc/correction.h>
#include <losc/correction.hpp>
#include <stdlib.h>

extern "C" {

losc_matrix *losc_ao_hamiltonian_correction(const losc_matrix *S,
                                            const losc_matrix *C_lo,
                                            const losc_matrix *Curvature,
                                            const losc_matrix *LocalOcc)
{
    const size_t nbasis = S->col_;
    losc_matrix *H_losc = losc_matrix_create(nbasis, nbasis);
    losc::C_API_ao_hamiltonian_correction(
        losc_matrix_to_eigen_const(S), losc_matrix_to_eigen_const(C_lo),
        losc_matrix_to_eigen_const(Curvature),
        losc_matrix_to_eigen_const(LocalOcc), losc_matrix_to_eigen(H_losc));
    return H_losc;
}

double energy_correction(const losc_matrix *Curvature,
                         const losc_matrix *LocalOcc)
{
    return losc::energy_correction(losc_matrix_to_eigen_const(Curvature),
                                   losc_matrix_to_eigen_const(LocalOcc));
}

double *orbital_energy_post_scf(const losc_matrix *H_dfa,
                                const losc_matrix *H_losc,
                                const losc_matrix *C_co)
{
    const size_t n = C_co->col_;
    double *eig = (double *)calloc(n, sizeof(double));
    losc::C_API_orbital_energy_post_scf(losc_matrix_to_eigen_const(H_dfa),
                                        losc_matrix_to_eigen_const(H_losc),
                                        losc_matrix_to_eigen_const(C_co), eig);
    return eig;
}
}
