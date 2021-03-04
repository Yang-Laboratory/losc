/**
 * @file
 * @brief The C interface for LOSC local occupation number.
 */

#ifndef _LOSC_INCLUDE_C_LOSC_LOCAL_OCCUPATION_H_
#define _LOSC_INCLUDE_C_LOSC_LOCAL_OCCUPATION_H_

#include <c_losc/matrix.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Create the local occupation matrix.
 * @param [in] C_lo: LO coefficient matrix with dimension `[nbasis, nlo]`.
 * @param [in] S: AO overlap matrix with dimension `[nbasis, nbasis]`.
 * @param [in] D: Spin density matrix under AO with dimension
 * `[nbasis, nbasis]`.
 * @return Create a local occupation matrix with dimension `[nlo, nlo]` and
 * return the pointer.
 * @see LoscLocalizerBase.lo: Obtain the LOs' coefficient matrix.
 * @note The local occupation matrix is defined as
 * \f$\displaystyle \lambda_{ij} = \langle \phi_i|\rho|\phi_j \rangle \f$,
 * where \f$ \lambda_{ij} \f$ is the local occupation matrix
 * element, \f$ \phi_i \f$ and \f$ \phi_j \f$ are the i-th and j-th localized
 * orbital and \f$ \rho \f$ is the spin density operator.
 * In matrix form, local occupation matrix \f$ L \f$ is expressed as
 * \f$ L = C_{\rm{LO}}^T S D S C_{\rm{LO}} \f$.
 */
losc_matrix *losc_local_occupation(const losc_matrix *C_lo,
                                   const losc_matrix *S, const losc_matrix *D);

#ifdef __cplusplus
}
#endif

#endif
