#ifndef __LOSC_INTERFACE_C_LOSC_LOCAL_OCCUPATION_H__
#define __LOSC_INTERFACE_C_LOSC_LOCAL_OCCUPATION_H__

#include "matrix.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Calculate the local occupation matrix.
 *
 * @details The local occupation matrix is defined as
 * \f$\lambda_{ij} = \langle \phi_i|\rho|\phi_j \rangle,\f$
 * where \f$\lambda_{ij}\f$ is the local occupation matrix
 * element, \f$\phi_i\f$ and \f$\phi_j\f$ are the i-th and j-th localized
 * orbital and \f$\rho\f$ is the spin density operator. Note that \f$\lambda\f$,
 * \f$\phi\f$ and \f$\rho\f$ are for the same spin. In matrix form,
 * local occupation matrix \f$L\f$ is expressed as
 * \f$L = C_{LO}^T * S * D * S * C_{LO}\f$.
 *
 * @param [in] C_lo: LO coefficient matrix with dimension [nbasis, nlo].
 * @param [in] S: AO overlap matrix with dimension [nbasis, nbasis].
 * @param [in] D: Spin density matrix under AO with dimension [nbasis, nbasis].
 *
 * @return losc_matrix*: The local occupation matrix with dimension [nlo, nlo].
 *
 * @see
 * LoscLocalizerBase.lo(): Obtain the LOs' coefficient matrix.
 */
losc_matrix *losc_local_occupation(const losc_matrix *C_lo,
                                   const losc_matrix *S, const losc_matrix *D);

#ifdef __cplusplus
}
#endif

#endif
