/**
 * @file local_occupation.h
 * @brief C++ interface for the local occupation matrix within LOSC.
 */

#ifndef _LOSC_INCLUDE_LOSC_LOCAL_OCCUPATION_HPP_
#define _LOSC_INCLUDE_LOSC_LOCAL_OCCUPATION_HPP_

#include <losc/eigen_def.hpp>

namespace losc {

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
 * @return LOSCMatrix: The local occupation matrix with dimension [nlo, nlo].
 *
 * @see
 * LocalizerBase::lo(): Obtain the LOs' coefficient matrix.
 */
LOSCMatrix local_occupation(ConstRefMat &C_lo, ConstRefMat &S, ConstRefMat &D);
void C_API_local_occupation(ConstRefMat &C_lo, ConstRefMat &S, ConstRefMat &D,
                            RefMat LocalOcc);

} // namespace losc
#endif // _LOSC_SRC_LOCAL_OCCUPATION_H_
