/**
 * @file
 * @brief C++ interface for the local occupation matrix within LOSC.
 */

#ifndef _LOSC_INCLUDE_LOSC_LOCAL_OCCUPATION_HPP_
#define _LOSC_INCLUDE_LOSC_LOCAL_OCCUPATION_HPP_

#include <losc/eigen_def.hpp>

namespace losc {

/**
 * @class __param__LocalOcc
 * @param [in] LocalOcc LOSC local occupation matrix with dimension of
 * `[nlo, nlo]`.
 */

/**
 * @class __param__local_occupation
 * @copydoc __param__C_lo
 * @param [in] S: AO overlap matrix with dimension `[nbasis, nbasis]`.
 * @param [in] D: Spin density matrix under AO with dimension
 * `[nbasis, nbasis]`.
 */

/**
 * @class __see__note__local_occupation
 * @see losc::LocalizerV2::lo(): Obtain the LOs' coefficient matrix.
 * @note The local occupation matrix is defined as
 * \f$\displaystyle \lambda_{ij} = \langle \phi_i|\rho|\phi_j \rangle \f$,
 * where \f$ \lambda_{ij} \f$ is the local occupation matrix
 * element, \f$ \phi_i \f$ and \f$ \phi_j \f$ are the i-th and j-th localized
 * orbital and \f$ \rho \f$ is the spin density operator.
 * In matrix form, local occupation matrix \f$ L \f$ is expressed as
 * \f$ L = C_{\rm{LO}}^T S D S C_{\rm{LO}} \f$.
 */

/**
 * @brief Calculate the local occupation matrix.
 * @copydoc __param__local_occupation
 * @return The local occupation matrix with dimension `[nlo, nlo]`.
 * @copydoc __see__note__local_occupation
 */
LOSCMatrix local_occupation(ConstRefMat &C_lo, ConstRefMat &S, ConstRefMat &D);

/**
 * @brief C interface to calculate the local occupation matrix.
 * @copydoc __param__local_occupation
 * @param [in, out] LocalOcc A matrix with dimension `[nlo, nlo]`. At entrance,
 * its data is ignored. At exit, it stores the local occupation matrix.
 * @copydoc __see__note__local_occupation
 */
void C_API_local_occupation(ConstRefMat &C_lo, ConstRefMat &S, ConstRefMat &D,
                            RefMat LocalOcc);
} // namespace losc
#endif // _LOSC_SRC_LOCAL_OCCUPATION_H_
