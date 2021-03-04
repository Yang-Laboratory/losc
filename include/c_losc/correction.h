/**
 * @file
 * @brief The C interface for LOSC corrections.
 */

#ifndef _LOSC_INCLUDE_C_LOSC_CORRECTION_H_
#define _LOSC_INCLUDE_C_LOSC_CORRECTION_H_

#include <c_losc/matrix.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @class __param__Curvature
 * @param [in] Curvature LOSC curvature matrix with dimension of `[nlo, nlo]`.
 */

/**
 * @class __param__LocalOcc
 * @param [in] LocalOcc Local occupation matrix with dimension of `[nlo, nlo]`.
 */

/**
 * @class __param__C_lo
 * @param [in] C_lo LOs coefficient matrix under AOs with dimension of
 * `[nbasis, nlo]`. The `i`-th column in `C_lo` matrix is the `i`-th LO.
 */

/**
 * @brief Calculate LOSC effective Hamiltonian under AO basis.
 * @details The LOSC effective Hamiltonian is constructed with LOs fixed. The
 * expression of the effective Hamiltonian is shown as Eq. S25 in the
 * supporting information of the original LOSC paper
 * (https://doi.org/10.1093/nsr/nwx111). This effective Hamiltonian is exact
 * in the developed version of SCF-LOSC (self-consistent LOSC). See reference
 * (J. Phys. Chem. Lett. 2020, 11, 23, 10269-10277) for more details about
 * how to perform reliable SCF-LOSC calculations.
 * @param [in] S AO overlap matrix with dimension [nbasis, nbasis].
 * @copydoc __param__C_lo
 * @copydoc __param__Curvature
 * @copydoc __param__LocalOcc
 * @return Create a LOSC effective Hamiltonian matrix with dimension
 * `[nbasis, nbasis]` and return the matrix pointer.
 * @see LoscCurvatureBase.kappa: Obtain the curvature matrix.
 * @see losc_local_occupation(): Obtain the local occupation matrix.
 * @note Make sure all the input matrices are corrresponding to the same spin.
 */
losc_matrix *losc_ao_hamiltonian_correction(const losc_matrix *S,
                                            const losc_matrix *C_lo,
                                            const losc_matrix *Curvature,
                                            const losc_matrix *LocalOcc);

/**
 * @brief Calculate the total energy correction from LOSC.
 * @details This is just the energy correction from LOSC, NOT the total energy
 * of LOSC-DFA. Total energy of LOSC-DFA is: ``E_losc_dfa = E_dfa + E_losc``.
 * @copydoc __param__Curvature
 * @copydoc __param__LocalOcc
 * @return The total energy correction from LOSC.
 * @see LoscCurvatureBase.kappa: obtain the LOSC curvature matrix.
 * @see losc_local_occupation(): obtain the LOSC local occupation matrix.
 */
double energy_correction(const losc_matrix *Curvature,
                         const losc_matrix *LocalOcc);

/**
 * @brief Calculate corrected orbital energy from LOSC in a post-SCF approach.
 * @details This function gives the final orbital energies WITH the correction
 * from LOSC. Note the difference to the function energy_correction() that only
 * calculates the energy correction. The corrected orbital energies are the
 * expectation values of converged DFA's COs on the LOSC-DFA Hamiltonian,
 * that is,
 * \f[
 * \epsilon_i = \langle \psi_i | H_{\rm{dfa}} + H_{\rm{losc}} | \psi_i \rangle.
 * \f]
 * @param [in] H_dfa The DFA Hamiltonian under AOs with dimension of
 * `[nbasis, nbasis]`.
 * @param [in] H_losc The LOSC effective Hamiltonian under AOs with dimension
 * of `[nbasis, nbasis]`.
 * @param [in] C_co The coefficient matrix of converged DFA COs under AOs
 * with dimension of `[nbasis, n]` (`n <= nbasis`, which is the number of COs).
 * You can choose whatever number of COs you want.
 * @return Create an array of the corrected orbital energies from LOSC with
 * size of `n` and return the pointer. The order of orbital energies match the
 * order of input COs (order of columns in `C_co` matrix).
 * @see losc_ao_hamiltonian_correction(): obtain the LOSC effective Hamiltonian
 * under AOs.
 * @note This function only calculates the orbital energies for the single spin
 * associated with the input matrices.
 * @note This function is just one of the ways to construct the LOSC corrected
 * orbital energies in the post-SCF LOSC calculations. It is the way we used to
 * produce results in the published paper for the post-SCF LOSC calculations.
 * Besides this way, there are another two ways to calculate corrected orbital
 * energies: (1) diagonalize the corrected LOSC-DFA Hamiltonian;
 * (2) follow Eq. 11 of the original paper of LOSC to calculate the corrections
 * to orbital energies. These three ways usually produce very similar results.
 */
double *orbital_energy_post_scf(const losc_matrix *H_dfa,
                                const losc_matrix *H_losc,
                                const losc_matrix *C_co);

#ifdef __cplusplus
}
#endif

#endif
