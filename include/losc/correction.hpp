/**
 * @file
 * @brief C++ interface for the correction from LOSC.
 */

#ifndef _LOSC_INCLUDE_LOSC_CORRECTION_HPP_
#define _LOSC_INCLUDE_LOSC_CORRECTION_HPP_

#include <losc/eigen_def.hpp>
#include <vector>

namespace losc {

using std::vector;

/**
 * @class __param__ao_hamiltonian_correction
 * @param [in] S AOs overlap matrix with dimension of `[nbasis, nbasis]`.
 * @copydoc __param__C_lo
 * @copydoc __param__Curvature
 * @copydoc __param__LocalOcc
 */

/**
 * @class __return__losc_Hamiltonian
 * @return the LOSC effective Hamiltonian with dimension of `[nbasis, nbasis]`.
 */

/**
 * @class __see__ao_hamiltonian_correction
 * @see losc::LocalizerV2::lo(): return the LOs coefficient matrix.
 * @see losc::CurvatureBase::kappa(): return the LOSC curvature matrix.
 * @see losc::local_occupation(): return the LOSC local occupation matrix.
 */

/**
 * @brief Calculate LOSC effective Hamiltonian under AO basis.
 * @details The LOSC effective Hamiltonian is used for self-consistent LOSC
 * (SCF-LOSC) calculations. This function constructs the LOSC effective
 * Hamiltonian with fixed LOs. This LOSC effective Hamiltonian is used as an
 * APPROXIMATION in the original LOSC paper
 * (see https://doi.org/10.1093/nsr/nwx111). In addition, this LOSC effective
 * Hamiltonian is used as the EXACT one in the developed version of the SCF-LOSC
 * (see J. Phys. Chem. Lett. 2020, 11, 23, 10269-10277).
 *
 * The expression of the LOSC effective Hamiltonian is shown as the Eq. S25 in
 * the supporting information of the paper of original LOSC, and the Eq. 12 in
 * the paper of the developed version.
 * @copydoc __param__ao_hamiltonian_correction
 * @copydoc __return__losc_Hamiltonian
 * @copydoc __see__ao_hamiltonian_correction
 */
LOSCMatrix ao_hamiltonian_correction(ConstRefMat &S, ConstRefMat &C_lo,
                                     ConstRefMat &Curvature,
                                     ConstRefMat &LocalOcc);

/**
 * @brief Calculate the total energy correction from LOSC.
 * @copydoc __param__Curvature
 * @copydoc __param__LocalOcc
 * @return The total energy correction from LOSC.
 * @see losc::CurvatureBase::kappa(): return the LOSC curvature matrix.
 * @see losc::local_occupation(): return the LOSC local occupation matrix.
 * @note This is just the energy correction from LOSC, NOT the total energy of
 * LOSC-DFA. Total energy of LOSC-DFA is ``E_losc_dfa = E_dfa + E_losc``.
 * @note This function only calculates the energy correction of the single spin
 * associated with the input matrices. The total LOSC energy correction involves
 * the contributions from both alpha and beta spins. Thus, you may need to call
 * this function twice (for unrestricted KS calculation definitely).
 */
double energy_correction(ConstRefMat &Curvature, ConstRefMat &LocalOcc);

/**
 * @class __param__orbital_energy_post_scf
   @param [in] H_dfa The DFA Hamiltonian under AOs with dimension of
 * `[nbasis, nbasis]`.
 * @param [in] H_losc The LOSC effective Hamiltonian under AOs with dimension
 * of `[nbasis, nbasis]`.
 * @param [in] C_co The coefficient matrix of converged DFA COs under AOs
 * with dimension of `[nbasis, n]` (`n <= nbasis`, which is the number of COs).
 * You can choose whatever number of COs you want.
 */

/**
 * @class __return__orbital_energy_post_scf
 * @return A vector of the corrected orbital energies from LOSC with size of
 * `n`. The order of orbital energies matches the order of the input COs
 * coefficient matrix (the order of columns in `C_co` matrix).
 */

/**
 * @class __note__orbital_energy_post_scf
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

/**
 * @brief Calculate corrected orbital energy from LOSC in a post-SCF approach.
 * @details This function gives the FINAL orbital energies. The corrections from
 * LOSC are applied at the return. This is the difference to the
 * losc::energy_correction().
 *
 * The corrected orbital energies are defined as the expectation values of
 * converged DFA COs on the LOSC-DFA Hamiltonian, that is,
 * \f$
 * \epsilon_i = \langle \psi_i | H_{\rm{dfa}} + H_{\rm{losc}} | \psi_i \rangle.
 * \f$
 * @copydoc __param__orbital_energy_post_scf
 * @copydoc __return__orbital_energy_post_scf
 * @see losc::ao_hamiltonian_correction(): return the LOSC effective
 * Hamiltonian under AOs.
 * @copydoc __note__orbital_energy_post_scf
 */
vector<double> orbital_energy_post_scf(ConstRefMat &H_dfa, ConstRefMat &H_losc,
                                       ConstRefMat &C_co);

/**
 * @brief The C interface to calculate LOSC effective Hamiltonian under AOs.
 * @copydoc __param__ao_hamiltonian_correction
 * @param [in, out] H_losc a matrix with dimension of [nbasis, nbasis] at the
 * input. At exit, it is updated to store the LOSC effective Hamiltonian.
 * @copydoc __see__ao_hamiltonian_correction
 * @see losc::ao_hamiltonian_correction()
 */
void C_API_ao_hamiltonian_correction(ConstRefMat &S, ConstRefMat &C_lo,
                                     ConstRefMat &Curvature,
                                     ConstRefMat &LocalOcc, RefMat H_losc);

/**
 * @brief C interface to calculate the corrected orbital energies from LOSC in a
 * post-SCF approach.
 * @copydoc __param__orbital_energy_post_scf
 * @param [in, out] eig a array with size of at least `n` (the number of input
 * COs). At exit, the first `n` number of elements of `eig` are updated to be
 * the LOSC-DFA orbital energies.
 * @see losc::orbital_energy_post_scf()
 * @copydoc __note__orbital_energy_post_scf
 */
void C_API_orbital_energy_post_scf(ConstRefMat &H_dfa, ConstRefMat &H_losc,
                                   ConstRefMat &C_co, double *eig);

} // namespace losc

#endif
