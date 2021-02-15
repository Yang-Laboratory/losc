/**
 * @file correction.hpp
 * @brief C++ interface for the correction from LOSC.
 */
#ifndef _LOSC_INCLUDE_LOSC_CORRECTION_HPP_
#define _LOSC_INCLUDE_LOSC_CORRECTION_HPP_

#include <losc/eigen_def.hpp>
#include <vector>

namespace losc {

using std::vector;
/**
 * @brief Calculate LOSC effective Hamiltonian under AO basis.
 *
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
 *
 * @param [in] S: AOs overlap matrix with dimension of [nbasis, nbasis].
 * @param [in] C_lo: LOs coefficient matrix under AOs with dimension of
 * [nbasis, nlo].
 * @param [in] Curvature: LOSC curvature matrix with dimension of [nlo, nlo].
 * @param [in] LocalOcc: LOSC local occupation matrix with dimension of
 * [nlo, nlo].
 *
 * @return LOSCMatrix: the LOSC effective Hamiltonian with dimension of
 * [nbasis, nbasis].
 *
 * @see
 * losc::LoscLocalizer::lo(): return the LOs coefficient matrix.
 * losc::CurvatureBase::kappa(): return the LOSC curvature matrix.
 * losc::local_occupation(): return the LOSC local occupation matrix.
 */
LOSCMatrix ao_hamiltonian_correction(ConstRefMat &S, ConstRefMat &C_lo,
                                     ConstRefMat &Curvature,
                                     ConstRefMat &LocalOcc);

/**
 * @brief The C interface to calculate LOSC effective Hamiltonian under AOs.
 *
 * @param H_losc [in, out]: a matrix with dimension of [nbasis, nbasis] at the
 * input. At exit, it is updated to store the LOSC effective Hamiltonian.
 *
 * @see
 * losc::ao_hamiltonian_correction(): requirements for other paramters.
 */
void C_API_ao_hamiltonian_correction(ConstRefMat &S, ConstRefMat &C_lo,
                                     ConstRefMat &Curvature,
                                     ConstRefMat &LocalOcc, RefMat H_losc);

/**
 * @brief Calculate the total energy correction from LOSC.
 *
 * This is just the energy correction from LOSC, NOT the total energy of
 * LOSC-DFA. Total energy of LOSC-DFA is: E_losc_dfa = E_dfa + E_losc.
 *
 * @param [in] Curvature: LOSC curvature matrix with dimension of [nlo, nlo].
 * @param [in] LocalOcc: LOSC local occupation matrix with dimension of
 * [nlo, nlo].
 *
 * @return double: the total energy correction from LOSC.
 *
 * @see
 * losc::CurvatureBase::kappa(): return the LOSC curvature matrix.
 * losc::local_occupation(): return the LOSC local occupation matrix.
 *
 * @note
 * This function only calculates the energy correction of the single spin
 * associated with the input matrices. The total LOSC energy correction involves
 * the contributions from both alpha and beta spins. Thus, you may need to call
 * this function twice (for unrestricted KS calculation definitely).
 */
double energy_correction(ConstRefMat &Curvature, ConstRefMat &LocalOcc);

/**
 * @brief Calculate corrected orbital energy from LOSC in a post-SCF approach.
 *
 * @details This function gives the FINAL orbital energies. The corrections from
 * LOSC are applied at the return. This is the difference to the
 * losc::energy_correction().
 *
 * The corrected orbital energies are the defined as the expectation values of
 * converged DFA COs on the LOSC-DFA Hamiltonian, that is,
 * \f[
 * \epsilon_i = \langle \psi_i | H_{\rm{dfa}} + H_{\rm{losc}} | \psi_i \rangle.
 * \f]
 *
 * @param [in] H_dfa: The DFA Hamiltonian under AOs with dimension of
 * [nbasis, nbasis].
 * @param [in] H_losc: The LOSC effective Hamiltonian under AOs with dimension
 * of [nbasis, nbasis].
 * @param [in] C_co: The coefficient matrix of converged DFA COs under AOs
 * with dimension of [nbasis, n] (n <= nbasis, which is the number of COs).
 * You can choose whatever number of COs you want.
 *
 * @return vector<double>: the corrected orbital energies from LOSC with size
 * of n. The order of orbital energies matches the order of the input COs
 * coefficient matrix (the order of columns in `C_co` matrix).
 *
 * @see
 * losc::ao_hamiltonian_correction(): return the LOSC effective Hamiltonian
 * under AOs.
 *
 * @note
 * This function only calculates the orbital energies for the single spin
 * associated with the input matrices.
 * @note
 * This function is just one of the ways to construct the LOSC corrected orbital
 * energies in the post-SCF LOSC calculations. It is the way we used to produce
 * results in the published paper for the post-SCF LOSC calculations. Besides
 * this way, there are another two ways to calculate corrected orbital energies:
 * (1) diagonalize the corrected LOSC-DFA Hamiltonian; (2) Follow Eq. 11 of
 * the original paper of LOSC to calculate the corrections to orbital energies.
 * These three ways usually produce very similar results. Here, we only provide
 * this approach approach and suggest to use it in practice.
 */
vector<double> orbital_energy_post_scf(ConstRefMat &H_dfa, ConstRefMat &H_losc,
                                       ConstRefMat &C_co);

/**
 * @brief C interface to calculate the corrected orbital energies from LOSC in a
 * post-SCF approach.
 *
 * @param eig [in, out]: a array with size of at least `n` (the number of input
 * COs). At exit, the first `n` number of elements of `eig` are updated to be
 * the LOSC-DFA orbital energies.
 *
 * @see
 * losc::orbital_energy_post_scf(): the requirements of other parameters.
 */
void C_API_orbital_energy_post_scf(ConstRefMat &H_dfa, ConstRefMat &H_losc,
                                   ConstRefMat &C_co, double *eig);

} // namespace losc

#endif
