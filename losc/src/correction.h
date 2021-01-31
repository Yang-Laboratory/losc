/**
 * @file correction.h
 * @brief C++ interface for the correction from LOSC.
 */
#ifndef __LOSC_SRC_CORRECTION_H__
#define __LOSC_SRC_CORRECTION_H__

#include "eigen_def.h"
#include <vector>

namespace losc {

using std::vector;
/**
 * @brief Calculate LOSC effective Hamiltonian under AO basis.
 *
 * @details The LOSC effective Hamiltonian is constructed with LOs fixed. The
 * expression of the effective Hamiltonian is shown as Eq. S25 in the
 * supporting information of the original LOSC paper
 * (https://doi.org/10.1093/nsr/nwx111). This effective Hamiltonian is exact
 * in the developed version of SCF-LOSC (self-consistent LOSC). See reference
 * (J. Phys. Chem. Lett. 2020, 11, 23, 10269-10277) for more details about
 * how to perform reliable SCF-LOSC calculations.
 *
 * @param [in] S: AO overlap matrix with dimension [nbasis, nbasis].
 * @param [in] C_lo: LO coefficient matrix under AO basis with dimension
 * [nbasis, nlo].
 * @param [in] Curvature: LOSC curvature matrix with dimension [nlo, nlo].
 * @param [in] LocalOcc: LOSC local occupation matrix with dimension [nlo, nlo].
 *
 * @return MatrixXd: the LOSC effective Hamiltonian with dimension
 * [nbasis, nbasis].
 *
 * @see
 * losc::CurvatureBase::kappa(): Obtain the curvature matrix.
 * losc::local_occupation(): Obtain the local occupation matrix.
 *
 * @note
 * Make sure all the input matrices are corrresponding to the same spin.
 */
MatrixXd ao_hamiltonian_correction(ConstRefMat &S, ConstRefMat &C_lo,
                                   ConstRefMat &Curvature,
                                   ConstRefMat &LocalOcc);

/**
 * @brief Calculate the total energy correction from LOSC.
 *
 * This is just the energy correction from LOSC, NOT the total energy of
 * LOSC-DFA. Total energy of LOSC-DFA is: E_losc_dfa = E_dfa + E_losc.
 *
 * @param [in] Curvature: LOSC curvature matrix with dimension [nlo, nlo].
 * @param [in] LocalOcc: LOSC local occupation matrix with dimension [nlo, nlo].
 *
 * @return double: the total energy correction from LOSC.
 *
 * @see
 * losc::CurvatureBase::kappa(): obtain the LOSC curvature matrix.
 * losc::local_occupation(): obtain the LOSC local occupation matrix.
 */
double energy_correction(ConstRefMat &Curvature, ConstRefMat &LocalOcc);

/**
 * @brief Calculate corrected orbital energy from LOSC in a post-SCF approach.
 *
 * @details This function gives the final orbital energies WITH the correction
 * from LOSC. Note the difference to the function energy_correction() that only
 * calculates the energy correction. The corrected orbital energies are the
 * expectation values of converged DFA's COs on the LOSC-DFA Hamiltonian,
 * that is,
 * \f[
 * \epsilon_i = \langle \psi_i | H_{\rm{dfa}} + H_{\rm{losc}} | \psi_i \rangle.
 * \f]
 *
 * @param [in] H_dfa: The DFA Hamiltonian under AOs with dimension of
 * [nbasis, nbasis].
 * @param [in] H_losc: The LOSC effective Hamiltonian under AOs with dimension
 * of [nbasis, nbasis].
 * @param [in] C_co: The coefficient matrix of converged DFA's COs under AOs
 * with dimension of [nbasis, n] (n <= nbasis, which is the number of COs).
 *
 * @return vector<double>: the corrected orbital energies from LOSC with size
 * of n. The order of orbital energies match the order of input COs (order
 * of columns in C_co matrix).
 *
 * @see
 * losc::ao_hamiltonian_correction(): obtain the LOSC effective Hamiltonian
 * under AOs.
 *
 * @note
 * This function is just one of the ways to construct the LOSC corrected orbital
 * energy in a post-SCF LOSC calculation. It is the way we used to produce
 * results in the published paper for the post-SCF LOSC calculations. Besides
 * this way, there are another two ways to calculate corrected orbital energies:
 * (1) diagonalize the corrected LOSC-DFA Hamiltonian; (2) Follow Eq. 11 to
 * calculate the corrections to orbital energies. These three ways usually
 * prouce very similar results. Here, we only provide the commonly used
 * approach.
 */
vector<double> orbital_energy_post_scf(ConstRefMat &H_dfa, ConstRefMat &H_losc,
                                       ConstRefMat &C_co);

} // namespace losc

#endif
