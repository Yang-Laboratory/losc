#ifndef _LOSC_CORRECTION_H_
#define _LOSC_CORRECTION_H_

#include <matrix/matrix.h>
#include <memory>
#include <vector>

namespace losc {

using matrix::Matrix;
using std::vector;
using SharedMatrix = std::shared_ptr<Matrix>;

/**
 * Calculate losc correcting Hamiltonian under AO basis.
 *
 * The Losc correcting Hamiltonian is constructed with LO frozen. The formula is
 * expressed as Eq. S25 in the supporting information of the original Losc paper
 * (https://doi.org/10.1093/nsr/nwx111).
 *
 * @ param [in] S: AO overlap matrix with dimension [nbasis, nbasis].
 * @ param [in] C_lo: LO coefficient matrix under AO basis with dimension [nlo, nbasis].
 * @ param [in] Curvature: Losc curvature matrix with dimension [nlo, nlo].
 * @ param [in] LocalOcc: Losc local occupation matrix with dimension [nlo, nlo].
 * @ return: the Losc correcting Hamiltonian with dimension [nbasis, nbasis].
*
 * Note:
 * 1. Make sure all the matrices have the same spin, if that is the case.
 */
SharedMatrix losc_hamiltonian_correction(const Matrix &S, const Matrix &C_lo,
                                         const Matrix &Curvature, const Matrix &LocalOcc);

/**
 * Calculate Losc correction to total energy.
 *
 * @ param [in] Curvature: Losc curvature matrix with dimension [nlo, nlo].
 * @ param [in] LocalOcc: Losc local occupation matrix with dimension [nlo, nlo].
 * @ return: the Losc correction to total energy.
 *
 * Note:
 * 1. Make sure all the input matrices have the same spin.
 * 2. The total Losc correction to the total energy is the summation of both alpha and beta spin.
 */
double losc_total_energy_correction(const Matrix &Curvature, const Matrix &LocalOcc);

/**
 * Calculate Losc correction to orbital energy with direct correction.
 *
 * The dirrect correction $\Delta \epsilon$ is defined in Eq. 11 in the original
 * Losc paper (https://doi.org/10.1093/nsr/nwx111).
 *
 * @ param [in] S: AO overlap matrix with dimension [nbasis, nbasis].
 * @ param [in] C_co: CO coefficient matrix under AO basis with dimension [nlo, nbasis].
 * @ param [in] C_lo: LO coefficient matrix under AO basis with dimension [nlo, nbasis].
 * @ param [in] Curvature: Losc curvature matrix with dimension [nlo, nlo].
 * @ param [in] LocalOcc: Losc local occupation matrix with dimension [nlo, nlo].
 * @ return: the Losc correction to orbital energy.
 *
 * Note:
 * 1. Make sure all the input matrices have the same spin, if that is the case.
 * 2. If the `nlo` is not equal the `nbasis` (in the case of selecting window for
 *    orbitals in the Losc localization process), the input `C_co` is NO longer a
 *    square matrix. So make sure you slice the square CO coefficient matrix before
 *    feed it in this function.
 * 3. The size of returned vector is the number of the LOs. The order of the correction
 *    is aligned to the COs' order in `C_co` matrix.
 */
vector<double> losc_orbital_energy_correction(const Matrix &S, const Matrix &C_co,
                                              const Matrix &C_lo, const Matrix &Curvature,
                                              const Matrix &LocalOcc);

/**
 * Calculate Losc corrected orbital energies with projection to CO coefficient.
 *
 * @ param [in] H_dfa: DFA Hamiltonian under AO with dimension [nbasis, nbasis].
 * @ param [in] H_losc: Losc correcting Hamiltonian under AO with dimension [nbasis, nbasis].
 * @ param [in] C_co: CO coefficient matrix under AO basis with dimension [n, nbasis].
 * @ return: the Losc corrected orbital energies for the `n` COs.
 *
 * Note:
 * 1. Make sure all the input matrices have the same spin.
 * 2. The return orbital energies have already been corrected by Losc. No further steps to take.
 * 3. The dimension of `C_co` is [n, nbasis], where n <= nbasis. Usually, you will will choose
 *    n=nlo for all the localized orbitals, but this is required.
 */
vector<double> losc_corrected_orbital_energy_by_projection(const Matrix &H_dfa, const Matrix &H_losc,
                                                           const Matrix &C_co);

/**
 * Calculate Losc corrected orbital energies by diagonalizing the Losc corrected total Hamiltonian matrix.
 *
 * @ param [in] H_dfa: DFA Hamiltonian under AO with dimension [nbasis, nbasis].
 * @ param [in] H_losc: Losc correcting Hamiltonian under AO with dimension [nbasis, nbasis].
 * @ param [in] S_neg_half: AO overlap matrix to the power of (-1/2), that is S^(-1/2).
 * @ return: the Losc corrected orbital energies for the all the COs.
 *
 * Note:
 * 1. Make sure all the input matrices have the same spin.
 * 2. The return orbital energies have already been corrected by Losc. No further steps to take.
 */
vector<double> losc_corrected_orbital_energy_by_diagonalize(const Matrix &H_dfa, const Matrix &H_losc,
                                                            const Matrix &S_neg_half);

}

#endif // _LOSC_CORRECTION_H_
