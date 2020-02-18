/**
 * @file
 * @brief declaration relates to Losc correction.
 */
#ifndef _LOSC_SRC_CORRECTION_H_
#define _LOSC_SRC_CORRECTION_H_

#include <matrix/matrix.h>
#include <memory>
#include <vector>

/**
 * @brief top-level losc namespace.
 */
namespace losc {

using matrix::Matrix;
using std::vector;
using SharedMatrix = std::shared_ptr<Matrix>;

/**
 * @brief Calculate Losc correcting Hamiltonian under AO basis.
 *
 * @details The Losc correcting Hamiltonian is constructed with LO frozen. The
 * formula is expressed as Eq. S25 in the supporting information of the
 * original Losc paper (https://doi.org/10.1093/nsr/nwx111).
 *
 * @param [in] S: AO overlap matrix with dimension [nbasis, nbasis].
 * @param [in] C_lo: LO coefficient matrix under AO basis with dimension
 * [nlo, nbasis]. see losc::LocalizerBase::compute().
 * @param [in] Curvature: Losc curvature matrix with dimension [nlo, nlo].
 * see losc::CurvatureBase::compute().
 * @param [in] LocalOcc: Losc local occupation matrix with dimension [nlo, nlo].
 * see losc::local_occupation_matrix().
 * @return SharedMatrix: the Losc correcting Hamiltonian with dimension
 * [nbasis, nbasis].
 *
 * @note Make sure all the input matrices have the same spin. The returned
 * Losc correcting Hamiltonian matrix is only for the input spin.
 */
SharedMatrix losc_hamiltonian_correction(const Matrix &S, const Matrix &C_lo,
                                         const Matrix &Curvature,
                                         const Matrix &LocalOcc);

/**
 * @brief Calculate Losc correction to total energy.
 *
 * @param [in] Curvature: Losc curvature matrix with dimension [nlo, nlo].
 * see losc::CurvatureBase::compute().
 * @param [in] LocalOcc: Losc local occupation matrix with dimension [nlo, nlo].
 * see losc::local_occupation_matrix().
 * @return double: the Losc correction to total energy.
 *
 * @note Make sure all the input matrices have the same spin. The return total
 * energy correction is only for the input spin.
 * @note The total Losc correction to the total energy is the summation of both
 * alpha and beta spin.
 */
double losc_total_energy_correction(const Matrix &Curvature,
                                    const Matrix &LocalOcc);

/**
 * @brief Calculate Losc correction to orbital energy with direct correction.
 *
 * @details The dirrect correction $\Delta \epsilon$ is defined in Eq. 11 in the
 * original Losc paper (https://doi.org/10.1093/nsr/nwx111).
 *
 * @param [in] S: AO overlap matrix with dimension [nbasis, nbasis].
 * @param [in] C_co: CO coefficient matrix under AO basis with dimension [nlo,
 * nbasis].
 * @param [in] C_lo: LO coefficient matrix under AO basis with dimension [nlo,
 * nbasis]. see losc::LocalizerBase::compute().
 * @param [in] Curvature: Losc curvature matrix with dimension [nlo, nlo].
 * see losc::CurvatureBase::compute().
 * @param [in] LocalOcc: Losc local occupation matrix with dimension [nlo, nlo].
 * see losc::local_occupation_matrix().
 * @return std::vector<double>: the Losc correction to orbital energy.
 *
 * @note Make sure all the input matrices have the same spin. The returned
 * orbital energy correction is only for input spin.
 * @note If you select only parts of orbitals (not all the orbitals, that is
 * `nlo` is less than `nbasis`) to do localization, the input `C_co` should NO
 * longer be a square matrix. In this case, you need to slice the original
 * squared CO coefficient matrix, and make sure the used \p C_co matrix
 * corresponds to all the COs you used in the localization process.
 * @note The size of returned vector is the number of the LOs. The order of the
 * correction is aligned to the COs' order in `C_co` matrix.
 */
vector<double> losc_orbital_energy_correction(const Matrix &S,
                                              const Matrix &C_co,
                                              const Matrix &C_lo,
                                              const Matrix &Curvature,
                                              const Matrix &LocalOcc);

/**
 * @brief Calculate Losc corrected orbital energies with projection to CO
 * coefficient.
 *
 * @param [in] H_dfa: DFA Hamiltonian under AO with dimension [nbasis, nbasis].
 * @param [in] H_losc: Losc correcting Hamiltonian under AO with dimension
 * [nbasis, nbasis]. See losc::losc_hamiltonian_correction().
 * @param [in] C_co: CO coefficient matrix under AO basis with dimension [n,
 * nbasis].
 * @return std::vector<double>: the Losc corrected orbital energies for the `n`
 * COs.
 *
 * @note Make sure all the input matrices have the same spin. The returned
 * orbital corrections is only for the input spin.
 * @note The return orbital energies have already been corrected by Losc. No
 * further steps to take.
 * @note The dimension of `C_co` is [n, nbasis], where n <= nbasis. Usually, you
 * will will choose n=nlo for all the localized orbitals, but this is required.
 */
vector<double> losc_corrected_orbital_energy_by_projection(const Matrix &H_dfa,
                                                           const Matrix &H_losc,
                                                           const Matrix &C_co);

/**
 * @brief Calculate Losc corrected orbital energies by diagonalizing the Losc
 * corrected total Hamiltonian matrix.
 *
 * @param [in] H_dfa: DFA Hamiltonian under AO with dimension [nbasis, nbasis].
 * @param [in] H_losc: Losc correcting Hamiltonian under AO with dimension
 * [nbasis, nbasis]. See losc::losc_hamiltonian_correction().
 * @param [in] S_neg_half: AO overlap matrix to the power of (-1/2), that is
 * S^(-1/2).
 * @return std::vector<double>: the Losc corrected orbital energies for the all
 * the COs.
 *
 * @note Make sure all the input matrices have the same spin. The return orbital
 * corrections are only for the input spin.
 * @note The return orbital energies have already been corrected by Losc. No
 * further steps to take.
 */
vector<double> losc_corrected_orbital_energy_by_diagonalize(
    const Matrix &H_dfa, const Matrix &H_losc, const Matrix &S_neg_half);

} // namespace losc

#endif // _LOSC_SRC_CORRECTION_H_
