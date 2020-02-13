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
 * calculate losc correction to Hamiltonian.
 *
 * H_eff = S * C^T  * A * C * S:
 * S: AO overlap matrix.
 * C: LO coef matrix.
 * A: A_ij = \delta_{ij} 1/2 K_ii - K_ij \lambda_ij, K is kappa, \lambda is local occ.
 */
SharedMatrix losc_hamiltonian_matrix(const Matrix &S, const Matrix &C_lo, const Matrix &Curvature,
                                     const Matrix &LocalOcc);

/**
 * calculate losc correction to total energy.
 */
double losc_energy(const Matrix &Curvature, const Matrix &LocalOcc);

/**
 * calculate losc correction to orbital energy with direct correction.
 */
vector<double> losc_orbital_energy_direct_correction(const Matrix &S, const Matrix &C_co,
                                                     const Matrix &C_lo, const Matrix &Curvature,
                                                     const Matrix &LocalOcc);

/**
 * calculate orbital energy by projection.
 */
vector<double> losc_orbital_energy_by_projection(const Matrix &H_dfa, const Matrix &H_losc,
                                                 const Matrix &C_co);

/**
 * calculate orbital energy by diagnoalization.
 */
vector<double> losc_orbital_energy_by_diagonalize(const Matrix &Shalf, const Matrix &H_dfa,
                                                  const Matrix &H_losc);
}


#endif // _LOSC_CORRECTION_H_
