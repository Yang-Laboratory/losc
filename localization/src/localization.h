/**
 * file: localization.h
 */
#ifndef _LOCALIZATION_LOCALIZATION_H_
#define _LOCALIZATION_LOCALIZATION_H_

#include <matrix/matrix.h>
#include <vector>
#include <iostream>
#include <cstddef>

namespace localization {

using matrix::Matrix;
using std::vector;

/**
 * Losc2 localization data interface.
 *
 * Prepare all the data specified in this struct before
 * calling Losc2 localization function, localization::LocalizerLosc2(),
 * to create the localization orbital coefficient matrix.
 *
 * Note: All the matrices are treated as general matrix.
 */
struct LocalizerLosc2InputData
{
    /**
     * The number of localized orbitals.
     * Default is 0.
     */
    size_t nlmo = 0;

    /**
     * The number of atomic basis set.
     * Default is 0.
     */
    size_t nbasis = 0;

    /**
     * If or not to use the random permutation in Jacobi-Sweep
     * algorithm for localization. Default is true.
     */
    bool use_js_random_permutation = true;

    /**
     * The maximum iteration number for Jacobi-Sweep
     * algorithm. Default is 1000.
     */
    size_t js_max_iter = 1000;

    /**
     * Jacobi-Sweep localization convergence tolerance.
     * Default is 1e-10.
     */
    double js_tol = 1e-10;

    /**
     * Losc2 localization parameter C.
     * Default is 1000, that is the default value for Losc2.
     */
    double para_c = 1000;

    /**
     * Losc2 localization parameter gamma.
     * Default is 0.78, that is the optimized value for Losc2.
     */
    double para_gamma = 0.78;

    /**
     * u_matrix points to a unitary matrix, dimension (nlmo, nlmo), that
     * transfers the localized orbital (LO) basis to localized orbital.
     *
     * Each row of the matrix refers to one LO that expanded on the LO basis.
     * Default u_matrix is a null pointer.
     */
    Matrix *u_matrix = nullptr;

    /**
     * lo_basis_coef points to a matrix, dimension (nlmo, nbasis), that refers
     * to the LO basis coefficient matrix expanded on atomic basis set.
     *
     * Each row of the matrix refers to one LO basis that expanded on the atomic basis.
     * Default is a null pointer.
     */
    const Matrix *lo_basis_coef = nullptr;

    /**
     * hamiltonian_ao points to a matrix, dimension (nbasis, nbasis),
     * that is the Hamiltonian matrix under atomic basis used in the Losc2 localization.
     *
     * Default is a null pointer.
     */
    const Matrix *hamiltonian_ao = nullptr;

    /**
     * dipole_ao is a vector of matrix pointer.
     * The x, y and z component of the dipole matrix pointer are stored in order in the vector.
     * The dimension of each dipole component matrix is (nbasis, nbasis).
     *
     * Default all the dipole matrix pointers are null pointer.
     */
    vector<const Matrix *> dipole_ao = {nullptr, nullptr, nullptr};
};

/**
 * Create a localized orbital coefficient matrix that expanded on atomic basis set
 * based on the input localization data.
 *
 * @ param[in] input_data: the input data for localization. On exit, the input_data.u_matrix
 *  is updated with final U matrix.
 * @ return Matrix: the LO coefficient matrix on atomic basis.
 */
Matrix *LocalizerLosc2(LocalizerLosc2InputData &input_data);


}

#endif
