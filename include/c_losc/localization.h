/**
 * @file
 * @brief The C interface for LOSC localization.
 */

#ifndef _LOSC_INCLUDE_C_LOSC_LOCALIZATION_H_
#define _LOSC_INCLUDE_C_LOSC_LOCALIZATION_H_

#include <c_losc/matrix.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

//**********************************************
// ==> Interface for `losc::LocalizerBase`
//**********************************************

typedef struct _LocalizerBase _LocalizerBase;

/**
 * @brief The C interface for LOSC localizer base.
 * @note This struct exports all the functions for LOSC localization base class.
 * This struct servers as helper class to avoid code duplication
 * for the implementation of C interface. There is no need for C users to
 * create a such struct in practice.
 */
typedef struct LoscLocalizerBase {
    _LocalizerBase *_p_base;

    // List interface of methods exposed to C users bellow.

    /**
     * @brief Set the maximum iteration number for localization.
     * @param [in] max_iter The maximum number of iteration.
     * @par Example
     * @code
     * // Assume `pbase` is a existing variable with type `LoscLocalizerBase *`.
     * // Demo: set max iteration number to be 100.
     * pbase->set_max_iter(pbase, 100);
     * @endcode
     */
    void (*set_max_iter)(const LoscLocalizerBase *self, size_t max_iter);

    /**
     * @brief Set the convergence of localization.
     * @param [in] tol The localization convergence tolerance.
     * @par Example
     * @code
     * // Assume `pbase` is a existing variable with type `LoscLocalizerBase *`.
     * // Demo: set the tolerance to be 1e-6.
     * pbase->set_convergence(pbase, 1.0e-6);
     * @endcode
     */
    void (*set_convergence)(const LoscLocalizerBase *self, double tol);

    /**
     * Set flag for doing random permutation or not for Jacobi-Sweep algorithm.
     * @param [in] flag `True` for doing random permutation `False` otherwise.
     * @par Example
     * @code
     * // Assume `pbase` is a existing variable with type `LoscLocalizerBase *`.
     * // Demo: set to false.
     * pbase->set_random_permutation(pbase, false);
     * @endcode
     */
    void (*set_random_permutation)(const LoscLocalizerBase *self, bool flag);

    /**
     * @brief Calculate the LOs and the unitary
     * @param [in, out] L A matrix with dimension of `[nbasis, nlo]`.
     * The data of `L` is ignored at entrace. At exit, it stores the LO
     * coefficient matrix.
     * @param [in, out] U A matrix with dimension of `[nlo, nlo]`.
     * At entrance, it stores a unitary matrix as the initial guess for the
     * localization. The unitarity of the input U matrix is not
     * verified. At exit, it is updated by the localization process.
     * @note The dimensions of input matrices will be check. It will throw an
     * exception if they do not match the expectation.
     * @note The relation between the LOs and AOs is
     * \f$ \displaystyle
     * \psi_i = \sum_\mu C_{\mu i} \phi_{\mu},
     * \f$
     * in which \f$ \psi_i \f$ is the LO, \f$ C_{\mu i} \f$ is the LO
     * coefficient matrix and \f$ \phi_i \f$ is the AO.
     * The dimension of LO coefficient matrix is `[nbasis, nlo]`.
     * The relation between the LOs and LO basis is
     * \f$ \displaystyle
     * \psi_i = \sum_j U_{\mu i} \phi_\mu,
     * \f$
     * in which \f$ \psi_i \f$ is the LO, \f$ U_{\mu i} \f$ is the U matrix
     * and \f$ \phi_\mu \f$ is the LO basis.
     * The dimension of the U matrix is `[nlo, nlo]`.
     */
    void (*lo_U)(const LoscLocalizerBase *self, losc_matrix *L, losc_matrix *U);
} LoscLocalizerBase;

//**********************************************
// ==> Interface for `losc::LoscLocalizerV2`
//**********************************************

/**
 * @class __code__demo_using_localizer_pbase
 * @code
 * // Demo: Set maximum iteration number to be 100.
 * pLocalizer->p_base->set_max_iter(pLocalizer, 100);
 *
 * // Demo: calculate LOs and U matrix.
 * // allocate LO coefficient matrix and U guess.
 * losc_matrix *L = losc_matrix_create(nbasis, nlo);
 * // assume U_data stores a unitary matrix.
 * losc_matrix *U = losc_matrix_create(nlo, nlo, U_data);
 * // calculate LO coefficient matrix and update U matrix.
 * pLocalizer->p_base->lo_U(pLocalizer->p_base, L, U);
 * @endcode
 * @see LoscLocalizerBase: see all the functions provided in the localizer
 * base structure.
 */

typedef struct _LoscLocalizerV2 _LoscLocalizerV2;

/**
 * @brief The C interface for LOSC localizer version 2.
 * @details LOSC localization version 2 is defined in Eq. 7 in
 * the LOSC2 paper (J. Phys. Chem. Lett. 2020, 11, 4, 1528–1535).
 * @see losc_localizer_v2_create(): create a LoscLocalizerV2.
 * @see losc_localizer_v2_create(): free a LoscLocalizerV2.
 */
typedef struct LoscLocalizerV2 {
    /**
     * This variable points to a real C++ `losc::LoscLocalizerV2` class object.
     * @note For Internal usage only.
     */
    _LoscLocalizerV2 *_p_v2;

    /**
     * @brief A pointer to a LoscLocalizerBase variable. Use this variable
     * to call functions that are provided in the LoscLocalizerBase.
     * @par Example
     * Demo of calling functions that provided in the LoscLocalizerBase.
     * Assume a variable `pLocalizer` exists with type of `LoscLocalizerV1 *`.
     * @copydoc __code__demo_using_localizer_pbase
     */
    LoscLocalizerBase *p_base;

    // @note
    // If `losc::LoscLocalizerV2` has new public functions (new functions added
    // in `losc::LoscLocalizerV2` not in its base class) in the future, add the
    // corresponding C interface with new function pointers below.
} LoscLocalizerV2;

/**
 * @brief Constructor of LoscLoscLocalizerV2.
 * @param [in] C_lo_basis LO basis coefficient matrix under AO with
 * dimension `[nbasis, nlo]`.
 * @param [in] H_ao Hamiltonian matrix under AO used in localization with
 * dimension `[nbasis, nbasis]`.
 * @param [in] D_ao dipole matrix under AO in x, y and z order with
 * dimension `[nbasis, nbasis]`.
 * @return create a LoscLoscLocalizerV1 struct and return the pointer.
 */
LoscLocalizerV2 *losc_localizer_v2_create(
    const losc_matrix *C_lo_basis, const losc_matrix *H_ao,
    const losc_matrix *D_ao[3]);

/**
 * Free a LoscLoscLocalizerV2.
 * @param [in, out] ptr a pointer with type `LoscLocalizerV2 *`. At exit,
 * `ptr` is set to null.
 */
#define losc_localizer_v2_free(ptr) _losc_localizer_v2_free(&(ptr))
void *_losc_localizer_v2_free(LoscLocalizerV2 **pptr_self);

#ifdef __cplusplus
}
#endif

#endif
