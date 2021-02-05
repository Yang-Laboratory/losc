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
 * @details This struct exports the functions in `losc::LocalizerBase` in
 * the C++ code. This struct servers as helper class to avoid code duplication
 * for the implementation of C interface. C users have no access to create
 * a such struct in practice, because we do not provide a constructor function.
 */
typedef struct LoscLocalizerBase {
    _LocalizerBase *_p_base;

    // List interface of methods exposed to C users bellow.

    /**
     * Set the maximum iteration number for localization.
     */
    void (*set_max_iter)(const LoscLocalizerBase *self, size_t max_iter);

    /**
     * Set the convergence of localization.
     */
    void (*set_convergence)(const LoscLocalizerBase *self, double tol);

    /**
     * Set flag for doing random permutation or not for Jacobi-Sweep algorithm.
     */
    void (*set_random_permutation)(const LoscLocalizerBase *self, bool flag);

    /**
     * @brief Calculate the LOs and the unitary
     * @param L [in, out]: Dimension of [nbasis, nlo]. The content of L is
     * ignored. At exit, it is overwrote to be the LO coefficient matrix.
     * @param U [in, out]: Dimension of [nlo, nlo]. At entrance, it stores a
     * unitary matrix as the initial guess for the localization. The unitarity
     * of the input U matrix is not verified. At exit, it is updated by the
     * localization process.
     *
     * @note
     * 1. The dimensions of input matrices will be check. It will throw an
     * exception if they do not match the expectation.
     */
    void (*lo_U)(const LoscLocalizerBase *self, losc_matrix *L, losc_matrix *U);
} LoscLocalizerBase;

//**********************************************
// ==> Interface for `losc::LoscLocalizerV2`
//**********************************************

typedef struct _LoscLocalizerV2 _LoscLocalizerV2;

/**
 * @brief The C interface for LOSC localizer version 2.
 * @details LOSC localization version 2 is defined in Eq. 7 in
 * the LOSC2 paper (J. Phys. Chem. Lett. 2020, 11, 4, 1528–1535).
 */
typedef struct LoscLocalizerV2 {
    /**
     * This variable points to a real C++ `losc::LoscLocalizerV2` class object.
     *
     * @note
     * For Internal usage only.
     */
    _LoscLocalizerV2 *_p_v2;

    /**
     * @brief This variable stores a bunch of function pointers that point to
     * `losc::LocalizerBase` public interface.
     *
     * @example
     * How to call the functions to do localization:
     * @code
     * // First create pV1 that points to a LoscLoscLocalizerV1 struct by using
     * // losc_localizer_v1_create function.
     * // Set maximum iteration number to be 100.
     * pV1->p_base->set_max_iter(pV1, 100);
     * // allocate LO coefficient matrix and U guess.
     * losc_matrix *L = losc_matrix_create(nbasis, nlo);
     * // assume U_data stores a unitary matrix.
     * losc_matrix *U = losc_matrix_create(nlo, nlo, U_data);
     * // calculate LO coefficient matrix and update U matrix.
     * pV1->p_base(pV1, L, U);
     * @endcode
     *
     * @see
     * LoscLocalizerBase: For how to all the supported functions in the base.
     *
     * @example
     * How to use base functions.
     */
    LoscLocalizerBase *p_base;

    /**
     * @note
     * If `losc::LoscLocalizerV2` has new public functions (new functions added
     * in `losc::LoscLocalizerV2` not in its base class) in the future, add the
     * corresponding C interface with new function pointers below.
     */
} LoscLocalizerV2;

/**
 * Constructor of LoscLoscLocalizerV2.
 *
 * @param [in] C_lo_basis: LO basis coefficient matrix under AO with
 * dimension [nbasis, nlo].
 * @param [in] H_ao: Hamiltonian matrix under AO used in localization with
 * dimension [nbasis, nbasis].
 * @param [in] Dipole_ao: dipole matrix under AO in x, y and z order with
 * dimension [nbasis, nbasis].
 *
 * @return a pointer to LoscLoscLocalizerV1 struct.
 */
LoscLocalizerV2 *losc_localizer_v2_create(
    const losc_matrix *C_lo_basis, const losc_matrix *H_ao,
    const losc_matrix *D_ao[3]);

/**
 * Destructor of LoscLoscLocalizerV2.
 */
void *_losc_localizer_v2_free(LoscLocalizerV2 **pptr_self);
#define losc_localizer_v2_free(ptr_self) _losc_localizer_v2_free(&(ptr_self))

#ifdef __cplusplus
}
#endif

#endif
