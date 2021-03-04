/**
 * @file
 * @brief The C interface for LOSC curvature matrix.
 */

#ifndef _LOSC_INCLUDE_C_LOSC_CURVATURE_H_
#define _LOSC_INCLUDE_C_LOSC_CURVATURE_H_

#include <c_losc/matrix.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

//**********************************************
// ==> Interface for `losc::DFAInfo`
//**********************************************

/**
 * @brief The DFA information.
 * @note This is a incomplete declaration. You can not create a LoscDFAInfo
 * variable by youself.
 * @see losc_dfa_info_create(): create a LoscDFAInfo
 * @see losc_dfa_info_free(): free a LoscDFAInfo
 */
typedef struct LoscDFAInfo LoscDFAInfo;

/**
 * @brief Create a struct LoscDFAInfo.
 * @param [in] gga_x: The total weights of all GGA and LDA type exchanges.
 * @param [in] hf_x: The total weights of HF exchanges.
 * @param [in] name: The name of the DFA. Default to an empty string.
 * @par Example
 * Taking B3LYP functional as an example. The B3LYP functional is
 * \f[
 *      E^{\rm{B3LYP}}_{\rm{xc}} = E^{\rm{LDA}}_{\rm{x}} +
 *      a_0 (E^{\rm{HF}}_{\rm{x}} - E^{\rm{LDA}}_{\rm{x}}) +
 *      a_x (E^{\rm{GGA}}_{\rm{x}} - E^{\rm{LDA}}_{\rm{x}}) +
 *      E^{\rm{LDA}}_{\rm{c}} + a_c (E^{\rm{GGA}}_{\rm{c}} -
 *      E^{\rm{LDA}}_{\rm{c}}),
 * \f]
 * in which exchanges end with suffix ``_x`` and correlations end with
 * suffix ``_c``, \f$ a_0=0.20 \f$, \f$ a_x=0.72 \f$ and \f$ a_c=0.81 \f$.
 * Only the exchanges are considered and the correlations are ignored.
 * The GGA and LDA exchanges are viewed as the same type. Therefore,
 * the total weights of GGA and LDA exchanges are
 * \f$ 1 - a_0 + a_x \times (1 - 1) = 1 - a_0 = 0.80 \f$,
 * and the total weights of HF exchanges is
 * \f$ a_0 = 0.20 \f$. To construct a ``DFAInfo`` object for B3LYP
 * functional, one should do the following
 * @code
 * b3lyp = losc_dfa_info_create(0.80, 0.20, "B3LYP")
 * @endcode
 */
LoscDFAInfo *losc_dfa_info_create(double gga_x, double hf_x, const char *name);

/**
 * @brief Free a LoscDFAInfo struct.
 * @param [in, out] ptr a pointer with type `LoscDFAInfo *`. At exit,
 * `ptr` is set to null.
 */
#define losc_dfa_info_free(ptr) _losc_DFAInfo_free(&(ptr))
void _losc_dfa_info_free(LoscDFAInfo **pptr_m);

//**********************************************
// ==> Interface for `losc::CurvatureBase`
//**********************************************

// The incomplete struct in C-side that is `losc::CurvatureBase` in C++ side.
typedef struct _LoscCurvatureBase _LoscCurvatureBase;

/**
 * @brief The C interface for LOSC curvature base class.
 * @note This struct exports all the functions for LOSC curvature base class.
 * This struct servers as helper class to avoid code duplication
 * for the implementation of C interface. There is no need for C users to
 * create a such struct in practice.
 */
typedef struct LoscCurvatureBase {
    // points to a real `losc::CurvatureBase` object.
    _LoscCurvatureBase *_p_base;

    // List interface of methods exposed to C users bellow.

    /**
     * @brief Return the number of LOs.
     * @par Example
     * @code
     * // Assume `pbase` is a existing variable with type `LoscCurvatureBase *`.
     * // Demo: get the number of LOs.
     * size_t n = pbase->nlo(pbase);
     * @endcode
     */
    size_t (*nlo)(const LoscCurvatureBase *self);

    /**
     * @brief Return the number of nfitbasis.
     * @par Example
     * @code
     * // Assume `pbase` is a existing variable with type `LoscCurvatureBase *`.
     * // Demo: get the number of fit basis.
     * size_t n = pbase->nfitbasis(pbase);
     * @endcode
     */
    size_t (*nfitbasis)(const LoscCurvatureBase *self);

    /**
     * @brief Return the number of grid points.
     * @par Example
     * @code
     * // Assume `pbase` is a existing variable with type `LoscCurvatureBase *`.
     * // Demo: get the number of grid points.
     * size_t n = pbase->npts(pbase);
     * @endcode
     */
    size_t (*npts)(const LoscCurvatureBase *self);

    /**
     * @brief Compute the LOSC curvature matrix with an in-out way.
     * @param [in, out] K curvature matrix with dimension `[nlo, nlo]`.
     * @note The curvature matrix should be allocated in advance with an
     * expected dimension. Otherwise, it will throw an exception.
     * @par Example
     * @code
     * // Assume `pbase` is a existing variable with type `LoscCurvatureBase *`.
     * // Demo: calculate curvature matrix.
     * nlo = pbase->nlo(pbase);
     * losc_matrix *K = losc_create_matrix(nlo, nlo);
     * pbase->kappa(pbase, K);
     * @endcode
     */
    void (*kappa)(const LoscCurvatureBase *self, losc_matrix *K);
} LoscCurvatureBase;

//**********************************************
// ==> Interface for `losc::CurvatureV1`
//**********************************************

// The incomplete struct in C-side that is `losc::CurvatureV1` in C++-side.
typedef struct _LoscCurvatureV1 _LoscCurvatureV1;

/**
 * @class __code__demo_using_curvature_pbase
 * @code
 * // Demo: get number of AOs by calling LoscCurvatureV1.nbasis
 * size_t nbasis = pCurvature->p_base->nbasis(pCurvature);
 *
 * // Demo: get number of LOs by calling LoscCurvatureV1.nlo
 * size_t nlo = pCurvature->p_base->nlo(pCurvature);
 *
 * // Demo: calculate the curvature V1 matrix.
 * losc_matrix *K = losc_matrix_create(nlo, nlo);
 * pCurvature->p_base->kappa(pCurvature->p_base, K);
 * @endcode
 * @see LoscCurvatureBase: see all the functions provided in the curvature
 * base structure.
 */

/**
 * @brief The C interface for LOSC curvature class of version 1.
 * @details Curvature version 1 is defined as \f$ \kappa \f$ in
 * Eq. 3 in the original LOSC paper (https://doi.org/10.1093/nsr/nwx11).
 * @see losc_curvature_v1_create(): create a LosCurvatureV1 object.
 * @see losc_curvature_v1_free(): free a LosCurvatureV1 object.
 */
typedef struct LoscCurvatureV1 {
    // points to a real `losc::CurvatureV1` object.
    _LoscCurvatureV1 *_p_v1;

    /**
     * @brief A pointer to a LoscCurvatureBase variable. Use this variable
     * to call functions that are provided in the LoscCurvatureBase.
     * @par Example
     * Demo of calling functions that provided in the LoscCurvatureBase.
     * Assume a variable `pCurvature` exists with type of `LoscCurvatureV1 *`.
     * @copydoc __code__demo_using_curvature_pbase
     */
    LoscCurvatureBase *p_base;

    // Note:
    // If `losc::CurvatureV1` has new public functions (new functions added in
    // `losc::CurvatureV1` not in its base class) in the future, add the
    // corresponding C interface with new function pointers below.

    /**
     * @brief Set parameter tau.
     * @param [in] tau parameter tau.
     * @par Example
     * @code
     * // Assume `pV1` is a existing variable with type `LoscCurvatureV1 *`.
     * // Demo: set parameter tau with 0.8.
     * pV1->set_tau(pV1, 0.8);
     * @endcode
     */
    void (*set_tau)(const LoscCurvatureV1 *self, double tau);
} LoscCurvatureV1;

/**
 * @class __param__LoscCurvatureBase
 * @param [in] dfa Type of DFA.
 * @param [in] df_pii Three-body integral \f$ \langle p |ii \rangle \f$
 * used in density fitting. Index `p` is for fitbasis and index `i`
 * is for LOs. The dimension of `df_pii` is `[nfitbasis, nlo]`.
 * @param [in] df_Vpq_inv Inverse of \f$ \langle p | 1/\mathbf{r}
 * | q \rangle \f$ matrix used in density fitting. Index `p` and `q`
 * are for LOs. The dimension of `df_Vpq_inv` is `[nfitbasis, nfitbasis]`.
 * @param [in] grid_lo LOs' value on grid points with dimension of
 * `[npts, nlo]`.
 * @param [in] grid_weight Coefficient vector for numerical integral on
 * grid with size of `npts`.
 */

/**
 * @class __note__LoscCurvatureBase
 * @note This function requires all the matrices are allocated in
 * advance. This could be memory consuming. For example, the `grid_lo`
 * matrix can be very large.
 */

/**
 * @brief Create a struct LoscCurvatureV1.
 * @copydoc __param__LoscCurvatureBase
 * @return Create a struct LoscCurvatureV1 and return the pointer.
 * @copydoc __note__LoscCurvatureBase
 */
LoscCurvatureV1 *losc_curvature_v1_create(const LoscDFAInfo *dfa,
                                          const losc_matrix *df_pii,
                                          const losc_matrix *df_Vpq_inv,
                                          const losc_matrix *grid_lo,
                                          const double *grid_weight);

/**
 * @brief Free a struct LoscCurvatureV1.
 * @param [in, out] ptr a pointer with type `LoscCurvatureV1 *`. At exit,
 * `ptr` is set to null.
 */
#define losc_curvature_v1_free(ptr) _losc_curvature_v1_free(&(ptr))
void *_losc_curvature_v1_free(LoscCurvatureV1 **pptr_self);

//**********************************************
// ==> Interface for `losc::CurvatureV1`
//**********************************************

// The incomplete struct in C-side that is `losc::CurvatureV2` in C++-side.
typedef struct _LoscCurvatureV2 _LoscCurvatureV2;

/**
 * @brief The C interface for LOSC curvature class of version 2.
 * @details The curvature matrix is defined as \f$\tilde{\kappa}\f$ in
 * Eq. 8 of the LOSC2 paper (J. Phys. Chem. Lett. 2020, 11, 4, 1528-1535).
 * @see losc_curvature_v2_create(): create a LosCurvatureV2 object.
 * @see losc_curvature_v2_free(): free a LosCurvatureV2 object.
 */
typedef struct LoscCurvatureV2 {
    // points to a real `losc::CurvatureV2` object.
    _LoscCurvatureV2 *_p_v2;

    /**
     * @brief A pointer to a LoscCurvatureBase variable. Use this variable
     * to call functions that are provided in the LoscCurvatureBase.
     * @par Example
     * Demo of calling functions that provided in the LoscCurvatureBase.
     * Assume a variable `pCurvature` exists with type of `LoscCurvatureV2 *`.
     * @copydoc __code__demo_using_curvature_pbase
     */
    LoscCurvatureBase *p_base;

    // Note
    // If `losc::CurvatureV2` has new public functions (new functions added in
    // `losc::CurvatureV2` not in its base class) in the future, add the
    // corresponding C interface with new function pointers below.

    /**
     * @brief Set parameter tau.
     * @par Example
     * @code
     * // Assume `pV2` is a existing variable with type `LoscCurvatureV2 *`.
     * // Demo: set parameter tau to be 0.8.
     * pV2->set_tau(pV2, 0.8);
     * @endcode
     */
    void (*set_tau)(const LoscCurvatureV2 *self, double tau);

    /**
     * @brief Set parameter zeta.
     * @par Example
     * @code
     * // Assume `pV2` is a existing variable with type `LoscCurvatureV2 *`.
     * // Demo: set parameter zeta to be 0.8.
     * pV2->set_zeta(pV2, 0.8);
     * @endcode
     */
    void (*set_zeta)(const LoscCurvatureV2 *self, double zeta);
} LoscCurvatureV2;

/**
 * @brief Constructor of LoscCurvatureV2.
 * @copydoc __param__LoscCurvatureBase
 * @return Create a struct LoscCurvatureV2 and return the pointer.
 * @copydoc __note__LoscCurvatureBase
 */
LoscCurvatureV2 *losc_curvature_v2_create(const LoscDFAInfo *dfa,
                                          const losc_matrix *df_pii,
                                          const losc_matrix *df_Vpq_inv,
                                          const losc_matrix *grid_lo,
                                          const double *grid_weight);

/**
 * @brief Free a struct LoscCurvatureV2.
 * @param [in, out] ptr a pointer with type `LoscCurvatureV1 *`. At exit,
 * `ptr` is set to null.
 */
#define losc_curvature_v2_free(ptr) _losc_curvature_v2_free(&(ptr))
void *_losc_curvature_v2_free(LoscCurvatureV2 **pptr_self);

#ifdef __cplusplus
}
#endif

#endif
