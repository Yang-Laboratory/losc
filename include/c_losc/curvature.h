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

// The incomplete struct in C-side that is `losc::DFAInfo` in C++-side.
typedef struct LoscDFAInfo LoscDFAInfo;

/**
 * Constructor of DFAInfo struct that stores the information
 * of associated DFA.
 *
 * @param gga_x: the total weight of GGA and LDA type exchange.
 * @param hf_x: the weight of HF exchange.
 * @param name: the description of the DFA.
 *
 * @return a pointer to DFAInfo struct.
 */
LoscDFAInfo *losc_dfa_info_create(double gga_x, double hf_x, const char *name);

/**
 * Destructor of DFAInfo struct.
 *
 * @param pptr_m: a pointer that points to a DFAInfo pointer. At exit,
 * `*pptr_m` is set to null.
 */
void _losc_dfa_info_free(LoscDFAInfo **pptr_m);
#define losc_dfa_info_free(ptr_m) _losc_DFAInfo_free(&(ptr_m))

//**********************************************
// ==> Interface for `losc::CurvatureBase`
//**********************************************

// The incomplete struct in C-side that is `losc::CurvatureBase` in C++-side.
typedef struct _LoscCurvatureBase _LoscCurvatureBase;

/**
 * @brief The C interface for LOSC curvature base class.
 * @details This struct exports the functions in `losc::CurvatureBase` in
 * the C++ code. This struct servers as helper class to avoid code duplication
 * for the implementation of C interface. There is no need for C users to
 * create a such struct in practice.
 */
typedef struct LoscCurvatureBase {
    // points to a real `losc::CurvatureBase` object.
    _LoscCurvatureBase *_p_base;

    // List interface of methods exposed to C users bellow.

    /**
     * Return the number of LOs.
     */
    size_t (*nlo)(const LoscCurvatureBase *self);

    /**
     * Return the number of nfitbasis.
     */
    size_t (*nfitbasis)(const LoscCurvatureBase *self);

    /**
     * Return the number of npts.
     */
    size_t (*npts)(const LoscCurvatureBase *self);

    /**
     * Compute the LOSC curvature matrix with an in-out way.
     *
     * @param K [in, out]: curvature matrix with dimension [nlo, nlo].
     * @note
     * The curvature matrix should be allocated in advance with an
     * expected dimension. Otherwise, it will throw an exception.
     */
    void (*kappa)(const LoscCurvatureBase *self, losc_matrix *K);
} LoscCurvatureBase;

//**********************************************
// ==> Interface for `losc::CurvatureV1`
//**********************************************

// The incomplete struct in C-side that is `losc::CurvatureV1` in C++-side.
typedef struct _LoscCurvatureV1 _LoscCurvatureV1;

/**
 * @brief The C interface for LOSC curvature class of version 1.
 * @details Curvature version 1 is defined as \f$ \kappa \f$ in
 * Eq. 3 in the original LOSC paper (https://doi.org/10.1093/nsr/nwx11).
 */
typedef struct LoscCurvatureV1 {
    // points to a real `losc::CurvatureV1` object.
    _LoscCurvatureV1 *_p_v1;

    /**
     * @brief This variable stores a bunch of function pointers that point to
     * `losc::CurvatureBase` public interface.
     *
     * @code
     * // First create pV1 that points to a LoscCurvatureV1 struct by using
     * // losc_curvature_v1_create function.
     * size_t nbasis = pV1->p_base->nbasis(pV1); // Get number of AOs.
     * size_t nlo = pV1->p_base->nlo(pV1); // Get number of LOs.
     * losc_matrix *K = losc_matrix_create(nlo, nlo); // allocate curvature
     * matrix. pV1->p_base(pV1, K); // calculate curvature matrix and store it
     * in K.
     * @endcode
     *
     * @see
     * LoscCurvatureBase: For how to all the supported functions in the base.
     *
     * @example
     * How to use base functions.
     */
    LoscCurvatureBase *p_base;

    /**
     * @note
     * If `losc::CurvatureV1` has new public functions (new functions added in
     * `losc::CurvatureV1` not in its base class) in the future, add the
     * corresponding C interface with new function pointers below.
     */

    /**
     * @brief Set parameter tau.
     */
    void (*set_tau)(const LoscCurvatureV1 *self, double tau);
} LoscCurvatureV1;

/**
 * Constructor of LoscCurvatureV1.
 *
 * @param dfa_info: The associated DFA information.
 * @param df_pii: The three-body integral \f$\langle p|ii \rangle\f$ matrix
 * used in density fitting. Dimension is [nfitbasis, nlo].
 * The integral is defined as the following
 * \f[
 * \langle \phi_p | \phi_i \phi_i \rangle
 * = \int \frac{\phi_p(\rm{r}) \phi_i(\rm{r'}) \phi_i(\rm{r'})}{|\rm{r}
 *   - \rm{r'}|} \mbox{d} \rm{r} \mbox{d} \rm{r'},
 * \f]
 * in which index \f$p\f$ is for fitbasis and index \f$i\f$ is for LOs.
 * @param df_Vpq_inv: Inverse of \f$V_{pq}\f$ matrix used in density fitting.
 * Dimension is [nfitbasis, nfitbasis].
 * @param grid_lo: Value of LOs on grid points. Dimension is
 * [npts, nlo]. You have to build the matrix by yourself. To construct
 * the matrix, you can do as the following,
 * @code
 * for (size_t ip = 0; ip < npts; ++ip) {
 *     for (size_t i = 0; i < nlo; ++i) {
 *         // set the value of i-th LO
 *         // basis on grid point ip as zero.
 *         (*grid_lo)(ip, i) = 0;
 *     }
 * }
 * @endcode
 * @param grid_weight: Weights of grid points. The size of the array should be
 * at least `npts`. You have to build the vector by yourself. To construct the
 * vector, you can do as the following,
 * @code
 * for (size_t ip = 0; ip < npts; ++ip) {
 *     (*grid_weight_)[ip] = 0.0; // set each grid coefficient as zero.
 * }
 * @endcode
 *
 * @return a pointer to LoscCurvatureV1 struct.
 */
LoscCurvatureV1 *losc_curvature_v1_create(const LoscDFAInfo *dfa_info,
                                          const losc_matrix *df_pii,
                                          const losc_matrix *df_Vpq_inv,
                                          const losc_matrix *grid_lo,
                                          const double *grid_weight);

/**
 * Destructor of LoscCurvatureV1.
 */
void *_losc_curvature_v1_free(LoscCurvatureV1 **pptr_self);
#define losc_curvature_v1_free(ptr_self) _losc_curvature_v1_free(&(ptr_self))

//**********************************************
// ==> Interface for `losc::CurvatureV1`
//**********************************************

// The incomplete struct in C-side that is `losc::CurvatureV2` in C++-side.
typedef struct _LoscCurvatureV2 _LoscCurvatureV2;

/**
 * @brief The C interface for LOSC curvature class of version 2.
 * @details The curvature matrix is defined as \f$\tilde{\kappa}\f$ in
 * Eq. 8 of the LOSC2 paper (J. Phys. Chem. Lett. 2020, 11, 4, 1528-1535).
 */
typedef struct LoscCurvatureV2 {
    // points to a real `losc::CurvatureV2` object.
    _LoscCurvatureV2 *_p_v2;

    /**
     * @brief This variable stores a bunch of function pointers that point to
     * `losc::CurvatureBase` public interface.
     *
     * @see
     * LoscCurvatureBase: For how to all the supported functions in the base.
     * LoscCurvatureV1: See similar examples of how to use base functions.
     */
    LoscCurvatureBase *p_base;

    /**
     * @note
     * If `losc::CurvatureV2` has new public functions (new functions added in
     * `losc::CurvatureV2` not in its base class) in the future, add the
     * corresponding C interface with new function pointers below.
     */

    /**
     * @brief Set parameter tau.
     */
    void (*set_tau)(const LoscCurvatureV2 *self, double tau);

    /**
     * @brief Set parameter zeta.
     */
    void (*set_zeta)(const LoscCurvatureV2 *self, double zeta);
} LoscCurvatureV2;

/**
 * Constructor of LoscCurvatureV2.
 *
 * @see
 * losc_curvature_v1_create: The interface of arguments is the same as
 * losc_curvature_v1_create function.
 */
LoscCurvatureV2 *losc_curvature_v2_create(const LoscDFAInfo *dfa_info,
                                          const losc_matrix *df_pii,
                                          const losc_matrix *df_Vpq_inv,
                                          const losc_matrix *grid_lo,
                                          const double *grid_weight);

/**
 * Destructor of LoscCurvatureV2.
 */
void *_losc_curvature_v2_free(LoscCurvatureV2 **pptr_self);
#define losc_curvature_v2_free(ptr_self) _losc_curvature_v2_free(&(ptr_self))

#ifdef __cplusplus
}
#endif

#endif
