/**
 * @file losc.h
 * @brief C interface for localized orbital scaling correction (LOSC) library.
 *
 */

#ifndef __LOSC_INTERFACE_C_LOSC_LOSC_H__
#define __LOSC_INTERFACE_C_LOSC_LOSC_H__

#ifdef __cplusplus
extern "C" {
#endif

#include "matrix.h"

typedef struct DFAInfo DFAInfo;
typedef struct CurvatureV1 CurvatureV1;
typedef struct CurvatureV2 CurvatureV2;

//***********************
// ==> Class DFAInfo <==
//***********************

/**
 * Constructor of DFAInfo struct that stores the information
 * of associated DFA.
 *
 * @param gga_x_wt: the weight of GGA type exchange.
 * @param hf_x_wt: the weight of HF exchange.
 * @param name: the description of the DFA.
 *
 * @return a pointer to DFAInfo struct.
 */
DFAInfo *losc_DFAInfo_create(double gga_x_wt, double hf_x_wt, const char *name);

/**
 * Destructor of DFAInfo struct.
 *
 * @param pptr_m: a pointer that points to a DFAInfo pointer. At exit,
 * `*pptr_m` is set to null.
 */
void _losc_DFAInfo_free(DFAInfo **pptr_m);
#define losc_DFAInfo_free(ptr_m) _losc_DFAInfo_free(&(ptr_m))

//***********************
// ==> Class CurvatureV1
//***********************

/**
 * Constructor of CurvatureV1 struct.
 *
 * @param dfa_info: The associated DFA information.
 * @param C_lo: LO coefficient matrix under AO. Dimension is [nbasis, nlo].
 * The relation between LOs and AOs is
 * \f[
 * \psi_p = C_{\mu p} \phi_{\mu}
 * \f]
 * in which $\psi_p$ is the \f$p\f$-th LO, \f$C_{\mu p}\f$ is the LO
 * coefficient matrix and \f$\phi_{\mu}\f$ is the \f$\mu\f$-th AO.
 * @param df_pii: The three-body integral \f$\langle p|ii \rangle\f$ matrix
 * used in density fitting. Dimension is [nfitbasis, nlo].
 * The integral is defined as the following
 * \f[
 * \langle \phi_p | \phi_i \phi_i \rangle
 * = \int \frac{\phi_p(\rm{r}) \phi_i(\rm{r'}) \phi_i(\rm{r'})}{|\rm{r}
 *   - \rm{r'}|} \mbox{d} \rm{r} \mbox{d} \rm{r'},
 * \f]
 * in which index \f$p\f$ is for fitbasis and index \f$i\f$ is for LOs.
 * @param df_Vpq_inverse: Inverse of \f$V_{pq}\f$ matrix used in density
 * fitting. Dimension is [nfitbasis, nfitbasis].
 * @param grid_basis_value: Value of AO basis on grid. Dimension is
 * [npts, nbasis]. You have to build the matrix by yourself. To construct
 * the matrix, you can do as the following,
 * @code
 * for (size_t ip = 0; ip < npts; ++ip) {
 *     for (size_t i = 0; i < nbasis; ++i) {
 *         // set the value of i-th AO
 *         // basis on grid point ip as zero.
 *         (*grid_basis_value_)(ip, i) = 0;
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
 * @return a pointer to CurvatureV1 struct.
 */
CurvatureV1 *losc_CurvatureV1_create(const DFAInfo *dfa_info,
                                     const losc_matrix *C_lo,
                                     const losc_matrix *df_pii,
                                     const losc_matrix *df_Vpq_inv,
                                     const losc_matrix *grid_basis_value,
                                     const double *grid_weight);

/**
 * Create a curvature matrix version 1.
 *
 * @param cur: the curvature version 1 struct object.
 * @return a pointer to curvature matrix version 1.
 */
losc_matrix *losc_CurvatureV1_kappa(const CurvatureV1 *cur);

//***********************
// ==> Class CurvatureV2
//***********************

/**
 * Constructor of CurvatureV2 struct.
 *
 * @note
 * The interface of arguments are the same as losc_CurvatureV1_create.
 */
CurvatureV2 *losc_CurvatureV2_create(const DFAInfo *dfa_info,
                                     const losc_matrix *C_lo,
                                     const losc_matrix *df_pii,
                                     const losc_matrix *df_Vpq_inv,
                                     const losc_matrix *grid_basis_value,
                                     const double *grid_weight);

/**
 * Create a curvature matrix version 2.
 *
 * @param cur: the curvature version 2 struct object.
 * @return a pointer to curvature matrix version 2.
 */
losc_matrix *losc_CurvatureV2_kappa(const CurvatureV2 *cur);

#ifdef __cplusplus
}
#endif

#endif
