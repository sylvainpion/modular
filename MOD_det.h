#ifndef MOD_DET_H
#define MOD_DET_H

/* Select the sign reconstruction method you want to use for the determinant.
 * Lagrange's method is set by default.
 *
 * If you want to use Newton's [probabilistic or not] method.
 *  -1 -> deterministic
 *   i -> probability to fall ~= 2^(1-27i)
 */

#define MOD_LAGRANGE
/* #define MOD_NEWTON 2 */
/* #define MOD_RELAX_MUL */


/*************************************
 * End of user-configurable options. *
 *************************************/

int MOD_sign_det (double **mat, const int dim);
int MOD_det_is_null (double **mat, const int dim);
void MOD_det_init(int dim);
void MOD_det_clear();

extern int MOD_dim_max;

#endif /* MOD_DET_H */
