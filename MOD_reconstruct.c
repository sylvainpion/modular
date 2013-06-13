/* Operations on several residues (sign reconstruction, etc...) */

#include "MOD.h"
#include "MOD_reconstruct.h"

MOD_type *MOD_tree_base, *MOD_tree_sum;

int MOD_lagrange(void)
{
#ifdef MOD_THREADS
    MOD_type MOD_current, MOD_current_inv;
#endif
    MOD_type sum=0.0;
    int j, i=MOD_nb_modulis;

    while ( (i!=0) && (MOD_ABS(sum) <= MOD_epsilon_lagrange[i]) )
    {
        for (j=0; j<i; j++)
        {
            MOD_current = MOD_primes[j];
            MOD_current_inv = MOD_primes_inv[j];
            MOD_tree_sum[j] = MOD_mul(MOD_residues[j], MOD_inverses[i][j], MOD_current, MOD_current_inv) * MOD_current_inv;
        };
        j=0; --i;
        while (j<i)
        {
            MOD_tree_sum[1+i+j] = MOD_tree_sum[2*j] + MOD_tree_sum[2*j+1];
            j++;
        };
        sum = MOD_tree_sum[2*i] - MOD_round(MOD_tree_sum[2*i]);
    };
    return MOD_SGN(sum);
}

int MOD_relax_mul(void)
{
#ifdef MOD_THREADS
    MOD_type MOD_current, MOD_current_inv;
#endif
    MOD_type sum = 0.0;
    int j, flag_is_zero = MOD_TRUE;

    for (j=0; j<MOD_nb_modulis; j++)
    {
	MOD_current = MOD_primes[j];
	MOD_current_inv = MOD_primes_inv[j];
	MOD_tree_base[j] = MOD_mul(MOD_residues[j],
				MOD_inverses[MOD_nb_modulis][j],
				MOD_current, MOD_current_inv);
	MOD_tree_sum[j] = MOD_tree_base[j] * MOD_current_inv;
	flag_is_zero = (MOD_residues[j] == 0) && flag_is_zero;
    };

    if (flag_is_zero) return 0;

    while (MOD_TRUE)
    {
        j=0;
        while (j<(MOD_nb_modulis-1))
        {
            MOD_tree_sum[MOD_nb_modulis+j] = MOD_tree_sum[2*j] + MOD_tree_sum[2*j+1];
            j++;
        };
        sum = MOD_tree_sum[2*MOD_nb_modulis-2] -
		MOD_round(MOD_tree_sum[2*MOD_nb_modulis-2]);

        if (MOD_ABS(sum) > MOD_epsilon_lagrange[MOD_nb_modulis])
            return MOD_SGN(sum);

        for (j=0; j<MOD_nb_modulis; j++)
        {
            MOD_current = MOD_primes[j];
            MOD_current_inv = MOD_primes_inv[j];
            MOD_tree_base[j] = MOD_mul(MOD_tree_base[j],
				MOD_mul_csts[j][MOD_nb_modulis],
				MOD_current, MOD_current_inv);
            MOD_tree_sum[j] = MOD_tree_base[j] * MOD_current_inv;
        };
    };
}


/* Computes (and stores), the i'th coeff of Newton/Knuth, the one before
 * being already computed (and the residues stored in MOD_residues[]).
 * Computes the i'th Newton v_i.
 */
void MOD_newton (const int i)
{
    int j=0;
    MOD_newton_coef[i]=MOD_residues[i];

    while (j<i)
    {
        MOD_newton_coef[i] = MOD_reduce((MOD_newton_coef[i] - MOD_newton_coef[j]) * MOD_single_inv[i][j], MOD_current, MOD_current_inv);
        j++;
    };
}


/* Returns a (good, but maybe not the best) double approximation of the value,
 * represented by its residues.  Simple Horner scheme.
 * NB: The Newton's coefficients have to be computed before.
 * NB2: Overflow is not managed.
 */
double MOD_to_double ()
{
  double v = 0.0;
  int i;

  for (i=MOD_nb_modulis-2; i>=0; i--)
  {
    printf("i= %d\n", i);
    v = v*MOD_primes[i] + MOD_newton_coef[i];
  };
  return v;
}


/* Computes MOD_residues[i] from MOD_newton_coef[0->i-1].
 * (assuming det is < prod(m0...m_i-1))
 * Deduce "expr mod m_i" (MOD_residues[i]), from MOD_newton_coef[0 -> i-1].
 */
void MOD_deduction (const int i)
{
    int j=i-1;

    MOD_residues[i] = MOD_newton_coef[j];

    while ((j--)>0)
        MOD_residues[i] = MOD_reduce(MOD_residues[i] * MOD_primes[j] + MOD_newton_coef[j], MOD_current, MOD_current_inv);
}


/* Given a positive bound, compute the necessary number of moduli. */
int MOD_fp_needs_mods (double bound)
{
        int nb=0;
        double tmp = bound;

        while ( tmp * MOD_epsilon_orig[nb] > .5 )
        {
                tmp *= MOD_primes_inv_maj[nb++];
                if (nb > MOD_nb_mod_max)
                {
                   MOD_clear();
                   MOD_init(MOD_nb_mod_max+10);
                   return MOD_fp_needs_mods (bound);
                };
        };

        return nb;
}
