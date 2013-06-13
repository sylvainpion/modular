#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "_FPU.h"
#include "MOD.h"
#include "MOD_primes.h"
#include "MOD_is_prime.h"


#ifndef MOD_THREADS
MOD_type MOD_current, MOD_current_inv;
#endif

MOD_type **MOD_mul_csts, **MOD_inverses, **MOD_single_inv;
MOD_type *MOD_primes, *MOD_primes_inv, *MOD_primes_inv_maj;
MOD_type *MOD_epsilon_lagrange, *MOD_epsilon_orig;
MOD_type *MOD_newton_coef, *MOD_max_li, *MOD_max_co, *MOD_residues;
// MOD_type *MOD_tree_base, *MOD_tree_sum;
int MOD_nb_mod_max, MOD_nb_modulis;


/* This function must be called before any use of the modular package.
 * It does allocate memory, and a lot of precomputations.
 * You must tell the number of modulis, that you plan to use.
 */

void MOD_init (int nb_mod)
{
    int i, j, k, p;
    MOD_type tmp;
#ifdef MOD_THREADS
    MOD_type MOD_current, MOD_current_inv;
#endif

    /* Set the global variable. */

    MOD_nb_mod_max = nb_mod;

    /* Fill the arrays. */

    /* Computes "MOD_nb_mod_max" prime numbers inferior to MOD_MAX_PRIME.
     * And store them in MOD_primes[].
     */

    MOD_primes = (MOD_type *) malloc(sizeof(MOD_type) * MOD_nb_mod_max);
    for(i=0; i<MOD_MIN(HOWMANY_MODULIS,MOD_nb_mod_max); i++)
	MOD_primes[i] = MOD_pre_primes[i];

    p = (int) MOD_primes[i-1]-1;
    while (i<MOD_nb_mod_max)
    {
	if (MOD_is_prime(p) == 1)
	{
	    /* Quick test to see if it can be used for "quick" MOD_reduce.*/
	    MOD_type delta= ((MOD_type) p) * MOD_CST_2_26 + (p-1)/2;
	    MOD_current = (MOD_type) p;
	    MOD_current_inv = 1/MOD_current;
	    if ((delta - MOD_round(delta/MOD_current    )*MOD_current) !=
		(delta - MOD_round(delta*MOD_current_inv)*MOD_current))
		printf("%d is not good...\n", p);
	    MOD_primes[i] = (MOD_type) p;
	    ++i;
	};
	p--;
    };


    /* Compute the inverses of the m_i, rounding to nearest. */

    MOD_primes_inv = (MOD_type *) malloc(sizeof(MOD_type) * MOD_nb_mod_max);
    for(i=0; i<MOD_nb_mod_max; i++)
	MOD_primes_inv[i] = 1/MOD_primes[i];


    /* Compute the inverses of the m_i, rounding to infinity. */

    MOD_primes_inv_maj = (MOD_type *) malloc(sizeof(MOD_type) * MOD_nb_mod_max);
    _FPU_set_rounding_to_infinity();
    for(i=0; i<MOD_nb_mod_max; i++)
	MOD_primes_inv_maj[i] = 1/MOD_primes[i];
    _FPU_set_rounding_to_nearest();


    /* Compute the $w_i^{(j)} = ( \prod_{k!=j}^i m_k )^{-1} \mod m_j$. */

    MOD_inverses = (MOD_type **) malloc (sizeof(MOD_type *) * (1+MOD_nb_mod_max));
    MOD_inverses[0] = NULL;		/* It must never be accessed !!! */
    MOD_inverses[1] = (MOD_type *) malloc (sizeof(MOD_type) * (((MOD_nb_mod_max+2) * (MOD_nb_mod_max+1)) / 2));
    for (i=2; i<=MOD_nb_mod_max; i++)
	MOD_inverses[i] = MOD_inverses[i-1] + i-1 ;
    for (i=1; i<=MOD_nb_mod_max;i++)	/* When we use i modulis. */
	for (j=0; j<i; j++)	/* Compute w_i^{(j)}. */
	{
	    MOD_current = MOD_primes[j];
	    MOD_current_inv = MOD_primes_inv[j];
	    for(k=0,tmp=1; k<i; k++)
		tmp = (k==j)? tmp : MOD_mul(tmp, MOD_primes[k], MOD_current, MOD_current_inv);
	    MOD_inverses[i][j] = MOD_inv(tmp, MOD_current, MOD_current_inv);
	};
    

    /* Computes the m_j^{-1} [m_i]. (used by Newton, and MOD_deduction())
     * The simple inverses MOD_single_inv[i][j] = m_j^{-1} [m_i].
     * Only computed for j<i and 0<i<MOD_nb_mod_max.
     */

    MOD_single_inv = (MOD_type **) malloc (sizeof(MOD_type *) * MOD_nb_mod_max);
    MOD_single_inv[0] = NULL;	/* It must never be accessed !!! */
    MOD_single_inv[1] = (MOD_type *) malloc (sizeof(MOD_type) * (((MOD_nb_mod_max+1) * MOD_nb_mod_max)/2));
    for (i=2; i<MOD_nb_mod_max; i++)
	MOD_single_inv[i] = MOD_single_inv[i-1] + i-1;
    for (i=1; i<MOD_nb_mod_max;i++)
    {
	MOD_current = MOD_primes[i];
	MOD_current_inv = MOD_primes_inv[i];
	for (j=0; j<i; j++)
	    MOD_single_inv[i][j] = MOD_inv(MOD_primes[j], MOD_current, MOD_current_inv);
    };
    

    /* Computes the error bounds for Lagrange's method. */

    MOD_epsilon_lagrange = (MOD_type *) malloc (sizeof(MOD_type) * (MOD_nb_mod_max+1));
    for (i=0; i<=MOD_nb_mod_max; i++)
	MOD_epsilon_lagrange[i] = i*MOD_MANTISSA_ERROR;


    /* Computes the necessary bounds on det/m: MOD_epsilon_orig[]. */

    MOD_epsilon_orig = (MOD_type *) malloc (sizeof(MOD_type) * (MOD_nb_mod_max+1));
    for (i=0; i<=MOD_nb_mod_max; i++)
	MOD_epsilon_orig[i] = 1-i*MOD_MANTISSA_ERROR;


    /* Computes the constants by what we can multiply, 
     * in the MOD_relax_mul method: MOD_mul_csts[][].
     */
    
    MOD_mul_csts = (MOD_type **) malloc(sizeof(MOD_type *) * MOD_nb_mod_max);
    MOD_mul_csts[0] = (MOD_type *) malloc(sizeof(MOD_type) * MOD_nb_mod_max * (MOD_nb_mod_max+1));
    for (i=1; i<MOD_nb_mod_max; i++)
	MOD_mul_csts[i] = MOD_mul_csts[i-1] + (MOD_nb_mod_max+1);
    _FPU_set_rounding_to_zero();
    for (j=0; j<MOD_nb_mod_max; j++)
    {
	MOD_current = MOD_primes[j];
	MOD_current_inv = MOD_primes_inv[j];
	for (i=0; i<=MOD_nb_mod_max; i++)
	    MOD_mul_csts[j][i] = MOD_reduce(MOD_round(MOD_MAX_MUL/( (i==0)?1:i ) -2), MOD_current, MOD_current_inv);
    };
    _FPU_set_rounding_to_nearest();


    /* Now, we ask for memory space, without initilizing it. */
    
    /* To do the sums in a "tree-fashion". */
    /* MOD_tree_sum[] and MOD_tree_base[] */

    MOD_tree_sum  = (MOD_type *) malloc(sizeof(MOD_type) * MOD_nb_mod_max * 2);
    MOD_tree_base = (MOD_type *) malloc(sizeof(MOD_type) * MOD_nb_mod_max);
    

    /* To store the expr \mod m_i. MOD_residues[] */

    MOD_residues = (MOD_type *) malloc(sizeof(MOD_type) * MOD_nb_mod_max);


    /* To compute the coefficients of Newton/Knuth MOD_newton_coef[]. */

    MOD_newton_coef = (MOD_type *) malloc(sizeof(MOD_type) * MOD_nb_mod_max);


    if (MOD_newton_coef == NULL ||
	MOD_residues == NULL ||
	MOD_tree_sum == NULL ||
	MOD_tree_base == NULL ||
	MOD_mul_csts == NULL ||
	MOD_mul_csts[0] == NULL ||
	MOD_epsilon_orig == NULL ||
	MOD_epsilon_lagrange == NULL ||
	MOD_single_inv == NULL ||
	MOD_single_inv[1] == NULL ||
	MOD_inverses == NULL ||
	MOD_inverses[1] == NULL ||
	MOD_primes_inv_maj == NULL ||
	MOD_primes_inv == NULL ||
	MOD_primes == NULL)
    {
    	printf("Not enough memory available for modular precomputations.\n");
	exit(-1);
    };
}

/* This function must be called, if you want to free memory space used by
 * the modular package. You can then recall MOD_init() with different values.
 */

void MOD_clear (void)
{
    free(MOD_primes);
    free(MOD_primes_inv);
    free(MOD_primes_inv_maj);
    free(MOD_inverses[1]);
    free(MOD_inverses);	
    free(MOD_single_inv[1]);
    free(MOD_single_inv);
    free(MOD_epsilon_lagrange);
    free(MOD_epsilon_orig);
    free(MOD_mul_csts[0]);
    free(MOD_mul_csts);	
    free(MOD_tree_sum);	
    free(MOD_tree_base);
    free(MOD_residues);	
    free(MOD_newton_coef);	
}
