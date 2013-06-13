/* This file contains 2 main functions MOD_sign_det() and MOD_det_is_null().
 * They return the sign/nullity of the determinant of the matrix of doubles
 * (containing integer coefficients <2^54).
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef MOD_THREADS
#include <pthread.h>
pthread_t MOD_threads[MOD_THREADS];
#endif

#include "_FPU.h"
#include "MOD.h"
#include "MOD_det.h"

#ifndef MOD_THREADS
#define MOD_det2x2(a,b,c)	MOD_det2x2(a)
#define MOD_det3x3(a,b,c)	MOD_det3x3(a)
#define MOD_det4x4(a,b,c)	MOD_det4x4(a)
#define MOD_gauss(a,b,c,d)	MOD_gauss(a,b)
#endif

inline MOD_type MOD_det2x2(MOD_type **mat, MOD_type, MOD_type);
MOD_type MOD_det3x3(MOD_type **mat, MOD_type, MOD_type);
MOD_type MOD_det4x4(MOD_type **mat, MOD_type, MOD_type);
MOD_type MOD_gauss (MOD_type **mat, int dim, MOD_type, MOD_type);
int      MOD_det_needs_mods (double **mat, const int dim);

static MOD_type **MOD_mat_tmp, *MOD_mat_work_space;
static MOD_type *MOD_max_li, *MOD_max_co;
static double *MOD_det_bounds;
static int *MOD_nb_dim_mod, MOD_dim_max;
static int MOD_best_det_bounds[]={1,1,2,4,16,48,160,576,4096,14336,73728,327680,2985984,14929920,77635584};


/* One way to compute dynamicaly the number of modulis needed: Hadamard */
/* NB: problem: to compute the number of moduli needed, we need some
       precomputations to be done... so it's incremental (not optimal) */

int MOD_det_needs_mods (double **mat, const int dim)
{
    double tmp, co=MOD_det_bounds[dim], li=co;
    int i, j;
    int nb_col=MOD_nb_dim_mod[dim];
    int nb_lig=MOD_nb_dim_mod[dim];

    if (dim>MOD_dim_max)
    {
    	MOD_dim_max = dim;
    	MOD_det_clear();
	MOD_det_init(MOD_dim_max);
    };

    for (i=0; i<dim; i++)
        MOD_max_li[i] = MOD_max_co[i] = 0;

    for (i=0; i<dim; i++)
        for (j=0; j<dim; j++)
        {
            tmp = MOD_ABS(mat[i][j]);
            MOD_max_li[i] = MOD_MAX(MOD_max_li[i], tmp);
            MOD_max_co[j] = MOD_MAX(MOD_max_co[j], tmp);
        };

    _FPU_set_rounding_to_infinity();

    for (i=0; i<dim; i++)
    {
        li *= MOD_max_li[i];
	while (li*MOD_epsilon_orig[nb_lig] > .5)
        {
	    li *= MOD_primes_inv_maj[nb_lig++];
            if (nb_lig > MOD_nb_mod_max)
            {
              MOD_clear();
	      _FPU_set_rounding_to_nearest();
              MOD_init(MOD_nb_mod_max+10);
              return MOD_det_needs_mods (mat,dim);
            };
        };
        co *= MOD_max_co[i];
	while (co*MOD_epsilon_orig[nb_col] > .5)
        {
	    co *= MOD_primes_inv_maj[nb_col++];
            if (nb_col > MOD_nb_mod_max)
            {
              MOD_clear();
	      _FPU_set_rounding_to_nearest();
              MOD_init(MOD_nb_mod_max+10);
              return MOD_det_needs_mods (mat,dim);
            };
        };
    };

    _FPU_set_rounding_to_nearest();

    if ((li == 0.0) || (co == 0.0))
        return 0;

    return MOD_MIN(nb_col, nb_lig);
}

inline MOD_type MOD_det2x2(MOD_type **mat, MOD_type MOD_current, MOD_type MOD_current_inv)
{
    return MOD_det2(    mat[0][0], mat[0][1],
			mat[1][0], mat[1][1], MOD_current, MOD_current_inv);
}

MOD_type MOD_det3x3(MOD_type **mat, MOD_type MOD_current, MOD_type MOD_current_inv)
{
    MOD_type a = MOD_det2(  mat[1][1], mat[1][2],
			    mat[2][1], mat[2][2], MOD_current, MOD_current_inv);
    MOD_type b = MOD_det2(  mat[0][1], mat[0][2],
			    mat[2][1], mat[2][2], MOD_current, MOD_current_inv);
    MOD_type c = MOD_det2(  mat[0][1], mat[0][2],
			    mat[1][1], mat[1][2], MOD_current, MOD_current_inv);
    MOD_type d = MOD_det2(  mat[0][0], mat[1][0],
			    b, a, MOD_current, MOD_current_inv);
    a = MOD_mul(mat[2][0], c, MOD_current, MOD_current_inv);
    return MOD_add(a, d, MOD_current, MOD_current_inv);
}

MOD_type MOD_det4x4(MOD_type **mat, MOD_type MOD_current, MOD_type MOD_current_inv)
{
    MOD_type p0 = MOD_det2( mat[0][0], mat[1][0],
			    mat[0][1], mat[1][1], MOD_current, MOD_current_inv);
    MOD_type p1 = MOD_det2( mat[0][0], mat[2][0],
			    mat[0][1], mat[2][1], MOD_current, MOD_current_inv);
    MOD_type p2 = MOD_det2( mat[0][0], mat[3][0],
			    mat[0][1], mat[3][1], MOD_current, MOD_current_inv);
    MOD_type p3 = MOD_det2( mat[1][0], mat[2][0],
			    mat[1][1], mat[2][1], MOD_current, MOD_current_inv);
    MOD_type p4 = MOD_det2( mat[1][0], mat[3][0],
			    mat[1][1], mat[3][1], MOD_current, MOD_current_inv);
    MOD_type p5 = MOD_det2( mat[2][0], mat[3][0],
			    mat[2][1], mat[3][1], MOD_current, MOD_current_inv);
    MOD_type q0 = MOD_det2( mat[0][2], mat[1][2],
			    mat[0][3], mat[1][3], MOD_current, MOD_current_inv);
    MOD_type q1 = MOD_det2( mat[0][2], mat[2][2],
			    mat[0][3], mat[2][3], MOD_current, MOD_current_inv);
    MOD_type q2 = MOD_det2( mat[0][2], mat[3][2],
			    mat[0][3], mat[3][3], MOD_current, MOD_current_inv);
    MOD_type q3 = MOD_det2( mat[1][2], mat[2][2],
			    mat[1][3], mat[2][3], MOD_current, MOD_current_inv);
    MOD_type q4 = MOD_det2( mat[1][2], mat[3][2],
			    mat[1][3], mat[3][3], MOD_current, MOD_current_inv);
    MOD_type q5 = - MOD_det2( mat[2][2], mat[3][2], /* This one is opposed. */
			    mat[2][3], mat[3][3], MOD_current, MOD_current_inv);

    /* Then do: p0.(-q5) - p1.q4 + p2.q3 + p3.q2 - p4.q1 + p5.q0. */

    MOD_type a = MOD_det2(p2,p1,q4,q3, MOD_current, MOD_current_inv);
    MOD_type b = MOD_det2(p3,p0,q5,q2, MOD_current, MOD_current_inv);
    MOD_type c = MOD_det2(p5,p4,q1,q0, MOD_current, MOD_current_inv);

    return MOD_add(MOD_add(a,b, MOD_current, MOD_current_inv),c, MOD_current, MOD_current_inv);
}


/* Implementation of the Gauss algorithm, using RNS.
 *
 * NB1: Does not check if dimension == 1 or 0 (maybe crash :-).
 * NB2: Works only if MOD_current is prime.
 */

MOD_type MOD_gauss (MOD_type **mat, int dim, MOD_type MOD_current, MOD_type MOD_current_inv)
{
    switch (dim)
    {
    	case 2:	return MOD_det2x2(mat, MOD_current, MOD_current_inv);
    	case 3:	return MOD_det3x3(mat, MOD_current, MOD_current_inv);
    	case 4:	return MOD_det4x4(mat, MOD_current, MOD_current_inv);
    	default:
	{
	int i, j, pivot;
	MOD_type numer = 1.0;
	MOD_type denom = 1.0;
	MOD_type denom_partial = 1.0;	/* partial product of the pivots. */
	MOD_type pivot_value;		/* caches the value of the pivot. */

	/* Gauss is "recursif terminal", so I prefer a loop. */

	while (dim != 4)
	{
	    /* looking for the pivot. */

	    pivot = 0;
	    while (( pivot_value = mat[pivot][0] ) == 0)
	    {
		if ((++pivot) == dim) return 0;
		denom = -denom;
	    };

	    denom_partial = MOD_mul(denom_partial, pivot_value, MOD_current, MOD_current_inv);
	    denom = MOD_mul(denom, denom_partial, MOD_current, MOD_current_inv);

	    for (i=0; i<dim; i++)
		if (i != pivot)
		    for (j=1; j<dim; j++)
			mat[i][j] = MOD_det2(   pivot_value, mat[pivot][j],
						mat[i][0]  , mat[i][j], MOD_current, MOD_current_inv);
	    --dim;
	    for (i=0; i<dim; i++)
		mat[i] = 1 + mat[i + (i>=pivot)];
	};

	numer = MOD_det4x4(mat, MOD_current, MOD_current_inv);
	denom_partial = MOD_mul(denom_partial,denom_partial, MOD_current, MOD_current_inv);
	denom = MOD_mul (denom, denom_partial, MOD_current, MOD_current_inv);

	/* Returns the result of the final division. */
	return MOD_div (numer, denom, MOD_current, MOD_current_inv);
	};
    };
}

#ifdef MOD_THREADS
typedef struct {
  int thr, dim;
  MOD_type **mat;
} pass_thr;

void *MOD_gauss_thr (void *pass_arg)
{
  int i,j,k;
  int thr = ((pass_thr *) pass_arg)->thr;
  int dim = ((pass_thr *) pass_arg)->dim;
  MOD_type **mat = ((pass_thr *) pass_arg)->mat;
  MOD_type MOD_current, MOD_current_inv;

  for (i=thr; i<MOD_nb_modulis; i+=MOD_THREADS)
  {
    MOD_mat_tmp[thr*MOD_THREADS] = MOD_mat_work_space + dim*dim*thr;
    for (j=1+thr*MOD_THREADS; j<dim+thr*MOD_THREADS; j++)
      MOD_mat_tmp[j] = MOD_mat_tmp[j-1] + dim;
    MOD_current = MOD_primes[i];
    MOD_current_inv = MOD_primes_inv[i];

    for (j=0; j<dim; j++)
      for (k=0; k<dim; k++)
        MOD_mat_tmp[j+thr*MOD_THREADS][k] = MOD_reduce(mat[j][k], MOD_current, MOD_current_inv);

    MOD_residues[i] = MOD_gauss (MOD_mat_tmp + thr*MOD_THREADS, dim, MOD_current, MOD_current_inv);

// Sometimes it works, sometimes it doesn't...  I disbeleive it !!!

  };
  return NULL;
}
#endif

int MOD_sign_det (double **mat, const int dim)
{
#ifndef MOD_THREADS
  int i,j,k;
#else
  int thr;
  pass_thr tab_thr[MOD_THREADS];
#endif
#ifdef MOD_NEWTON
  int newton_sig=0;
  int newton_current_prob = MOD_NEWTON;
#endif

/* Determine a bound on the deterninant => number of modulis needed. */

  MOD_nb_modulis = MOD_det_needs_mods(mat,dim);
  if (MOD_nb_modulis == 0) return 0;

/* Computes all (det % m_i). */

#ifndef MOD_THREADS
  for (i=0; i<MOD_nb_modulis; i++)
  {
    MOD_mat_tmp[0] = MOD_mat_work_space;
    for (j=1; j<dim; j++)
      MOD_mat_tmp[j] = MOD_mat_tmp[j-1] + dim;
    MOD_current = MOD_primes[i];
    MOD_current_inv = MOD_primes_inv[i];

    /* Compute the matrix of coefs % MOD_current */
    for (j=0; j<dim; j++)
      for (k=0; k<dim; k++)
	MOD_mat_tmp[j][k] = MOD_reduce(mat[j][k], MOD_current, MOD_current_inv); 
    /* Use Gauss to compute the determinant modulo MOD_current. */
    MOD_residues[i] = MOD_gauss (MOD_mat_tmp,dim, MOD_current, MOD_current_inv);

#ifdef MOD_NEWTON
    /* Compute Newton's coefficient (only useful for Newton's method). */
    MOD_newton(i);
    if (MOD_newton_coef[i] != 0)
    {
      newton_current_prob = MOD_NEWTON;
      newton_sig = (MOD_newton_coef[i] > 0) ? 1 : -1;
    }
    else	/* Probabilistic variant. */
    if ((--newton_current_prob) == 0)
      return newton_sig;
#endif /* MOD_NEWTON */
  };
#else /* MOD_THREADS */
#ifdef MOD_NEWTON
#  error "Newton is not compatible with Threads. Sorry !"
#endif

	/* Let the threads go ! */
  for (thr=0; thr<MOD_THREADS; thr++)
  {
    tab_thr[thr].thr = thr;
    tab_thr[thr].dim = dim;
    tab_thr[thr].mat = mat;
    pthread_create(&MOD_threads[thr], NULL, *MOD_gauss_thr, (void *) (tab_thr + thr));
  };

        /* Wait for the threads to finish. */
  for (thr=0; thr<MOD_THREADS; thr++)
    pthread_join(MOD_threads[thr],NULL);

#endif /* MOD_THREADS */

    /* Computes the sign. */
#ifdef MOD_NEWTON
  return newton_sig;
#endif
#ifdef MOD_LAGRANGE
  return MOD_lagrange();
#endif
#ifdef MOD_RELAX_MUL
  return MOD_relax_mul();
#endif
}


/* Same, but check only if ==0 or not.
 * It's far quicker in average, if det !=0.
 * Be carefull, it returns MOD_TRUE or MOD_FALSE, not 0 or 1,
 * whether det==0 or not.
 */

int MOD_det_is_null (double **mat, const int dim)
{
    int i,j,k;
#ifdef MOD_THREADS
    MOD_type MOD_current, MOD_current_inv;
#endif

/* Determine a bound on the deterninant => number of modulis needed. */

    MOD_nb_modulis = MOD_det_needs_mods(mat,dim);
    if (MOD_nb_modulis == 0) return 0;

/* Computes all (det % m_i). */

    for (i=0; i<MOD_nb_modulis; i++)
    {
        MOD_mat_tmp[0] = MOD_mat_work_space;
	for (j=1; j<dim; j++)
	    MOD_mat_tmp[j] = MOD_mat_tmp[j-1] + dim;
	MOD_current = MOD_primes[i];
	MOD_current_inv = MOD_primes_inv[i];
    /* Compute the matrix of coefs % MOD_current */
	for (j=0; j<dim; j++)
	    for (k=0; k<dim; k++)
		MOD_mat_tmp[j][k] = MOD_reduce(mat[j][k], MOD_current, MOD_current_inv);
    /* Use Gauss to compute the determinant modulo MOD_current.
     * And end if possible.
     */
	if (MOD_gauss (MOD_mat_tmp,dim, MOD_current, MOD_current_inv) != 0) return MOD_FALSE;
    };
    return MOD_TRUE;	/* All det \mod m_i are null. */
}


void MOD_det_init(int dim)
{
    int i, j;

    MOD_dim_max = dim;

    /* Reserve some memory to compute the determinants modulo m_i.
     * MOD_mat_work_space[] and MOD_mat_tmp[].
     */

#ifdef MOD_THREADS
    MOD_mat_work_space = (MOD_type *) malloc(sizeof(MOD_type) *MOD_dim_max*MOD_dim_max*MOD_THREADS);
    MOD_mat_tmp = (MOD_type **) malloc(sizeof(MOD_type *) * MOD_dim_max*MOD_THREADS);
#else
    MOD_mat_work_space = (MOD_type *) malloc(sizeof(MOD_type) *MOD_dim_max*MOD_dim_max);
    MOD_mat_tmp = (MOD_type **) malloc(sizeof(MOD_type *) * MOD_dim_max);
#endif
    if ( (MOD_mat_work_space == NULL) || (MOD_mat_tmp == NULL) ) exit(-1);

    /* Needed by function MOD_det_needs_mods(). */

    MOD_max_li = (MOD_type *) malloc (sizeof(MOD_type) * MOD_dim_max);
    MOD_max_co = (MOD_type *) malloc (sizeof(MOD_type) * MOD_dim_max);
    if ( (MOD_max_li == NULL) || (MOD_max_co == NULL) ) exit(-1);

    _FPU_set_rounding_to_nearest();

    /* Computes the best bounds known for determinants.
     * It's predefined til dimension 14, then Hadamard's bound is used.
     * NB: MOD_nb_dim_mod contains the number of moduli needed,
     * and MOD_det_bounds[] contains the rest quantity.
     */

    MOD_det_bounds = (double *) malloc(sizeof(double) * (MOD_dim_max+1));
    MOD_nb_dim_mod = (int *) malloc(sizeof(int) * (MOD_dim_max+1));
    if ( (MOD_det_bounds == NULL) || (MOD_nb_dim_mod == NULL) ) exit(-1);
    _FPU_set_rounding_to_infinity();

    for (i=0; i<MOD_MIN(sizeof(MOD_best_det_bounds)/sizeof(int),MOD_dim_max+1); i++)
    {
        MOD_det_bounds[i] = MOD_best_det_bounds[i];
        MOD_nb_dim_mod[i] = 0;
    };
    while (i<=MOD_dim_max)
    {
    /* Compute Hadamard's bound, ie d^(d/2). */

        MOD_nb_dim_mod[i] = 0;
        MOD_det_bounds[i] = 1.0;
        for(j=1; j<=(i/2);j++)
        {
            MOD_det_bounds[i] *= i;
            if (MOD_det_bounds[i]/MOD_primes[MOD_nb_dim_mod[i]] > 1.0)
            {
                MOD_det_bounds[i] /= MOD_primes[MOD_nb_dim_mod[i]];
                ++MOD_nb_dim_mod[i];
		if (MOD_nb_dim_mod[i] >= MOD_nb_mod_max)
		{
		    MOD_clear();
		    _FPU_set_rounding_to_nearest();
		    MOD_init(MOD_nb_mod_max+10);
		    _FPU_set_rounding_to_infinity();
		};
            };
        };
        if (i%2 != 0)
        {
            MOD_det_bounds[i] *= sqrt(i);
            if (MOD_det_bounds[i]/MOD_primes[MOD_nb_dim_mod[i]] > 1.0)
            {
                MOD_det_bounds[i] /= MOD_primes[MOD_nb_dim_mod[i]];
                ++MOD_nb_dim_mod[i];
                if (MOD_nb_dim_mod[i] >= MOD_nb_mod_max)
                {
                    MOD_clear();
                    _FPU_set_rounding_to_nearest();
                    MOD_init(MOD_nb_mod_max+10);
                    _FPU_set_rounding_to_infinity();
                };
            };
        };
        ++i;
    };

    _FPU_set_rounding_to_nearest();
}

void MOD_det_clear()
{
    free(MOD_mat_work_space);
    free(MOD_mat_tmp);
    free(MOD_max_li);
    free(MOD_max_co);
    free(MOD_det_bounds);
    free(MOD_nb_dim_mod);
}
