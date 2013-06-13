#ifndef MOD_INIT_H
#define MOD_INIT_H

/* Main high level functions. */

void MOD_init (int nb_mod);
void MOD_clear (void);

/* Some global variables. */

extern int MOD_nb_mod_max, MOD_nb_modulis;
extern MOD_type MOD_current, MOD_current_inv;
extern MOD_type **MOD_mul_csts, **MOD_inverses, **MOD_single_inv;
extern MOD_type *MOD_primes, *MOD_primes_inv, *MOD_primes_inv_maj;
extern MOD_type *MOD_epsilon_lagrange, *MOD_epsilon_orig;
extern MOD_type *MOD_tree_sum, *MOD_tree_base;
extern MOD_type *MOD_residues, *MOD_newton_coef;

#endif /* MOD_INIT_H */
