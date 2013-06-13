#ifndef MOD_RECONSTRUCT_H
#define MOD_RECONSTRUCT_H

/* Operations on several residues (sign reconstruction, etc...) */

extern MOD_type *MOD_tree_base, *MOD_tree_sum;

/* The different methods used to compute the sign. */

int MOD_lagrange (void);
int MOD_relax_mul (void);
void MOD_newton (const int i);
void MOD_deduction (const int i);
int MOD_fp_needs_mods (double bound);

/* Simple deduction method */
double MOD_to_double ();

#endif /* MOD_RECONSTRUCT_H */
