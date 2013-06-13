#include <stdio.h>
#include <math.h>

#define MOD_RECONSTRUCT_H

#include "MOD.h"

/* Return 1 if p is prime, otherwise the minimal factor of p. */

int MOD_is_prime (const int p)
{
    int k, sq;

    if (p%2 == 0) return 2;		/* Test if p is even. */

    sq = 2 + (int) sqrt((double) p);	/* Compute sqrt(p) by excess. */

    for (k=3; k<=sq; k+=2)             	/* Divide by odd numbers. */
        if (p%k == 0) return k;

    return 1;
}
