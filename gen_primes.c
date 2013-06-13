
/* Output some primes to MOD_primes.h
 * It is used to reduce the computation time of MOD_init().
 *
 * The program takes an argument: the number of primes to generate.
 * They are generated down from the #define MOD_MAX_PRIME.
 *
 * Author: Sylvain.Pion@sophia.inria.fr, 1997.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define MOD_RECONSTRUCT_H

#include "MOD.h"
#include "MOD_consts.h"
#include "MOD_is_prime.h"

static MOD_type MOD_current, MOD_current_inv;

int main(int argc, char **argv)
{
    int i, p, howmany;

    if (argc < 2)
    {
        printf( "Usage: %s n\n"
                "n is the number of modulis you wanna generate.\n", argv[0]);
        exit(-1);
    };
    howmany = atoi(argv[1]);

    printf("#ifndef MOD_PRIMES_H\n"
           "#define MOD_PRIMES_H\n\n"
           "/* This file is generated automaticaly.\n"
           " * See gen_primes.c for comments.\n"
           " * %d modulis computed.\n"
	   " */\n\n"
           "#define HOWMANY_MODULIS %d\n\n"
           "static int MOD_pre_primes[] = {\n", howmany, howmany);

    /* Gives "howmany" prime numbers inferior to MOD_MAX_PRIME. */

    for(p=MOD_MAX_PRIME, i=howmany; (p>0) && (i>0); p--)
        if (MOD_is_prime(p) == 1)
        {
	    /* Little test to see if it can be used for "quick" MOD_reduce. */
	    MOD_type delta= ((MOD_type) p) * MOD_CST_2_26 + (p-1)/2;
	    MOD_current = (MOD_type) p;
	    MOD_current_inv = 1/MOD_current;
	    if ((delta - MOD_round(delta/MOD_current    )*MOD_current) !=
		(delta - MOD_round(delta*MOD_current_inv)*MOD_current))
		printf("%d is not good... rejected.\n", p);

            printf("%d", p);
            --i;
            printf( (i>0) ? ",\n" : "};\n#endif /* MOD_PRIMES_H */\n");
        };

    return 0;
}
