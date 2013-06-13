/*
 * Benchmarks (and compares) the different sign reconstruction methods.
 * And maybe approximate value reconstruction methods.
 *
 * Sylvain.Pion@sophia.inria.fr, '97, '98.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MOD.h"

// #define MOD_NEWTON -1	/* probability for the newton method */
#define MOD_NEWTON 2	/* probability for the newton method */

void residues_gen (void);
void usage (char **argv);
MOD_type random_b (int b);
void print_maple (void);

int nb_random;

int main(int argc, char **argv)
{
    char method;
    int newton_sig;
    unsigned int i;
    unsigned long long iter;

    if (argc < 5) usage(argv);

    MOD_nb_modulis = atoi(argv[1]);
    nb_random = atoi(argv[2]);
//    iter = 1 + 1000 * atoi(argv[3]);
    iter = 1 + atoi(argv[3]);
    method = argv[4][0];

    MOD_init (MOD_nb_modulis);

    residues_gen ();
    print_maple ();

    newton_sig = MOD_SGN(MOD_newton_coef[0]);
    for (i=1;i<MOD_nb_modulis;i++)
    	if (MOD_newton_coef[i] != 0)
	    newton_sig = MOD_SGN(MOD_newton_coef[i]);

    printf("sign given by method Newton = %d\n", newton_sig);
    printf("sign given by method lagrange = %d\n", MOD_lagrange());
    printf("sign given by method relax_mul= %d\n", MOD_relax_mul());

    printf("approximate value= %g\n", MOD_to_double());

    switch (method) {
    case 'l': while ((--iter)>1) MOD_lagrange(); break;
    case 'r': while ((--iter)>1) MOD_relax_mul(); break;
    case 'n': while ((--iter)>1)
    {
#if (MOD_NEWTON != -1)
	int newton_current_prob = MOD_NEWTON;
#endif
    	newton_sig=0;
	for (i=0; i<MOD_nb_modulis; i++)
	{
	    MOD_current = MOD_primes [i];
	    MOD_current_inv = MOD_primes_inv[i];
	    MOD_newton(i);
	    if (MOD_newton_coef[i] != 0)
	    {
#if (MOD_NEWTON != -1)
		newton_current_prob = MOD_NEWTON;
#endif
		newton_sig = (MOD_newton_coef[i] > 0) ? 1 : -1;
	    }
#if (MOD_NEWTON != -1)
            else --newton_current_prob;
	    if (newton_current_prob == 0)
		i = MOD_nb_modulis;
#endif
	};
    };
    break;
    default: printf("The method must be l[agrange], r[elax_mul] or n[ewton]\n");
    };

    MOD_clear();
    return 0;
}

void usage(char **argv)
{
	printf("Usage: %s nb_mods nb_randoms iter method\n\n"
		"nb_mods    : number of modulis\n"
		"nb_randoms : number of independant modulis\n"
		"iter       : number of iterations\n"
		"method     : l[agrange], r[elax_mul] or n[ewton]\n", argv[0]);
	exit(-1);
}

/* Print a Maple-compatible output, to compute the exact value. */

void print_maple()
{
    int i;
    // Le mieux est d'utiliser Newton_deduction() pour afficher...

    printf("evalf(");
    for (i=0; i<MOD_nb_modulis-1; i++)
	printf("%d + %d*(", (int) MOD_newton_coef[i], (int) MOD_primes[i]);
    printf("0");
    for (i=0; i<MOD_nb_modulis; i++) printf(")");
    printf(";\n");
}


void residues_gen ()
{
    unsigned int i;

/* The first nb_random ones are randomly chosen. */

    for (i=0; i<nb_random; i++)
    {
    	MOD_current = MOD_primes [i];
	MOD_current_inv = MOD_primes_inv[i];
    	MOD_residues[i] = MOD_reduce( random_b(53), MOD_current, MOD_current_inv );
	MOD_newton(i);
    };

/* The following ones are computed using the deduction method. */

    if (nb_random == 0)
    {
        MOD_residues[0] = MOD_newton_coef[0] = 0;
	++i;
    };

    while (i<MOD_nb_modulis)
    {
    	MOD_current = MOD_primes [i];
	MOD_current_inv = MOD_primes_inv[i];
	MOD_deduction(i);
	MOD_newton(i);
	i++;
    };
}

/* Returns a random integer on b bits, stored in a MOD_type.
 * NB: This is not optimized at all.
 */

MOD_type random_b (int b)
{
        MOD_type tmp=0.0;

        while (b>0)
        {
                tmp *= 2;
                tmp += (rand() / 3) % 2;
                --b;
        };

        tmp *= (((rand() / 3) % 2) == 0) ? 1 : -1;

        return tmp;
}
