
/* Compute computation times using MOD-method.
 *
 * Computes a matrix of dimension "dim", with random coefficients on the
 * firsts "a" columns, the rest being the Identity plus randomly chosen
 * vectors from the "a" firsts columns: the determinant must be about 2^(a*b)
 * All coeffs are integers of b bits (b<54).
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MOD.h"
#include "MOD_det.h"

MOD_type *mat_space, **mat_ptr;

void usage(char **argv);
MOD_type **mat_gen(int dim, int a, int b);
MOD_type random_b (int b);
void print_maple(MOD_type **mat, int dim);


int main(int argc, char **argv)
{
	int dim, bits, iter, cols, test;
	MOD_type **mat_test;

	if (argc < 5) usage(argv);

	dim = atoi(argv[1]);
	cols = atoi(argv[2]);
	bits = atoi(argv[3]);
	iter = atoi(argv[4]);
	test = (argc == 6) ? atoi(argv[5]) : 1;

	mat_space = (MOD_type *) malloc (sizeof(MOD_type) * dim * dim);
	mat_ptr = (MOD_type **) malloc (sizeof(MOD_type*) * dim);

	if (mat_space == NULL || mat_ptr == NULL)
	{
	  printf("No memory to play with.\n");
	  exit(-1);
	};

	MOD_init(100);
	MOD_det_init(dim);

	mat_test = mat_gen(dim, cols, bits);
	// print_maple(mat_test, dim);

	if (test == 1)
	{
		printf("sign given by method MOD = %d\n",
			MOD_sign_det((MOD_type **) mat_test,dim));
	
		while (iter>1)
		{
			MOD_sign_det((MOD_type **) mat_test,dim);
			--iter;
		};
	}
	else
	{
                printf("nullity given by method MOD = %d\n",
                        MOD_det_is_null((MOD_type **) mat_test,dim));

		while (iter>1)
		{
                        MOD_det_is_null((MOD_type **) mat_test,dim);
			--iter;
		};
	};

	MOD_clear();
	MOD_det_clear();
	return 0;
}

void usage(char **argv)
{
	printf("Usage: %s dim cols bits iter [test]\n"
			"\n"
			"dim : dimension of the matrix\n"
			"cols: number of randomly chosen columns\n"
			"bits: number of bits of the entrees\n"
			"iter: number of iterations to have a usefull time\n"
			"test: 1 for MOD_sign_det() (default), and \n"
			"      0 for MOD_det_is_null()\n\n"
			"N.B.: The determinant is expected to have cols*bits bits.\n",
			argv[0]);
	exit(-1);
}

/* Print an output Maple-compatible, to compute the determinant. */

void print_maple(MOD_type **mat, int dim)
{
        int i,j;
        printf("det([");
        for (i=0; i<dim; i++)
        {
                printf("[");
                for (j=0; j<dim; j++)
                {
                        printf("%.0f", mat[i][j]);
                        printf( (j == (dim-1)) ? "]" : ", ");
                };
                printf( (i == (dim-1)) ? "]);\n" : ",\n");
        };
}

MOD_type **mat_gen(int dim, int a, int b)
{
        int i, j, k, diag;

        if (a==0) {a=dim-1; diag=0;}
        else diag=1;

        for (i=0; i<dim; i++)   mat_ptr[i] = mat_space + (i * dim);
        for (i=0; i<a; i++)
                for (j=0; j<dim; j++)
                        mat_ptr[i][j] = random_b (b);

        while (i<dim)
        {
                k = ((int) random_b (10)) % a;
                if (k<0) k += a;
                for (j=0; j<dim; j++)
            mat_ptr[i][j] = ((i==j) ? diag : 0) + mat_ptr[k][j];
                i++;
        };

        return mat_ptr;
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
