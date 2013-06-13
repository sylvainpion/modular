#ifndef MOD_BASE_OPS_H
#define MOD_BASE_OPS_H

/* Basic operations concerning RNS. */

#include <math.h>	/* for fabs() */
#ifdef MOD_THREADS
#include <pthread.h>
#endif
#include "MOD_consts.h"	/* for MOD_CST_CUT */


/* Returns {-1,0,+1} depending on the sign of "a". */
#define MOD_SGN(a) ( ((a)>0) ? 1 : ( ((a)<0) ? -1 : 0) )

#define MOD_ABS(a) (fabs(a))
#define MOD_TRUE   (0==0)
#define MOD_FALSE  (0==1)

#define MOD_MIN(a,b) ((a)<(b) ? (a) : (b))
#define MOD_MAX(a,b) ((a)>(b) ? (a) : (b))

#ifndef MOD_THREADS
#define MOD_soft_reduce(a,b,c)		MOD_soft_reduce(a)
#define MOD_reduce(a,b,c)		MOD_reduce(a)
#define MOD_mul(a,b,c,d)		MOD_mul(a,b)
#define MOD_add(a,b,c,d)		MOD_add(a,b)
#define MOD_sub(a,b,c,d)		MOD_sub(a,b)
#define MOD_det2(a,b,c,d,e,f)		MOD_det2(a,b,c,d)
#define MOD_perm2(a,b,c,d,e,f)		MOD_perm2(a,b,c,d)
#define MOD_exp(a,b,c,d)		MOD_exp(a,b)
#define MOD_bezout(a,b,c)		MOD_bezout(a)
#define MOD_inv(a,b,c)			MOD_inv(a)
#define MOD_div(a,b,c,d)		MOD_div(a,b)
#endif

static inline MOD_type MOD_round (MOD_type a);
static inline MOD_type MOD_soft_reduce (MOD_type a, MOD_type, MOD_type);
static inline MOD_type MOD_reduce (MOD_type a, MOD_type, MOD_type);
static inline MOD_type MOD_mul (MOD_type a, MOD_type b, MOD_type, MOD_type);
static inline MOD_type MOD_add (MOD_type a, MOD_type b, MOD_type, MOD_type);
static inline MOD_type MOD_det2  (MOD_type a,MOD_type b,MOD_type c,MOD_type d, MOD_type, MOD_type);
static inline MOD_type MOD_perm2 (MOD_type a,MOD_type b,MOD_type c,MOD_type d, MOD_type, MOD_type);
static inline MOD_type MOD_exp (MOD_type a, int b, MOD_type, MOD_type);
static inline MOD_type MOD_bezout (MOD_type ri1, MOD_type, MOD_type);
static inline MOD_type MOD_inv (MOD_type a, MOD_type, MOD_type);
static inline MOD_type MOD_div (MOD_type a, MOD_type b, MOD_type, MOD_type);

/* Quick integer rounding, valid if a<2^51. */
static inline MOD_type MOD_round (MOD_type a)
{
	return (a + MOD_CST_CUT) - MOD_CST_CUT;	
}


/* Big modular reduction (e.g. after a det2). */
static inline MOD_type MOD_reduce (MOD_type a, MOD_type MOD_current, MOD_type MOD_current_inv)
{
/*	return a - MOD_current * MOD_round(a / MOD_current); */
	return a - MOD_current * MOD_round(a * MOD_current_inv);
}


/* Little modular reduction (e.g. after a simple addition). */
static inline MOD_type MOD_soft_reduce (MOD_type a, MOD_type MOD_current, MOD_type MOD_current_inv)
{
    MOD_type b = 2*a;
    return (b>MOD_current) ? a-MOD_current :
	  ((b<-MOD_current) ? a+MOD_current : a);
}


/* a*b */
static inline MOD_type MOD_mul (MOD_type a, MOD_type b, MOD_type MOD_current, MOD_type MOD_current_inv)
{
    MOD_type c = a*b;
    return MOD_reduce(c, MOD_current, MOD_current_inv);
}


/* a+b */
static inline MOD_type MOD_add (MOD_type a, MOD_type b, MOD_type MOD_current, MOD_type MOD_current_inv)
{
	MOD_type c = a+b;
	return MOD_soft_reduce(c, MOD_current, MOD_current_inv);
/*	return MOD_reduce(c, MOD_current, MOD_current_inv); */
}


/* a*d-b*c */
static inline MOD_type MOD_det2 (MOD_type a, MOD_type b, MOD_type c, MOD_type d, MOD_type MOD_current, MOD_type MOD_current_inv)
{
    MOD_type e = a*d-b*c;
    return MOD_reduce(e, MOD_current, MOD_current_inv);
}


/* a*d+b*c */
static inline MOD_type MOD_perm2 (MOD_type a, MOD_type b, MOD_type c, MOD_type d, MOD_type MOD_current, MOD_type MOD_current_inv)
{
    MOD_type e = a*d+b*c;
    return MOD_reduce(e, MOD_current, MOD_current_inv);
}


/* a^b	(0<=b<2^32) */
static inline MOD_type MOD_exp (MOD_type a, int b, MOD_type MOD_current, MOD_type MOD_current_inv)
{
    if (b == 0) { return 1.0; }
    else
    {
        MOD_type c = a;
        int i=27;
        int mask = 0x04000000;

        while ((b & mask) == 0) /* NB: not very good if b is small... */
        {
            i--;
            b <<= 1;
        };

    /* Taking c=a at the beginning vs c=1 => one multiplication less. */
        i--;
        b <<= 1;

        for (; i!=0; i--)
        {
            c = MOD_mul(c, c, MOD_current, MOD_current_inv);
            if ((b & mask) != 0)
                c = MOD_mul(c, a, MOD_current, MOD_current_inv);
            b <<= 1;
        };

        return c;
    };
}


/* a^-1, using Bezout (extended Euclidian algorithm). */
static inline MOD_type MOD_bezout (MOD_type ri1, MOD_type MOD_current, MOD_type MOD_current_inv)
{
    MOD_type bi = 0.0;
    MOD_type bi1 = 1.0;
    MOD_type ri = MOD_current;
    MOD_type p, tmp, tmp2;

    while (MOD_ABS(ri1) != 1.0)
    {
        p = MOD_round(ri/ri1);
        tmp = bi - p * bi1;
	tmp2 = ri - p * ri1;
        bi = bi1;
	ri = ri1;
	bi1 = tmp;
        ri1 = tmp2;
    };

    return ri1 * MOD_soft_reduce(bi1, MOD_current, MOD_current_inv);	/* Quicker !!!! */
/*    return (ri1>0) ? MOD_soft_reduce(bi1, MOD_current, MOD_current_inv) : -MOD_soft_reduce(bi1, MOD_current, MOD_current_inv); */
}


/* a^-1, computed with Bezout [or a^(p-2) (Fermat's little theorem).] */
static inline MOD_type MOD_inv (MOD_type a, MOD_type MOD_current, MOD_type MOD_current_inv)
{
/*  return MOD_exp(a, (int) MOD_current-2, MOD_current, MOD_current_inv); */
    return MOD_bezout(a, MOD_current, MOD_current_inv);
}


/* a/b */
static inline MOD_type MOD_div (MOD_type a, MOD_type b, MOD_type MOD_current, MOD_type MOD_current_inv)
{
    return MOD_mul(a, MOD_inv(b, MOD_current, MOD_current_inv), MOD_current, MOD_current_inv);
}

#endif /* MOD_BASE_OPS_H */
