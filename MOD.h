#ifndef MOD_H
#define MOD_H

typedef double MOD_type;

#define MOD_CST_CUT		MOD_CST_3_51
#define MOD_MAX_PRIME		134217728             /* sqrt(2^54) */
#define MOD_MANTISSA_ERROR	(1/MOD_CST_2_53)
#define MOD_MAX_MUL		MOD_CST_2_52

#include "MOD_init.h"
#include "MOD_base_ops.h"
#include "MOD_reconstruct.h"

#endif /* MOD_H */
