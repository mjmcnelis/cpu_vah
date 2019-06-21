
#ifndef SOURCETERMS_H_
#define SOURCETERMS_H_

#include "Precision.h"


void source_terms(precision * const __restrict__ S, const precision * const __restrict__ Q, precision e_s, precision t, const precision * const __restrict__ qi1, const precision * const __restrict__ qj1, const precision * const __restrict__ qk1, const precision * const __restrict__ e1, const precision * const __restrict__ ui1, const precision * const __restrict__ uj1, const precision * const __restrict__ uk1, precision ut, precision ux, precision uy, precision un, precision ut_p, precision ux_p, precision uy_p, precision un_p, precision dt, precision dx, precision dy, precision dn, precision etabar);


#endif
