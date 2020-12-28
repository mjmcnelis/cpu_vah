
#ifndef SOURCETERMS_H_
#define SOURCETERMS_H_

#include "Precision.h"
#include "AdaptiveTimeStep.h"
#include "Parameters.h"

void source_terms_aniso_hydro(precision * const __restrict__ S, const precision * const __restrict__ Q, precision e_s, precision lambda_s, precision aT_s, precision aL_s, precision t, const precision * const __restrict__ qi1, const precision * const __restrict__ qj1, const precision * const __restrict__ qk1, const precision * const __restrict__ e1, const precision * const __restrict__ ui1, const precision * const __restrict__ uj1, const precision * const __restrict__ uk1, precision ux, precision uy, precision un, precision ux_p, precision uy_p, precision un_p, precision dt_prev, precision dx, precision dy, precision dn, hydro_parameters hydro);

void source_terms_viscous_hydro(precision * const __restrict__ S, const precision * const __restrict__ Q, precision e_s, precision t, const precision * const __restrict__ qi1, const precision * const __restrict__ qj1, const precision * const __restrict__ qk1, const precision * const __restrict__ e1, const precision * const __restrict__ ui1, const precision * const __restrict__ uj1, const precision * const __restrict__ uk1, precision ux, precision uy, precision un, precision ux_p, precision uy_p, precision un_p, precision dt_prev, precision dx, precision dy, precision dn, hydro_parameters hydro);

#endif
