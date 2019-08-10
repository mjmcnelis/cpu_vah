#ifndef ADAPTIVETIMESTEP_H_
#define ADAPTIVETIMESTEP_H_

#include "Precision.h"

precision set_time_step(precision t_max, precision dt_min);

precision compute_max_time_step(precision t, const hydro_variables * const __restrict__ q, const precision * const __restrict__ e, const fluid_velocity * const __restrict__ u, const fluid_velocity * const __restrict__ up, int nx, int ny, int nz, precision dt, precision dt_prev, precision dx, precision dy, precision dn, precision etabar_const);

#endif