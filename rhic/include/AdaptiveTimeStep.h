#ifndef ADAPTIVETIMESTEP_H_
#define ADAPTIVETIMESTEP_H_

#include "Precision.h"
#include "DynamicalVariables.h"

typedef struct
{
	precision dt_CFL;		// CFL bound
	precision dt_micro;		// relaxation time scale
	precision dt_rate;		// hydro variables' evolution rate

} hydro_time_scales;


precision set_time_step(hydro_time_scales dt_hydro, precision dt_min);


hydro_time_scales compute_hydro_time_scales(precision t, const hydro_variables * const __restrict__ q, const precision * const __restrict__ e, const fluid_velocity * const __restrict__ u, const fluid_velocity * const __restrict__ up, int nx, int ny, int nz, precision dt, precision dt_prev, precision dx, precision dy, precision dn, precision etabar_const, precision dt_min);

#endif