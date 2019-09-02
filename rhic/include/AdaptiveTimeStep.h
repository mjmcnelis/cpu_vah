
#ifndef ADAPTIVETIMESTEP_H_
#define ADAPTIVETIMESTEP_H_

#include "Precision.h"
#include "DynamicalVariables.h"
#include "Parameters.h"

precision compute_adaptive_time_step(precision t, precision dt_CFL, precision dt_source, precision dt_min);

precision compute_dt_CFL(precision t, lattice_parameters lattice, hydro_parameters hydro);

precision compute_dt_source(precision t, const hydro_variables * const __restrict__ q_prev, const hydro_variables * const __restrict__ q, const hydro_variables * const __restrict__ f, precision dt_prev, lattice_parameters lattice);

#endif