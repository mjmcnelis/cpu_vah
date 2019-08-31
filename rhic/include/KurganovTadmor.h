
#ifndef KURGANOVTADMOR_H_
#define KURGANOVTADMOR_H_

#include "Precision.h"
#include "Parameters.h"

void euler_step(precision t, const hydro_variables * const __restrict__ Q_current, hydro_variables * const __restrict__ Q_update,
const precision * const __restrict__ e, const fluid_velocity * const __restrict__ u, const fluid_velocity * const __restrict__ up, precision dt, precision dt_prev, lattice_parameters lattice, hydro_parameters hydro, int update);

void evolve_hydro_one_time_step(int n, precision t, precision dt, precision dt_prev, lattice_parameters lattice, hydro_parameters hydro, int update);

#endif