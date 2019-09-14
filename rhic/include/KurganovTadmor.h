
#ifndef KURGANOVTADMOR_H_
#define KURGANOVTADMOR_H_

#include "Precision.h"
#include "Parameters.h"

void euler_step(precision t, const hydro_variables * const __restrict__ q_current, hydro_variables * const __restrict__ q_update,
const precision * const __restrict__ e_current, precision * const __restrict__ e_update, const fluid_velocity * const __restrict__ u_previous, const fluid_velocity * const __restrict__ u_current, fluid_velocity * const __restrict__ u_update, precision dt, precision dt_prev, lattice_parameters lattice, hydro_parameters hydro, int update, int RK2);

void evolve_hydro_one_time_step(int n, precision t, precision dt, precision dt_prev, lattice_parameters lattice, hydro_parameters hydro, int update);

#endif