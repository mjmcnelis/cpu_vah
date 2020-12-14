
#ifndef NEIGHBORCELLS_H_
#define NEIGHBORCELLS_H_

#include "Precision.h"
#include "DynamicalVariables.h"

// get neighbor cells of hydroydynamic variables along (x,y,n)
// to compute derivatives appearing in source and flux terms

// energy density (e)
void get_energy_density_neighbor_cells(const precision * const __restrict__ E, precision * const __restrict__ e1, int sim, int sip, int sjm, int sjp, int skm, int skp);


// fluid velocity (u)
void get_fluid_velocity_neighbor_cells(	fluid_velocity u_simm, fluid_velocity u_sim, fluid_velocity u_sip, fluid_velocity u_sipp,
										fluid_velocity u_sjmm, fluid_velocity u_sjm, fluid_velocity u_sjp, fluid_velocity u_sjpp,
										fluid_velocity u_skmm, fluid_velocity u_skm, fluid_velocity u_skp, fluid_velocity u_skpp,
										precision * const __restrict__ ui1, precision * const __restrict__ uj1, precision * const __restrict__ uk1,
										precision * const __restrict__ vxi, precision * const __restrict__ vyj, precision * const __restrict__ vnk, precision t2);


// hydro variables (q)
void get_hydro_variables_neighbor_cells(hydro_variables qm, hydro_variables qp, precision * const __restrict__ q);

#endif