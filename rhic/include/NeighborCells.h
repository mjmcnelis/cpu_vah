
#ifndef NEIGHBORCELLS_H_
#define NEIGHBORCELLS_H_

#include "Precision.h"

// get neighbor cells of hydroydynamic variables along (x,y,n)
// to compute derivatives appearing in source and flux terms 


// primary variables (e)
void get_primary_neighbor_cells(const precision * const __restrict__ E, precision * const __restrict__ e1, int sim, int sip, int sjm, int sjp, int skm, int skp);


void get_fluid_velocity_neighbor_cells(	fluid_velocity u_simm, fluid_velocity u_sim, fluid_velocity u_sip, fluid_velocity u_sipp, 
										fluid_velocity u_sjmm, fluid_velocity u_sjm, fluid_velocity u_sjp, fluid_velocity u_sjpp, 
										fluid_velocity u_skmm, fluid_velocity u_skm, fluid_velocity u_skp, fluid_velocity u_skpp,
										precision * const __restrict__ ui1, precision * const __restrict__ uj1, precision * const __restrict__ uk1,
										precision * const __restrict__ vxi, precision * const __restrict__ vyj, precision * const __restrict__ vnk, precision t2);

void get_conserved_neighbor_cells(	hydro_variables q_simm, hydro_variables q_sim, hydro_variables q_sip, hydro_variables q_sipp, 
									hydro_variables q_sjmm, hydro_variables q_sjm, hydro_variables q_sjp, hydro_variables q_sjpp, 
									hydro_variables q_skmm, hydro_variables q_skm, hydro_variables q_skp, hydro_variables q_skpp, 
									precision * const __restrict__ qi1, precision * const __restrict__ qj1, precision * const __restrict__ qk1, 
									precision * const __restrict__ qi2, precision * const __restrict__ qj2, precision * const __restrict__ qk2);

#endif