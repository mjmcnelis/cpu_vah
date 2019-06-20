
#ifndef NEIGHBORCELLS_H_
#define NEIGHBORCELLS_H_

#include "DynamicalVariables.h"

// get the neighbor cells of hydroydynamic variables along (x,y,n)
// used to compute the source and flux terms in the euler step

// primary variables of neighbor cells (e, p)
void get_primary_neighbor_cells(const PRECISION * const __restrict__ E, const PRECISION * const __restrict__ P, PRECISION * const __restrict__ e1, PRECISION * const __restrict__ p1, int sim, int sip, int sjm, int sjp, int skm, int skp);


// fluid velocity of neighbor cells (ut, ux, uy, un)
void get_u_neighbor_cells(const PRECISION * const __restrict__ ut, const PRECISION * const __restrict__ ux, const PRECISION * const __restrict__ uy, const PRECISION * const __restrict__ un, PRECISION * const __restrict__ ui1, PRECISION * const __restrict__ uj1, PRECISION * const __restrict__ uk1, int sim, int sip, int sjm, int sjp, int skm, int skp);


// spatial velocity of neighbor cells (vx, vy, vn)
void get_v_neighbor_cells(const PRECISION * const __restrict__ ut, const PRECISION * const __restrict__ ux, const PRECISION * const __restrict__ uy, const PRECISION * const __restrict__ un, PRECISION * const __restrict__ vxi, PRECISION * const __restrict__ vyj, PRECISION * const __restrict__ vnk, int simm, int sim, int sip, int sipp, int sjmm, int sjm, int sjp, int sjpp, int skmm, int skm, int skp, int skpp);


// conserved variables of neighbor cells (q)
void get_q_neighbor_cells(const PRECISION * const __restrict__ q_s, PRECISION * const __restrict__ qi1, PRECISION * const __restrict__ qj1, PRECISION * const __restrict__ qk1, PRECISION * const __restrict__ qi2, PRECISION * const __restrict__ qj2, PRECISION * const __restrict__ qk2, int * r, int simm, int sim, int sip, int sipp, int sjmm, int sjm, int sjp, int sjpp, int skmm, int skm, int skp, int skpp);

#endif