
#ifndef NEIGHBORCELLS_H_
#define NEIGHBORCELLS_H_

#include "Precision.h"

// get the neighbor cells of hydroydynamic variables along (x,y,n)
// used to compute the source and flux terms in the euler step


// primary variables of neighbor cells (e)
void get_primary_neighbor_cells(const precision * const __restrict__ E, precision * const __restrict__ e1, int sim, int sip, int sjm, int sjp, int skm, int skp);


// fluid velocity of neighbor cells (ut, ux, uy, un)
void get_u_neighbor_cells(const precision * const __restrict__ ut, const precision * const __restrict__ ux, const precision * const __restrict__ uy, const precision * const __restrict__ un, precision * const __restrict__ ui1, precision * const __restrict__ uj1, precision * const __restrict__ uk1, int sim, int sip, int sjm, int sjp, int skm, int skp);


// spatial velocity of neighbor cells (vx, vy, vn)
void get_v_neighbor_cells(const precision * const __restrict__ ut, const precision * const __restrict__ ux, const precision * const __restrict__ uy, const precision * const __restrict__ un, precision * const __restrict__ vxi, precision * const __restrict__ vyj, precision * const __restrict__ vnk, int simm, int sim, int sip, int sipp, int sjmm, int sjm, int sjp, int sjpp, int skmm, int skm, int skp, int skpp);


// conserved variables of neighbor cells (q)
void get_q_neighbor_cells(const precision * const __restrict__ q_s, precision * const __restrict__ qi1, precision * const __restrict__ qj1, precision * const __restrict__ qk1, precision * const __restrict__ qi2, precision * const __restrict__ qj2, precision * const __restrict__ qk2, int * r, int simm, int sim, int sip, int sipp, int sjmm, int sjm, int sjp, int sjpp, int skmm, int skm, int skp, int skpp);


#endif