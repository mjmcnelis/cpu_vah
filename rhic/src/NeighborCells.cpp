
#include "../include/NeighborCells.h"
#include "../include/Precision.h"


void get_primary_neighbor_cells(const precision * const __restrict__ E, precision * const __restrict__ e1, int sim, int sip, int sjm, int sjp, int skm, int skp)
{
	// energy density of neighbor cells stored in e1 = (e_sim, e_sip, e_sjm, e_sjp, e_skm, e_skp)
	e1[0] = E[sim];		// i-1
	e1[1] = E[sip];		// i+1
	e1[2] = E[sjm];		// j-1
	e1[3] = E[sjp];		// j+1
	e1[4] = E[skm];		// k-1
	e1[5] = E[skp];		// k+1
}


void get_u_neighbor_cells(const precision * const __restrict__ ut, const precision * const __restrict__ ux, const precision * const __restrict__ uy, const precision * const __restrict__ un, precision * const __restrict__ ui1, precision * const __restrict__ uj1, precision * const __restrict__ uk1, int sim, int sip, int sjm, int sjp, int skm, int skp)
{
	ui1[0] = ut[sim];		// fluid velocity of neighbor cells
	ui1[1] = ut[sip];		// along x [i-1, i+1] stored in ui1
							// (ut, ux, uy, un)
	ui1[2] = ux[sim];
	ui1[3] = ux[sip];

	ui1[4] = uy[sim];
	ui1[5] = uy[sip];

	ui1[6] = un[sim];
	ui1[7] = un[sip];

	//------------------------------------------------------------

	uj1[0] = ut[sjm];		// fluid velocity of neighbor cells
	uj1[1] = ut[sjp];		// along y [j-1, j+1] stored in uj1
							// (ut, ux, uy, un)
	uj1[2] = ux[sjm];
	uj1[3] = ux[sjp];

	uj1[4] = uy[sjm];
	uj1[5] = uy[sjp];

	uj1[6] = un[sjm];
	uj1[7] = un[sjp];

	//------------------------------------------------------------

	uk1[0] = ut[skm];		// fluid velocity of neighbor cells
	uk1[1] = ut[skp];		// along z [k-1, k+1] stored in uk1
							// (ut, ux, uy, un)
	uk1[2] = ux[skm];
	uk1[3] = ux[skp];

	uk1[4] = uy[skm];
	uk1[5] = uy[skp];

	uk1[6] = un[skm];
	uk1[7] = un[skp];
}


void get_v_neighbor_cells(const precision * const __restrict__ ut, const precision * const __restrict__ ux, const precision * const __restrict__ uy, const precision * const __restrict__ un, precision * const __restrict__ vxi, precision * const __restrict__ vyj, precision * const __restrict__ vnk, int simm, int sim, int sip, int sipp, int sjmm, int sjm, int sjp, int sjpp, int skmm, int skm, int skp, int skpp)
{
	vxi[0]  = ux[simm] / ut[simm];		// vx neighbor cells [i-2, i-1, i+1, i+2] stored in vxi
	vxi[1]  = ux[sim]  / ut[sim];
	vxi[2]  = ux[sip]  / ut[sip];
	vxi[3]  = ux[sipp] / ut[sipp];

	//------------------------------------------------------------

	vyj[0]  = uy[sjmm] / ut[sjmm];		// vy neighbor cells [j-2, j-1, j+1, j+2] stored in vyj
	vyj[1]  = uy[sjm]  / ut[sjm];
	vyj[2]  = uy[sjp]  / ut[sjp];
	vyj[3]  = uy[sjpp] / ut[sjpp];

	//------------------------------------------------------------

	vnk[0]  = un[skmm] / ut[skmm];		// vn neighbor cells [k-2, k-1, k+1, k+2] stored in vnk
	vnk[1]  = un[skm]  / ut[skm];
	vnk[2]  = un[skp]  / ut[skp];
	vnk[3]  = un[skpp] / ut[skpp];
}


void get_q_neighbor_cells(const precision * const __restrict__ q_s, precision * const __restrict__ qi1, precision * const __restrict__ qj1, precision * const __restrict__ qk1, precision * const __restrict__ qi2, precision * const __restrict__ qj2, precision * const __restrict__ qk2, int * r, int simm, int sim, int sip, int sipp, int sjmm, int sjm, int sjp, int sjpp, int skmm, int skm, int skp, int skpp)
{
	int n = *r;

	qi1[n]     = q_s[sim];		// q neighbor cells [i-1, i+1] stored in qi1
	qi1[n + 1] = q_s[sip];

	qi2[n]     = q_s[simm];		// q neighbor cells [i-2, i+2] stored in qi2
	qi2[n + 1] = q_s[sipp];

	//------------------------------------------------------------

	qj1[n]     = q_s[sjm];		// q neighbor cells [j-1, j+1] stored in qj1
	qj1[n + 1] = q_s[sjp];

	qj2[n]     = q_s[sjmm];		// q neighbor cells [j-2, j+2] stored in qj2
	qj2[n + 1] = q_s[sjpp];

	//------------------------------------------------------------

	qk1[n]     = q_s[skm];		// q neighbor cells [k-1, k+1] stored in qk1
	qk1[n + 1] = q_s[skp];

	qk2[n]     = q_s[skmm];		// q neighbor cells [k-2, k+2] stored in qk2
	qk2[n + 1] = q_s[skpp];

	*r += 2;
}




