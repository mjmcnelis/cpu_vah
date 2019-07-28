#include <math.h>
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


void get_u_neighbor_cells(const precision * const __restrict__ ux, const precision * const __restrict__ uy, const precision * const __restrict__ un, precision * const __restrict__ ui1, precision * const __restrict__ uj1, precision * const __restrict__ uk1, int sim, int sip, int sjm, int sjp, int skm, int skp)
{
	ui1[0] = ux[sim];	// fluid velocity of neighbor cells along x [i-1, i+1] stored in ui1
	ui1[1] = ux[sip];
	ui1[2] = uy[sim];
	ui1[3] = uy[sip];
	ui1[4] = un[sim];
	ui1[5] = un[sip];

	uj1[0] = ux[sjm];	// fluid velocity of neighbor cells along y [j-1, j+1] stored in uj1
	uj1[1] = ux[sjp];
	uj1[2] = uy[sjm];
	uj1[3] = uy[sjp];
	uj1[4] = un[sjm];
	uj1[5] = un[sjp];

	uk1[0] = ux[skm];	// fluid velocity of neighbor cells along z [k-1, k+1] stored in uk1
	uk1[1] = ux[skp];
	uk1[2] = uy[skm];
	uk1[3] = uy[skp];
	uk1[4] = un[skm];
	uk1[5] = un[skp];
}


void get_v_neighbor_cells(const precision * const __restrict__ ux, const precision * const __restrict__ uy, const precision * const __restrict__ un, precision * const __restrict__ vxi, precision * const __restrict__ vyj, precision * const __restrict__ vnk, int simm, int sim, int sip, int sipp, int sjmm, int sjm, int sjp, int sjpp, int skmm, int skm, int skp, int skpp, precision t2)
{
	precision ux_simm = ux[simm];
	precision uy_simm = uy[simm];
	precision un_simm = un[simm];
	precision ut_simm = sqrt(1.  +  ux_simm * ux_simm  +  uy_simm * uy_simm  +  t2 * un_simm * un_simm);

	precision ux_sim = ux[sim];
	precision uy_sim = uy[sim];
	precision un_sim = un[sim];
	precision ut_sim = sqrt(1.  +  ux_sim * ux_sim  +  uy_sim * uy_sim  +  t2 * un_sim * un_sim);

	precision ux_sip = ux[sip];
	precision uy_sip = uy[sip];
	precision un_sip = un[sip];
	precision ut_sip = sqrt(1.  +  ux_sip * ux_sip  +  uy_sip * uy_sip  +  t2 * un_sip * un_sip);

	precision ux_sipp = ux[sipp];
	precision uy_sipp = uy[sipp];
	precision un_sipp = un[sipp];
	precision ut_sipp = sqrt(1.  +  ux_sipp * ux_sipp  +  uy_sipp * uy_sipp  +  t2 * un_sipp * un_sipp);

	vxi[0]  = ux_simm / ut_simm;	// vx neighbor cells [i-2, i-1, i+1, i+2] stored in vxi
	vxi[1]  = ux_sim  / ut_sim;
	vxi[2]  = ux_sip  / ut_sip;
	vxi[3]  = ux_sipp / ut_sipp;


	precision ux_sjmm = ux[sjmm];
	precision uy_sjmm = uy[sjmm];
	precision un_sjmm = un[sjmm];
	precision ut_sjmm = sqrt(1.  +  ux_sjmm * ux_sjmm  +  uy_sjmm * uy_sjmm  +  t2 * un_sjmm * un_sjmm);

	precision ux_sjm = ux[sjm];
	precision uy_sjm = uy[sjm];
	precision un_sjm = un[sjm];
	precision ut_sjm = sqrt(1.  +  ux_sjm * ux_sjm  +  uy_sjm * uy_sjm  +  t2 * un_sjm * un_sjm);

	precision ux_sjp = ux[sjp];
	precision uy_sjp = uy[sjp];
	precision un_sjp = un[sjp];
	precision ut_sjp = sqrt(1.  +  ux_sjp * ux_sjp  +  uy_sjp * uy_sjp  +  t2 * un_sjp * un_sjp);

	precision ux_sjpp = ux[sjpp];
	precision uy_sjpp = uy[sjpp];
	precision un_sjpp = un[sjpp];
	precision ut_sjpp = sqrt(1.  +  ux_sjpp * ux_sjpp  +  uy_sjpp * uy_sjpp  +  t2 * un_sjpp * un_sjpp);

	vyj[0]  = uy_sjmm / ut_sjmm;	// vy neighbor cells [j-2, j-1, j+1, j+2] stored in vyj
	vyj[1]  = uy_sjm  / ut_sjm;
	vyj[2]  = uy_sjp  / ut_sjp;
	vyj[3]  = uy_sjpp / ut_sjpp;


	precision ux_skmm = ux[skmm];
	precision uy_skmm = uy[skmm];
	precision un_skmm = un[skmm];
	precision ut_skmm = sqrt(1.  +  ux_skmm * ux_skmm  +  uy_skmm * uy_skmm  +  t2 * un_skmm * un_skmm);

	precision ux_skm = ux[skm];
	precision uy_skm = uy[skm];
	precision un_skm = un[skm];
	precision ut_skm = sqrt(1.  +  ux_skm * ux_skm  +  uy_skm * uy_skm  +  t2 * un_skm * un_skm);

	precision ux_skp = ux[skp];
	precision uy_skp = uy[skp];
	precision un_skp = un[skp];
	precision ut_skp = sqrt(1.  +  ux_skp * ux_skp  +  uy_skp * uy_skp  +  t2 * un_skp * un_skp);

	precision ux_skpp = ux[skpp];
	precision uy_skpp = uy[skpp];
	precision un_skpp = un[skpp];
	precision ut_skpp = sqrt(1.  +  ux_skpp * ux_skpp  +  uy_skpp * uy_skpp  +  t2 * un_skpp * un_skpp);

	vnk[0]  = un_skmm / ut_skmm;	// vn neighbor cells [k-2, k-1, k+1, k+2] stored in vnk
	vnk[1]  = un_skm  / ut_skm;
	vnk[2]  = un_skp  / ut_skp;
	vnk[3]  = un_skpp / ut_skpp;
}


void get_q_neighbor_cells(const precision * const __restrict__ q, precision * const __restrict__ qi1, precision * const __restrict__ qj1, precision * const __restrict__ qk1, precision * const __restrict__ qi2, precision * const __restrict__ qj2, precision * const __restrict__ qk2, int * r, int simm, int sim, int sip, int sipp, int sjmm, int sjm, int sjp, int sjpp, int skmm, int skm, int skp, int skpp)
{
	int n = *r;

	qi1[n]     = q[sim];		// q neighbor cells [i-1, i+1] stored in qi1
	qi1[n + 1] = q[sip];

	qi2[n]     = q[simm];		// q neighbor cells [i-2, i+2] stored in qi2
	qi2[n + 1] = q[sipp];

	//------------------------------------------------------------

	qj1[n]     = q[sjm];		// q neighbor cells [j-1, j+1] stored in qj1
	qj1[n + 1] = q[sjp];

	qj2[n]     = q[sjmm];		// q neighbor cells [j-2, j+2] stored in qj2
	qj2[n + 1] = q[sjpp];

	//------------------------------------------------------------

	qk1[n]     = q[skm];		// q neighbor cells [k-1, k+1] stored in qk1
	qk1[n + 1] = q[skp];

	qk2[n]     = q[skmm];		// q neighbor cells [k-2, k+2] stored in qk2
	qk2[n + 1] = q[skpp];

	*r += 2;
}




