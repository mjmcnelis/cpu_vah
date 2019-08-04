#include <math.h>
#include "../include/DynamicalVariables.h"
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


void get_fluid_velocity_neighbor_cells(	fluid_velocity u_simm, fluid_velocity u_sim, fluid_velocity u_sip, fluid_velocity u_sipp,
										fluid_velocity u_sjmm, fluid_velocity u_sjm, fluid_velocity u_sjp, fluid_velocity u_sjpp,
										fluid_velocity u_skmm, fluid_velocity u_skm, fluid_velocity u_skp, fluid_velocity u_skpp,
										precision * const __restrict__ ui1, precision * const __restrict__ uj1, precision * const __restrict__ uk1,
										precision * const __restrict__ vxi, precision * const __restrict__ vyj, precision * const __restrict__ vnk, precision t2)
{
	precision ux_simm = u_simm.ux;
	precision uy_simm = u_simm.uy;
	precision un_simm = u_simm.un;
	precision ut_simm = sqrt(1.  +  ux_simm * ux_simm  +  uy_simm * uy_simm  +  t2 * un_simm * un_simm);

	precision ux_sim = u_sim.ux;
	precision uy_sim = u_sim.uy;
	precision un_sim = u_sim.un;
	precision ut_sim = sqrt(1.  +  ux_sim * ux_sim  +  uy_sim * uy_sim  +  t2 * un_sim * un_sim);

	precision ux_sip = u_sip.ux;
	precision uy_sip = u_sip.uy;
	precision un_sip = u_sip.un;
	precision ut_sip = sqrt(1.  +  ux_sip * ux_sip  +  uy_sip * uy_sip  +  t2 * un_sip * un_sip);

	precision ux_sipp = u_sipp.ux;
	precision uy_sipp = u_sipp.uy;
	precision un_sipp = u_sipp.un;
	precision ut_sipp = sqrt(1.  +  ux_sipp * ux_sipp  +  uy_sipp * uy_sipp  +  t2 * un_sipp * un_sipp);

	vxi[0] = ux_simm / ut_simm;		// vx neighbors [i-2, i-1, i+1, i+2] stored in vxi
	vxi[1] = ux_sim  / ut_sim;
	vxi[2] = ux_sip  / ut_sip;
	vxi[3] = ux_sipp / ut_sipp;

	ui1[0] = ux_sim;				// u neighbors [i-1, i+1] stored in ui1
	ui1[1] = ux_sip;
	ui1[2] = uy_sim;
	ui1[3] = uy_sip;
	ui1[4] = un_sim;
	ui1[5] = un_sip;

	//------------------------------------------------

	precision ux_sjmm = u_sjmm.ux;
	precision uy_sjmm = u_sjmm.uy;
	precision un_sjmm = u_sjmm.un;
	precision ut_sjmm = sqrt(1.  +  ux_sjmm * ux_sjmm  +  uy_sjmm * uy_sjmm  +  t2 * un_sjmm * un_sjmm);

	precision ux_sjm = u_sjm.ux;
	precision uy_sjm = u_sjm.uy;
	precision un_sjm = u_sjm.un;
	precision ut_sjm = sqrt(1.  +  ux_sjm * ux_sjm  +  uy_sjm * uy_sjm  +  t2 * un_sjm * un_sjm);

	precision ux_sjp = u_sjp.ux;
	precision uy_sjp = u_sjp.uy;
	precision un_sjp = u_sjp.un;
	precision ut_sjp = sqrt(1.  +  ux_sjp * ux_sjp  +  uy_sjp * uy_sjp  +  t2 * un_sjp * un_sjp);

	precision ux_sjpp = u_sjpp.ux;
	precision uy_sjpp = u_sjpp.uy;
	precision un_sjpp = u_sjpp.un;
	precision ut_sjpp = sqrt(1.  +  ux_sjpp * ux_sjpp  +  uy_sjpp * uy_sjpp  +  t2 * un_sjpp * un_sjpp);

	vyj[0] = uy_sjmm / ut_sjmm;		// vy neighbors [j-2, j-1, j+1, j+2] stored in vyj
	vyj[1] = uy_sjm  / ut_sjm;
	vyj[2] = uy_sjp  / ut_sjp;
	vyj[3] = uy_sjpp / ut_sjpp;

	uj1[0] = ux_sjm;				// u neighbors [j-1, j+1] stored in uj1
	uj1[1] = ux_sjp;
	uj1[2] = uy_sjm;
	uj1[3] = uy_sjp;
	uj1[4] = un_sjm;
	uj1[5] = un_sjp;

	//------------------------------------------------

	precision ux_skmm = u_skmm.ux;
	precision uy_skmm = u_skmm.uy;
	precision un_skmm = u_skmm.un;
	precision ut_skmm = sqrt(1.  +  ux_skmm * ux_skmm  +  uy_skmm * uy_skmm  +  t2 * un_skmm * un_skmm);

	precision ux_skm = u_skm.ux;
	precision uy_skm = u_skm.uy;
	precision un_skm = u_skm.un;
	precision ut_skm = sqrt(1.  +  ux_skm * ux_skm  +  uy_skm * uy_skm  +  t2 * un_skm * un_skm);

	precision ux_skp = u_skp.ux;
	precision uy_skp = u_skp.uy;
	precision un_skp = u_skp.un;
	precision ut_skp = sqrt(1.  +  ux_skp * ux_skp  +  uy_skp * uy_skp  +  t2 * un_skp * un_skp);

	precision ux_skpp = u_skpp.ux;
	precision uy_skpp = u_skpp.uy;
	precision un_skpp = u_skpp.un;
	precision ut_skpp = sqrt(1.  +  ux_skpp * ux_skpp  +  uy_skpp * uy_skpp  +  t2 * un_skpp * un_skpp);

	vnk[0] = un_skmm / ut_skmm;		// vn neighbors [k-2, k-1, k+1, k+2] stored in vnk
	vnk[1] = un_skm  / ut_skm;
	vnk[2] = un_skp  / ut_skp;
	vnk[3] = un_skpp / ut_skpp;

	uk1[0] = ux_skm;				// u neighbors [k-1, k+1] stored in uk1
	uk1[1] = ux_skp;
	uk1[2] = uy_skm;
	uk1[3] = uy_skp;
	uk1[4] = un_skm;
	uk1[5] = un_skp;
}


void get_hydro_neighbor_cells(hydro_variables qm, hydro_variables qp, precision * const __restrict__ q)
{
	int n = 0;

	// q neighbors [qm, qp] stored in q
	//------------------------------------------------

	q[n]	 = qm.ttt;
	q[n + 1] = qp.ttt;	n += 2;

	q[n]	 = qm.ttx;
	q[n + 1] = qp.ttx;	n += 2;

	q[n]	 = qm.tty;
	q[n + 1] = qp.tty;	n += 2;

	q[n]	 = qm.ttn;
	q[n + 1] = qp.ttn;	n += 2;

#ifdef ANISO_HYDRO
	q[n]	 = qm.pl;
	q[n + 1] = qp.pl;	n += 2;

#if (PT_MATCHING == 1)
	q[n]	 = qm.pt;
	q[n + 1] = qp.pt;	n += 2;
#endif
#endif

#ifdef PIMUNU
	q[n]	 = qm.pitt;
	q[n + 1] = qp.pitt;	n += 2;

	q[n]	 = qm.pitx;
	q[n + 1] = qp.pitx;	n += 2;

	q[n]	 = qm.pity;
	q[n + 1] = qp.pity;	n += 2;

	q[n]	 = qm.pitn;
	q[n + 1] = qp.pitn;	n += 2;

	q[n]	 = qm.pixx;
	q[n + 1] = qp.pixx;	n += 2;

	q[n]	 = qm.pixy;
	q[n + 1] = qp.pixy;	n += 2;

	q[n]	 = qm.pixn;
	q[n + 1] = qp.pixn;	n += 2;

	q[n]	 = qm.piyy;
	q[n + 1] = qp.piyy;	n += 2;

	q[n]     = qm.piyn;
	q[n + 1] = qp.piyn;	n += 2;

	q[n]	 = qm.pinn;
	q[n + 1] = qp.pinn;	n += 2;
#endif

#ifdef WTZMU
	q[n]  	 = qm.WtTz;
	q[n + 1] = qp.WtTz;	n += 2;

	q[n]	 = qm.WxTz;
	q[n + 1] = qp.WxTz;	n += 2;

	q[n]	 = qm.WyTz;
	q[n + 1] = qp.WyTz;	n += 2;

	q[n]	 = qm.WnTz;
	q[n + 1] = qp.WnTz;
#endif

#ifdef PI
	q[n]	 = qm.Pi;
	q[n + 1] = qp.Pi;
#endif
}









