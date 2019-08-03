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


void get_hydro_neighbor_cells(	hydro_variables q_simm, hydro_variables q_sim, hydro_variables q_sip, hydro_variables q_sipp,
								hydro_variables q_sjmm, hydro_variables q_sjm, hydro_variables q_sjp, hydro_variables q_sjpp,
								hydro_variables q_skmm, hydro_variables q_skm, hydro_variables q_skp, hydro_variables q_skpp,
								precision * const __restrict__ qi1, precision * const __restrict__ qj1, precision * const __restrict__ qk1,
								precision * const __restrict__ qi2, precision * const __restrict__ qj2, precision * const __restrict__ qk2)
{
	int n = 0;


	// q neighbors [i-1, i+1] stored in qi1
	//------------------------------------------------

	qi1[n]	   = q_sim.ttt;
	qi1[n + 1] = q_sip.ttt;		n += 2;

	qi1[n]	   = q_sim.ttx;
	qi1[n + 1] = q_sip.ttx;		n += 2;

	qi1[n]	   = q_sim.tty;
	qi1[n + 1] = q_sip.tty;		n += 2;

	qi1[n]	   = q_sim.ttn;
	qi1[n + 1] = q_sip.ttn;		n += 2;

#ifdef ANISO_HYDRO
	qi1[n]	   = q_sim.pl;
	qi1[n + 1] = q_sip.pl;		n += 2;

#if (PT_MATCHING == 1)
	qi1[n]	   = q_sim.pt;
	qi1[n + 1] = q_sip.pt;		n += 2;
#endif
#endif

#ifdef PIMUNU
	qi1[n]	   = q_sim.pitt;
	qi1[n + 1] = q_sip.pitt;	n += 2;

	qi1[n]	   = q_sim.pitx;
	qi1[n + 1] = q_sip.pitx;	n += 2;

	qi1[n]	   = q_sim.pity;
	qi1[n + 1] = q_sip.pity;	n += 2;

	qi1[n]	   = q_sim.pitn;
	qi1[n + 1] = q_sip.pitn;	n += 2;

	qi1[n]	   = q_sim.pixx;
	qi1[n + 1] = q_sip.pixx;	n += 2;

	qi1[n]	   = q_sim.pixy;
	qi1[n + 1] = q_sip.pixy;	n += 2;

	qi1[n]	   = q_sim.pixn;
	qi1[n + 1] = q_sip.pixn;	n += 2;

	qi1[n]	   = q_sim.piyy;
	qi1[n + 1] = q_sip.piyy;	n += 2;

	qi1[n]	   = q_sim.piyn;
	qi1[n + 1] = q_sip.piyn;	n += 2;

	qi1[n]	   = q_sim.pinn;
	qi1[n + 1] = q_sip.pinn;	n += 2;
#endif

#ifdef WTZMU
	qi1[n]	   = q_sim.WtTz;
	qi1[n + 1] = q_sip.WtTz;	n += 2;

	qi1[n]	   = q_sim.WxTz;
	qi1[n + 1] = q_sip.WxTz;	n += 2;

	qi1[n]	   = q_sim.WyTz;
	qi1[n + 1] = q_sip.WyTz;	n += 2;

	qi1[n]	   = q_sim.WnTz;
	qi1[n + 1] = q_sip.WnTz;
#endif

#ifdef PI
	qi1[n]	   = q_sim.Pi;
	qi1[n + 1] = q_sip.Pi;
#endif


	// q neighbors [i-2, i+2] stored in qi1
	//------------------------------------------------

	n = 0;

	qi2[n]	   = q_simm.ttt;
	qi2[n + 1] = q_sipp.ttt;	n += 2;

	qi2[n]	   = q_simm.ttx;
	qi2[n + 1] = q_sipp.ttx;	n += 2;

	qi2[n]	   = q_simm.tty;
	qi2[n + 1] = q_sipp.tty;	n += 2;

	qi2[n]	   = q_simm.ttn;
	qi2[n + 1] = q_sipp.ttn;	n += 2;

#ifdef ANISO_HYDRO
	qi2[n]	   = q_simm.pl;
	qi2[n + 1] = q_sipp.pl;		n += 2;

#if (PT_MATCHING == 1)
	qi2[n]	   = q_simm.pt;
	qi2[n + 1] = q_sipp.pt;		n += 2;
#endif
#endif

#ifdef PIMUNU
	qi2[n]	   = q_simm.pitt;
	qi2[n + 1] = q_sipp.pitt;	n += 2;

	qi2[n]	   = q_simm.pitx;
	qi2[n + 1] = q_sipp.pitx;	n += 2;

	qi2[n]	   = q_simm.pity;
	qi2[n + 1] = q_sipp.pity;	n += 2;

	qi2[n]	   = q_simm.pitn;
	qi2[n + 1] = q_sipp.pitn;	n += 2;

	qi2[n]	   = q_simm.pixx;
	qi2[n + 1] = q_sipp.pixx;	n += 2;

	qi2[n]	   = q_simm.pixy;
	qi2[n + 1] = q_sipp.pixy;	n += 2;

	qi2[n]	   = q_simm.pixn;
	qi2[n + 1] = q_sipp.pixn;	n += 2;

	qi2[n]	   = q_simm.piyy;
	qi2[n + 1] = q_sipp.piyy;	n += 2;

	qi2[n]	   = q_simm.piyn;
	qi2[n + 1] = q_sipp.piyn;	n += 2;

	qi2[n]	   = q_simm.pinn;
	qi2[n + 1] = q_sipp.pinn;	n += 2;
#endif

#ifdef WTZMU
	qi2[n]	   = q_simm.WtTz;
	qi2[n + 1] = q_sipp.WtTz;	n += 2;

	qi2[n]	   = q_simm.WxTz;
	qi2[n + 1] = q_sipp.WxTz;	n += 2;

	qi2[n]	   = q_simm.WyTz;
	qi2[n + 1] = q_sipp.WyTz;	n += 2;

	qi2[n]	   = q_simm.WnTz;
	qi2[n + 1] = q_sipp.WnTz;
#endif

#ifdef PI
	qi2[n]	   = q_simm.Pi;
	qi2[n + 1] = q_sipp.Pi;
#endif


	// q neighbors [j-1, j+1] stored in qj1
	//------------------------------------------------

	n = 0;

	qj1[n]	   = q_sjm.ttt;
	qj1[n + 1] = q_sjp.ttt;		n += 2;

	qj1[n]	   = q_sjm.ttx;
	qj1[n + 1] = q_sjp.ttx;		n += 2;

	qj1[n]	   = q_sjm.tty;
	qj1[n + 1] = q_sjp.tty;		n += 2;

	qj1[n]	   = q_sjm.ttn;
	qj1[n + 1] = q_sjp.ttn;		n += 2;

#ifdef ANISO_HYDRO
	qj1[n]	   = q_sjm.pl;
	qj1[n + 1] = q_sjp.pl;		n += 2;

#if (PT_MATCHING == 1)
	qj1[n]	   = q_sjm.pt;
	qj1[n + 1] = q_sjp.pt;		n += 2;
#endif
#endif

#ifdef PIMUNU
	qj1[n]	   = q_sjm.pitt;
	qj1[n + 1] = q_sjp.pitt;	n += 2;

	qj1[n]	   = q_sjm.pitx;
	qj1[n + 1] = q_sjp.pitx;	n += 2;

	qj1[n]	   = q_sjm.pity;
	qj1[n + 1] = q_sjp.pity;	n += 2;

	qj1[n]	   = q_sjm.pitn;
	qj1[n + 1] = q_sjp.pitn;	n += 2;

	qj1[n]	   = q_sjm.pixx;
	qj1[n + 1] = q_sjp.pixx;	n += 2;

	qj1[n]	   = q_sjm.pixy;
	qj1[n + 1] = q_sjp.pixy;	n += 2;

	qj1[n]	   = q_sjm.pixn;
	qj1[n + 1] = q_sjp.pixn;	n += 2;

	qj1[n]	   = q_sjm.piyy;
	qj1[n + 1] = q_sjp.piyy;	n += 2;

	qj1[n]	   = q_sjm.piyn;
	qj1[n + 1] = q_sjp.piyn;	n += 2;

	qj1[n]	   = q_sjm.pinn;
	qj1[n + 1] = q_sjp.pinn;	n += 2;
#endif

#ifdef WTZMU
	qj1[n]	   = q_sjm.WtTz;
	qj1[n + 1] = q_sjp.WtTz;	n += 2;

	qj1[n]	   = q_sjm.WxTz;
	qj1[n + 1] = q_sjp.WxTz;	n += 2;

	qj1[n]	   = q_sjm.WyTz;
	qj1[n + 1] = q_sjp.WyTz;	n += 2;

	qj1[n]	   = q_sjm.WnTz;
	qj1[n + 1] = q_sjp.WnTz;
#endif

#ifdef PI
	qj1[n]	   = q_sjm.Pi;
	qj1[n + 1] = q_sjp.Pi;
#endif


	// q neighbors [j-2, j+2] stored in qj2
	//------------------------------------------------

	n = 0;

	qj2[n]	   = q_sjmm.ttt;
	qj2[n + 1] = q_sjpp.ttt;	n += 2;

	qj2[n]	   = q_sjmm.ttx;
	qj2[n + 1] = q_sjpp.ttx;	n += 2;

	qj2[n]	   = q_sjmm.tty;
	qj2[n + 1] = q_sjpp.tty;	n += 2;

	qj2[n]	   = q_sjmm.ttn;
	qj2[n + 1] = q_sjpp.ttn;	n += 2;

#ifdef ANISO_HYDRO
	qj2[n]	   = q_sjmm.pl;
	qj2[n + 1] = q_sjpp.pl;		n += 2;

#if (PT_MATCHING == 1)
	qj2[n]	   = q_sjmm.pt;
	qj2[n + 1] = q_sjpp.pt;		n += 2;
#endif
#endif

#ifdef PIMUNU
	qj2[n]	   = q_sjmm.pitt;
	qj2[n + 1] = q_sjpp.pitt;	n += 2;

	qj2[n]	   = q_sjmm.pitx;
	qj2[n + 1] = q_sjpp.pitx;	n += 2;

	qj2[n]	   = q_sjmm.pity;
	qj2[n + 1] = q_sjpp.pity;	n += 2;

	qj2[n]	   = q_sjmm.pitn;
	qj2[n + 1] = q_sjpp.pitn;	n += 2;

	qj2[n]	   = q_sjmm.pixx;
	qj2[n + 1] = q_sjpp.pixx;	n += 2;

	qj2[n]	   = q_sjmm.pixy;
	qj2[n + 1] = q_sjpp.pixy;	n += 2;

	qj2[n]	   = q_sjmm.pixn;
	qj2[n + 1] = q_sjpp.pixn;	n += 2;

	qj2[n]	   = q_sjmm.piyy;
	qj2[n + 1] = q_sjpp.piyy;	n += 2;

	qj2[n]	   = q_sjmm.piyn;
	qj2[n + 1] = q_sjpp.piyn;	n += 2;

	qj2[n]	   = q_sjmm.pinn;
	qj2[n + 1] = q_sjpp.pinn;	n += 2;
#endif

#ifdef WTZMU
	qj2[n]	   = q_sjmm.WtTz;
	qj2[n + 1] = q_sjpp.WtTz;	n += 2;

	qj2[n]	   = q_sjmm.WxTz;
	qj2[n + 1] = q_sjpp.WxTz;	n += 2;

	qj2[n]	   = q_sjmm.WyTz;
	qj2[n + 1] = q_sjpp.WyTz;	n += 2;

	qj2[n]	   = q_sjmm.WnTz;
	qj2[n + 1] = q_sjpp.WnTz;
#endif

#ifdef PI
	qj2[n]	   = q_sjmm.Pi;
	qj2[n + 1] = q_sjpp.Pi;
#endif


	// q neighbors [k-1, k+1] stored in qk1
	//------------------------------------------------

	n = 0;

	qk1[n]	   = q_skm.ttt;
	qk1[n + 1] = q_skp.ttt;		n += 2;

	qk1[n]	   = q_skm.ttx;
	qk1[n + 1] = q_skp.ttx;		n += 2;

	qk1[n]	   = q_skm.tty;
	qk1[n + 1] = q_skp.tty;		n += 2;

	qk1[n]	   = q_skm.ttn;
	qk1[n + 1] = q_skp.ttn;		n += 2;

#ifdef ANISO_HYDRO
	qk1[n]	   = q_skm.pl;
	qk1[n + 1] = q_skp.pl;		n += 2;

#if (PT_MATCHING == 1)
	qk1[n]	   = q_skm.pt;
	qk1[n + 1] = q_skp.pt;		n += 2;
#endif
#endif

#ifdef PIMUNU
	qk1[n]	   = q_skm.pitt;
	qk1[n + 1] = q_skp.pitt;	n += 2;

	qk1[n]	   = q_skm.pitx;
	qk1[n + 1] = q_skp.pitx;	n += 2;

	qk1[n]	   = q_skm.pity;
	qk1[n + 1] = q_skp.pity;	n += 2;

	qk1[n]	   = q_skm.pitn;
	qk1[n + 1] = q_skp.pitn;	n += 2;

	qk1[n]	   = q_skm.pixx;
	qk1[n + 1] = q_skp.pixx;	n += 2;

	qk1[n]	   = q_skm.pixy;
	qk1[n + 1] = q_skp.pixy;	n += 2;

	qk1[n]	   = q_skm.pixn;
	qk1[n + 1] = q_skp.pixn;	n += 2;

	qk1[n]	   = q_skm.piyy;
	qk1[n + 1] = q_skp.piyy;	n += 2;

	qk1[n]	   = q_skm.piyn;
	qk1[n + 1] = q_skp.piyn;	n += 2;

	qk1[n]	   = q_skm.pinn;
	qk1[n + 1] = q_skp.pinn;	n += 2;
#endif

#ifdef WTZMU
	qk1[n]	   = q_skm.WtTz;
	qk1[n + 1] = q_skp.WtTz;	n += 2;

	qk1[n]	   = q_skm.WxTz;
	qk1[n + 1] = q_skp.WxTz;	n += 2;

	qk1[n]	   = q_skm.WyTz;
	qk1[n + 1] = q_skp.WyTz;	n += 2;

	qk1[n]	   = q_skm.WnTz;
	qk1[n + 1] = q_skp.WnTz;
#endif

#ifdef PI
	qk1[n]	   = q_skm.Pi;
	qk1[n + 1] = q_skp.Pi;
#endif


	// q neighbors [k-2, k+2] stored in qk2
	//------------------------------------------------

	n = 0;

	qk2[n]	   = q_skmm.ttt;
	qk2[n + 1] = q_skpp.ttt;	n += 2;

	qk2[n]	   = q_skmm.ttx;
	qk2[n + 1] = q_skpp.ttx;	n += 2;

	qk2[n]	   = q_skmm.tty;
	qk2[n + 1] = q_skpp.tty;	n += 2;

	qk2[n]	   = q_skmm.ttn;
	qk2[n + 1] = q_skpp.ttn;	n += 2;

#ifdef ANISO_HYDRO
	qk2[n]	   = q_skmm.pl;
	qk2[n + 1] = q_skpp.pl;		n += 2;

#if (PT_MATCHING == 1)
	qk2[n]	   = q_skmm.pt;
	qk2[n + 1] = q_skpp.pt;		n += 2;
#endif
#endif

#ifdef PIMUNU
	qk2[n]	   = q_skmm.pitt;
	qk2[n + 1] = q_skpp.pitt;	n += 2;

	qk2[n]	   = q_skmm.pitx;
	qk2[n + 1] = q_skpp.pitx;	n += 2;

	qk2[n]	   = q_skmm.pity;
	qk2[n + 1] = q_skpp.pity;	n += 2;

	qk2[n]	   = q_skmm.pitn;
	qk2[n + 1] = q_skpp.pitn;	n += 2;

	qk2[n]	   = q_skmm.pixx;
	qk2[n + 1] = q_skpp.pixx;	n += 2;

	qk2[n]	   = q_skmm.pixy;
	qk2[n + 1] = q_skpp.pixy;	n += 2;

	qk2[n]	   = q_skmm.pixn;
	qk2[n + 1] = q_skpp.pixn;	n += 2;

	qk2[n]	   = q_skmm.piyy;
	qk2[n + 1] = q_skpp.piyy;	n += 2;

	qk2[n]	   = q_skmm.piyn;
	qk2[n + 1] = q_skpp.piyn;	n += 2;

	qk2[n]	   = q_skmm.pinn;
	qk2[n + 1] = q_skpp.pinn;	n += 2;
#endif

#ifdef WTZMU
	qk2[n]	   = q_skmm.WtTz;
	qk2[n + 1] = q_skpp.WtTz;	n += 2;

	qk2[n]	   = q_skmm.WxTz;
	qk2[n + 1] = q_skpp.WxTz;	n += 2;

	qk2[n]	   = q_skmm.WyTz;
	qk2[n + 1] = q_skpp.WyTz;	n += 2;

	qk2[n]	   = q_skmm.WnTz;
	qk2[n + 1] = q_skpp.WnTz;
#endif

#ifdef PI
	qk2[n]	   = q_skmm.Pi;
	qk2[n + 1] = q_skpp.Pi;
#endif
}




