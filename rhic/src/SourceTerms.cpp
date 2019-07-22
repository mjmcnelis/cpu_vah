#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include "../include/SourceTerms.h"
#include "../include/FluxTerms.h"
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"
#include "../include/EquationOfState.h" // for bulk terms
#include "../include/TransportCoefficients.h"
#include "../include/Projections.h"

using namespace std;

#define A_1 -13.77
#define A_2 27.55
#define A_3 13.45

#define LAMBDA_1 0.9
#define LAMBDA_2 0.25
#define LAMBDA_3 0.9
#define LAMBDA_4 0.22

#define SIGMA_1 0.025
#define SIGMA_2 0.13
#define SIGMA_3 0.0025
#define SIGMA_4 0.022


inline precision bulkViscosityToEntropyDensity(precision T)   // need some hydro parameters here
{
	precision x = T / 1.01355;

	if(x > 1.05)
	{
		return LAMBDA_1*exp(-(x-1)/SIGMA_1) + LAMBDA_2*exp(-(x-1)/SIGMA_2)+0.001;
	}
	else if(x < 0.995)
	{
		return LAMBDA_3*exp((x-1)/SIGMA_3)+ LAMBDA_4*exp((x-1)/SIGMA_4)+0.03;
	}
	else
	{
		return A_1*x*x + A_2*x - A_3;
	}
}

// inline precision sign(precision x)
// {
// 	if(x < 0.0) return -1.0;
// 	else return 1.0;
// }

// inline precision minmod(precision x, precision y)
// {
// 	return (sign(x) + sign(y)) * fmin(fabs(x), fabs(y)) / 2.0;
// }

// inline precision minmod3(precision x, precision y, precision z)
// {
//    return minmod(x, minmod(y, z));
// }

inline precision spatial_derivative(precision qm, precision q, precision qp)
{
	//return minmod3(THETA * (q - qm), (qp - qm) / 2.0, THETA * (qp - q));
	return (qp - qm) / 2.0;
}



void source_terms(precision * const __restrict__ S, const precision * const __restrict__ q, precision e, precision t, const precision * const __restrict__ qi1, const precision * const __restrict__ qj1, const precision * const __restrict__ qk1, const precision * const __restrict__ e1, const precision * const __restrict__ ui1, const precision * const __restrict__ uj1, const precision * const __restrict__ uk1, precision ut, precision ux, precision uy, precision un, precision ut_p, precision ux_p, precision uy_p, precision un_p, precision dt, precision dx, precision dy, precision dn, precision etabar)
{
	precision utperp  = sqrt(1.0  +  ux * ux  +  uy * uy);
	
	precision zt = t * un / utperp;		// longitudinal basis vector 
	precision zn = ut / t / utperp;

	precision t2 = t * t;				// useful expressions
	precision t4 = t2 * t2;

	precision zt2  = zt * zt;
	precision ztzn = zt * zn;
	precision zn2  = zn * zn;
	precision t2zn = t2 * zn;

	precision ut2  = ut * ut;
	precision un2  = un * un;
	precision utun = ut * un;

	precision unzn = un * zn;

	precision utperp2 = utperp * utperp;

	precision p = equilibriumPressure(e);

	// conserved variables
	precision ttt = q[0];
	precision ttx = q[1];
	precision tty = q[2];
	precision ttn = q[3];
	precision pl  = q[4];

	int a = 5;

#if (PT_MATCHING == 1)
	precision pt  = q[a];	a++;
#else
	precision pt  = 0.5 * (e - pl);
#endif
#ifdef PIMUNU
	precision pitt = q[a];	a++;
	precision pitx = q[a];	a++;
	precision pity = q[a];	a++;
	precision pitn = q[a];	a++;
	precision pixx = q[a];	a++;
	precision pixy = q[a];	a++;
	precision pixn = q[a];	a++;
	precision piyy = q[a];	a++;
	precision piyn = q[a];	a++;
	precision pinn = q[a];	a++;
#else
	precision pitt = 0.0;
	precision pitx = 0.0;
	precision pity = 0.0;
	precision pitn = 0.0;
	precision pixx = 0.0;
	precision pixy = 0.0;
	precision pixn = 0.0;
	precision piyy = 0.0;
	precision piyn = 0.0;
	precision pinn = 0.0;
#endif
#ifdef WTZMU
	precision WtTz = q[a];	a++;
	precision WxTz = q[a];	a++;
	precision WyTz = q[a];	a++;
	precision WnTz = q[a];
#else
	precision WtTz = 0.0;
	precision WxTz = 0.0;
	precision WyTz = 0.0;
	precision WnTz = 0.0;
#endif


	// primary variable derivatives
#if (PT_MATCHING == 0)
	precision e_sim = e1[0];	// e [i-1, i+1]
	precision e_sip = e1[1];

	precision e_sjm = e1[2];	// e [j-1, j+1]
	precision e_sjp = e1[3];
	
	precision e_skm = e1[4];	// e [k-1, k+1]
	precision e_skp = e1[5];

	precision de_dx = spatial_derivative(e_sim, e, e_sip) / dx;
	precision de_dy = spatial_derivative(e_sjm, e, e_sjp) / dy;
	precision de_dn = spatial_derivative(e_skm, e, e_skp) / dn;
#endif

	// longitudinal pressure derivatives
	int n = 8;

	precision pl_sim = qi1[n];		// pl [i-1, i+1]	(sim = i-1)
	precision pl_sip = qi1[n + 1];	//				  	(sip = i+1)

	precision pl_sjm = qj1[n];		// pl [j-1, j+1]	(sjm = j-1)
	precision pl_sjp = qj1[n + 1];  //				 	(sjp = j+1)

	precision pl_skm = qk1[n];		// pl [k-1, k+1]  	(skm = k-1)
	precision pl_skp = qk1[n + 1];	//					(skp = k+1)

	precision dpl_dx = spatial_derivative(pl_sim, pl, pl_sip) / dx;
	precision dpl_dy = spatial_derivative(pl_sjm, pl, pl_sjp) / dy;
	precision dpl_dn = spatial_derivative(pl_skm, pl, pl_skp) / dn;

	// transverse pressure derivatives
#if (PT_MATCHING == 1)
	n += 2;

	precision pt_sim = qi1[n];		// pt [i-1, i+1]
	precision pt_sip = qi1[n + 1];

	precision pt_sjm = qj1[n];		// pt [j-1, j+1]
	precision pt_sjp = qj1[n + 1];

	precision pt_skm = qk1[n];		// pt [k-1, k+1]
	precision pt_skp = qk1[n + 1];

	precision dpt_dx = spatial_derivative(pt_sim, pt, pt_sip) / dx;
	precision dpt_dy = spatial_derivative(pt_sjm, pt, pt_sjp) / dy;
	precision dpt_dn = spatial_derivative(pt_skm, pt, pt_skp) / dn;
#else
	precision dpt_dx = 0.5 * (de_dx - dpl_dx);
	precision dpt_dy = 0.5 * (de_dy - dpl_dy);
	precision dpt_dn = 0.5 * (de_dn - dpl_dn);
#endif
#ifdef PIMUNU
	n += 2;

	precision pitt_sim = qi1[n];		// pitt [i-1, i+1]
	precision pitt_sip = qi1[n + 1];

	precision pitt_sjm = qj1[n];		// pitt [j-1, j+1]
	precision pitt_sjp = qj1[n + 1];

	precision pitt_skm = qk1[n];		// pitt [k-1, k+1]
	precision pitt_skp = qk1[n + 1];

	precision dpitt_dx = spatial_derivative(pitt_sim, pitt, pitt_sip) / dx;
	precision dpitt_dy = spatial_derivative(pitt_sjm, pitt, pitt_sjp) / dy;
	precision dpitt_dn = spatial_derivative(pitt_skm, pitt, pitt_skp) / dn;

	n += 2;

	precision pitx_sim = qi1[n];		// pitx [i-1, i+1]
	precision pitx_sip = qi1[n + 1];

	precision pitx_sjm = qj1[n];		// pitx [j-1, j+1]
	precision pitx_sjp = qj1[n + 1];

	precision pitx_skm = qk1[n];		// pitx [k-1, k+1]
	precision pitx_skp = qk1[n + 1];

	precision dpitx_dx = spatial_derivative(pitx_sim, pitx, pitx_sip) / dx;
	precision dpitx_dy = spatial_derivative(pitx_sjm, pitx, pitx_sjp) / dy;
	precision dpitx_dn = spatial_derivative(pitx_skm, pitx, pitx_skp) / dn;

	n += 2;

	precision pity_sim = qi1[n];		// pity [i-1, i+1]
	precision pity_sip = qi1[n + 1];

	precision pity_sjm = qj1[n];		// pity [j-1, j+1]
	precision pity_sjp = qj1[n + 1];

	precision pity_skm = qk1[n];		// pity [k-1, k+1]
	precision pity_skp = qk1[n + 1];

	precision dpity_dx = spatial_derivative(pity_sim, pity, pity_sip) / dx;
	precision dpity_dy = spatial_derivative(pity_sjm, pity, pity_sjp) / dy;
	precision dpity_dn = spatial_derivative(pity_skm, pity, pity_skp) / dn;

	n += 2;

	precision pitn_sim = qi1[n];		// pitn [i-1, i+1]
	precision pitn_sip = qi1[n + 1];

	precision pitn_sjm = qj1[n];		// pitn [j-1, j+1]
	precision pitn_sjp = qj1[n + 1];

	precision pitn_skm = qk1[n];		// pitn [k-1, k+1]
	precision pitn_skp = qk1[n + 1];

	precision dpitn_dx = spatial_derivative(pitn_sim, pitn, pitn_sip) / dx;
	precision dpitn_dy = spatial_derivative(pitn_sjm, pitn, pitn_sjp) / dy;
	precision dpitn_dn = spatial_derivative(pitn_skm, pitn, pitn_skp) / dn;

	n += 2;

	precision pixx_sim = qi1[n];		// pixx [i-1, i+1]
	precision pixx_sip = qi1[n + 1];

	precision dpixx_dx = spatial_derivative(pixx_sim, pixx, pixx_sip) / dx;

	n += 2;

	precision pixy_sim = qi1[n];		// pixy [i-1, i+1]
	precision pixy_sip = qi1[n + 1];

	precision pixy_sjm = qj1[n];		// pixy [j-1, j+1]
	precision pixy_sjp = qj1[n + 1];

	precision dpixy_dx = spatial_derivative(pixy_sim, pixy, pixy_sip) / dx;
	precision dpixy_dy = spatial_derivative(pixy_sjm, pixy, pixy_sjp) / dy;

	n += 2;

	precision pixn_sim = qi1[n];		// pixn [i-1, i+1]
	precision pixn_sip = qi1[n + 1];

	precision pixn_skm = qk1[n];		// pixn [k-1, k+1]
	precision pixn_skp = qk1[n + 1];

	precision dpixn_dx = spatial_derivative(pixn_sim, pixn, pixn_sip) / dx;
	precision dpixn_dn = spatial_derivative(pixn_skm, pixn, pixn_skp) / dn;

	n += 2;

	precision piyy_sjm = qj1[n];		// piyy [j-1, j+1]
	precision piyy_sjp = qj1[n + 1];

	precision dpiyy_dy = spatial_derivative(piyy_sjm, piyy, piyy_sjp) / dy;

	n += 2;

	precision piyn_sjm = qj1[n];		// piyn [j-1, j+1]
	precision piyn_sjp = qj1[n + 1];

	precision piyn_skm = qk1[n];		// piyn [k-1, k+1]
	precision piyn_skp = qk1[n + 1];

	precision dpiyn_dy = spatial_derivative(piyn_sjm, piyn, piyn_sjp) / dy;
	precision dpiyn_dn = spatial_derivative(piyn_skm, piyn, piyn_skp) / dn;

	n += 2;

	precision pinn_skm = qk1[n];		// pinn [k-1, k+1]
	precision pinn_skp = qk1[n + 1];

	precision dpinn_dn = spatial_derivative(pinn_skm, pinn, pinn_skp) / dn;
#else
	precision dpitt_dx = 0.0;
	precision dpitt_dy = 0.0;
	precision dpitt_dn = 0.0;

	precision dpitx_dx = 0.0;
	precision dpitx_dy = 0.0;
	precision dpitx_dn = 0.0;

	precision dpity_dx = 0.0;
	precision dpity_dy = 0.0;
	precision dpity_dn = 0.0;

	precision dpitn_dx = 0.0;
	precision dpitn_dy = 0.0;
	precision dpitn_dn = 0.0;

	precision dpixx_dx = 0.0;

	precision dpixy_dx = 0.0;
	precision dpixy_dy = 0.0;

	precision dpixn_dx = 0.0;
	precision dpixn_dn = 0.0;

	precision dpiyy_dy = 0.0;

	precision dpiyn_dy = 0.0;
	precision dpiyn_dn = 0.0;

	precision dpinn_dn = 0.0;
#endif
#ifdef WTZMU
	n += 2;

	precision WtTz_sim = qi1[n];		// WtTz [i-1, i+1]
	precision WtTz_sip = qi1[n + 1];

	precision WtTz_sjm = qj1[n];		// WtTz [j-1, j+1]
	precision WtTz_sjp = qj1[n + 1];

	precision WtTz_skm = qk1[n];		// WtTz [k-1, k+1]
	precision WtTz_skp = qk1[n + 1];

	precision dWtTz_dx = spatial_derivative(WtTz_sim, WtTz, WtTz_sip) / dx;
	precision dWtTz_dy = spatial_derivative(WtTz_sjm, WtTz, WtTz_sjp) / dy;
	precision dWtTz_dn = spatial_derivative(WtTz_skm, WtTz, WtTz_skp) / dn;

	n += 2;

	precision WxTz_sim = qi1[n];		// WxTz [i-1, i+1]
	precision WxTz_sip = qi1[n + 1];

	precision WxTz_sjm = qj1[n];		// WxTz [j-1, j+1]
	precision WxTz_sjp = qj1[n + 1];

	precision WxTz_skm = qk1[n];		// WxTz [k-1, k+1]
	precision WxTz_skp = qk1[n + 1];

	precision dWxTz_dx = spatial_derivative(WxTz_sim, WxTz, WxTz_sip) / dx;
	precision dWxTz_dy = spatial_derivative(WxTz_sjm, WxTz, WxTz_sjp) / dy;
	precision dWxTz_dn = spatial_derivative(WxTz_skm, WxTz, WxTz_skp) / dn;

	n += 2;

	precision WyTz_sim = qi1[n];		// WyTz [i-1, i+1]
	precision WyTz_sip = qi1[n + 1];

	precision WyTz_sjm = qj1[n];		// WyTz [j-1, j+1]
	precision WyTz_sjp = qj1[n + 1];

	precision WyTz_skm = qk1[n];		// WyTz [k-1, k+1]
	precision WyTz_skp = qk1[n + 1];

	precision dWyTz_dx = spatial_derivative(WyTz_sim, WyTz, WyTz_sip) / dx;
	precision dWyTz_dy = spatial_derivative(WyTz_sjm, WyTz, WyTz_sjp) / dy;
	precision dWyTz_dn = spatial_derivative(WyTz_skm, WyTz, WyTz_skp) / dn;

	n += 2;

	precision WnTz_sim = qi1[n];		// WnTz [i-1, i+1]
	precision WnTz_sip = qi1[n + 1];

	precision WnTz_sjm = qj1[n];		// WnTz [j-1, j+1]
	precision WnTz_sjp = qj1[n + 1];

	precision WnTz_skm = qk1[n];		// WnTz [k-1, k+1]
	precision WnTz_skp = qk1[n + 1];

	precision dWnTz_dx = spatial_derivative(WnTz_sim, WnTz, WnTz_sip) / dx;
	precision dWnTz_dy = spatial_derivative(WnTz_sjm, WnTz, WnTz_sjp) / dy;
	precision dWnTz_dn = spatial_derivative(WnTz_skm, WnTz, WnTz_skp) / dn;
#else
	precision dWtTz_dx = 0.0;
	precision dWtTz_dy = 0.0;
	precision dWtTz_dn = 0.0;

	precision dWxTz_dx = 0.0;
	precision dWxTz_dy = 0.0;
	precision dWxTz_dn = 0.0;

	precision dWyTz_dx = 0.0;
	precision dWyTz_dy = 0.0;
	precision dWyTz_dn = 0.0;

	precision dWnTz_dx = 0.0;
	precision dWnTz_dy = 0.0;
	precision dWnTz_dn = 0.0;
#endif


	// fluid velocity
	precision ut_sim = ui1[0];		// ut [i-1, i+1]
	precision ut_sip = ui1[1];
	precision ux_sim = ui1[2];		// ux [i-1, i+1]
	precision ux_sip = ui1[3];
	precision uy_sim = ui1[4];		// uy [i-1, i+1]
	precision uy_sip = ui1[5];
	precision un_sim = ui1[6];		// un [i-1, i+1]
	precision un_sip = ui1[7];

	precision ut_sjm = uj1[0];		// ut [j-1, j+1]
	precision ut_sjp = uj1[1];
	precision ux_sjm = uj1[2];		// ux [j-1, j+1]
	precision ux_sjp = uj1[3];
	precision uy_sjm = uj1[4];		// uy [j-1, j+1]
	precision uy_sjp = uj1[5];
	precision un_sjm = uj1[6];		// un [j-1, j+1]
	precision un_sjp = uj1[7];

	precision ut_skm = uk1[0];		// ut [k-1, k+1]
	precision ut_skp = uk1[1];
	precision ux_skm = uk1[2];		// ux [k-1, k+1]
	precision ux_skp = uk1[3];
	precision uy_skm = uk1[4];		// uy [k-1, k+1]
	precision uy_skp = uk1[5];
	precision un_skm = uk1[6];		// un [k-1, k+1]
	precision un_skp = uk1[7];

	// spatial derivatives
	precision dut_dx = spatial_derivative(ut_sim, ut, ut_sip) / dx;
	precision dut_dy = spatial_derivative(ut_sjm, ut, ut_sjp) / dy;
	precision dut_dn = spatial_derivative(ut_skm, ut, ut_skp) / dn;

	precision dux_dx = spatial_derivative(ux_sim, ux, ux_sip) / dx;
	precision dux_dy = spatial_derivative(ux_sjm, ux, ux_sjp) / dy;
	precision dux_dn = spatial_derivative(ux_skm, ux, ux_skp) / dn;

	precision duy_dx = spatial_derivative(uy_sim, uy, uy_sip) / dx;
	precision duy_dy = spatial_derivative(uy_sjm, uy, uy_sjp) / dy;
	precision duy_dn = spatial_derivative(uy_skm, uy, uy_skp) / dn;

	precision dun_dx = spatial_derivative(un_sim, un, un_sip) / dx;
	precision dun_dy = spatial_derivative(un_sjm, un, un_sjp) / dy;
	precision dun_dn = spatial_derivative(un_skm, un, un_skp) / dn;

	// time derivatives
	precision dut_dt = (ut - ut_p) / dt;
	precision dux_dt = (ux - ux_p) / dt;
	precision duy_dt = (uy - uy_p) / dt;
	precision dun_dt = (un - un_p) / dt;

	// radial velocity derivatives
	precision duT2_dt = ux * dux_dt  +  uy * duy_dt;
	precision duT2_dx = ux * dux_dx  +  uy * duy_dx;
	precision duT2_dy = ux * dux_dy  +  uy * duy_dy;
	precision duT2_dn = ux * dux_dn  +  uy * duy_dn;


	// scalar expansion rate
	precision theta = dut_dt  +  dux_dx  +  duy_dy  +  dun_dn  +  ut / t;


	// spatial velocity
	precision vx = ux / ut;
	precision vy = uy / ut;
	precision vn = un / ut;

	precision vx_sim = ux_sim / ut_sim;
	precision vx_sip = ux_sip / ut_sip;

	precision vy_sjm = uy_sjm / ut_sjm;
	precision vy_sjp = uy_sjp / ut_sjp;

	precision vn_skm = un_skm / ut_skm;
	precision vn_skp = un_skp / ut_skp;

	// spatial derivatives
	precision dvx_dx = spatial_derivative(vx_sim, vx, vx_sip) / dx;
	precision dvy_dy = spatial_derivative(vy_sjm, vy, vy_sjp) / dy;
	precision dvn_dn = spatial_derivative(vn_skm, vn, vn_skp) / dn;

	// divergence of v
	precision div_v = dvx_dx + dvy_dy + dvn_dn;
	//precision div_v = (dux_dx  +  duy_dy  +  dun_dn  -  vx * dut_dx  -  vy * dut_dy  -  vn * dut_dn) / ut;

	// other spatial velocity derivatives (move in if/else later)
	precision dvx_dn = (dux_dn  -  vx * dut_dn) / ut;
	precision dvy_dn = (duy_dn  -  vy * dut_dn) / ut;
	//:::::::::::::::::::::::::::::::::::::::::::::


	// longitudinal vector: spatial derivatives
	precision dzt_dt = t * (dun_dt  -  un * duT2_dt / utperp2) / utperp  +  zt / t;
	precision dzt_dx = t * (dun_dx  -  un * duT2_dx / utperp2) / utperp;
	precision dzt_dy = t * (dun_dy  -  un * duT2_dy / utperp2) / utperp;
	precision dzt_dn = t * (dun_dn  -  un * duT2_dn / utperp2) / utperp;

	precision dzn_dt = (dut_dt  -  ut * duT2_dt / utperp2) / (t * utperp)  -  zn / t;
	precision dzn_dx = (dut_dx  -  ut * duT2_dx / utperp2) / (t * utperp);
	precision dzn_dy = (dut_dy  -  ut * duT2_dy / utperp2) / (t * utperp);
	precision dzn_dn = (dut_dn  -  ut * duT2_dn / utperp2) / (t * utperp);


	// longitudinal and transverse expansion rates: thetaL = z_\mu Dz u^\mu, thetaT = NablaT_\mu u^\mu
	precision thetaL = - zt2 * dut_dt  +  t2 * zn2 * dun_dn  +  ztzn * (t2 * dun_dt  -  dut_dn)  +  t * zn2 * ut;
	precision thetaT = theta - thetaL;


#if (NUMBER_OF_RESIDUAL_CURRENTS != 0)
	transverse_projection Xi(ut, ux, uy, un, zt, zn, t2);	// Xi^{\mu\nu}
	double_transverse_projection Xi_2(Xi, t2, t4);			// Xi^{\mu\nu\alpha\beta}

	
	// acceleration = D u^\mu
	precision at = ut * dut_dt  +  ux * dut_dx  +  uy * dut_dy  +  un * dut_dn  +  t * un2;
	precision ax = ut * dux_dt  +  ux * dux_dx  +  uy * dux_dy  +  un * dux_dn;
	precision ay = ut * duy_dt  +  ux * duy_dx  +  uy * duy_dy  +  un * duy_dn;
	precision an = ut * dun_dt  +  ux * dun_dx  +  uy * dun_dy  +  un * dun_dn  +  2.0 * utun / t;


	// covariant time derivatives of z = D z^\mu
	precision D_zt = ut * dzt_dt  +  ux * dzt_dx  +  uy * dzt_dy  +  un * dzt_dn  +  t * unzn;
	precision D_zn = ut * dzn_dt  +  ux * dzn_dx  +  uy * dzn_dy  +  un * dzn_dn  +  (ut * zn  +  un * zt) / t;


	// longitudinal covariant derivative of u = Dz u^\mu
	precision Dz_ut = - zt * dut_dt  -  zn * dut_dn  -  t * unzn;
	precision Dz_ux = - zt * dux_dt  -  zn * dux_dn;
	precision Dz_uy = - zt * duy_dt  -  zn * duy_dn;
	precision Dz_un = - zt * dun_dt  -  zn * dun_dn  -  (ut * zn  +  un * zt) / t;


	// covariant transverse derivative of u contracted with z = z_\nu NablaT^\mu u^\nu (first compute z_\nu D^\mu u^\nu )
	precision z_NabTt_u =   zt * dut_dt  -  t2zn * dun_dt  -  t * unzn;
	precision z_NabTx_u = -(zt * dut_dx  -  t2zn * dun_dx);
	precision z_NabTy_u = -(zt * dut_dy  -  t2zn * dun_dy);
	precision z_NabTn_u = -(zt * dut_dn  -  t2zn * dun_dn  +  t * (un * zt  -  ut * zn)) / t2;

	Xi.transverse_project_vector(z_NabTt_u, z_NabTx_u, z_NabTy_u, z_NabTn_u);


	// transverse shear velocity tensor = sigmaT^{\mu\nu} (first compute D^{(\mu) u^{\nu)} )
	precision sTtt = dut_dt;
	precision sTtx = 0.5 * (dux_dt  -  dut_dx);
	precision sTty = 0.5 * (duy_dt  -  dut_dy);
	precision sTtn = 0.5 * (dun_dt  +  un / t  -  (dut_dn  +  t * un) / t2);
	precision sTxx = - dux_dx;
	precision sTxy = - 0.5 * (dux_dy  +  duy_dx);
	precision sTxn = - 0.5 * (dun_dx  +  dux_dn / t2);
	precision sTyy = - duy_dy;
	precision sTyn = - 0.5 * (dun_dy  +  duy_dn / t2);
	precision sTnn = - (dun_dn  +  ut / t) / t2;

	Xi_2.double_transverse_project_tensor(sTtt, sTtx, sTty, sTtn, sTxx, sTxy, sTxn, sTyy, sTyn, sTnn);
#endif


	// relaxation times
#ifdef CONFORMAL_EOS
	precision T = powf(e / EOS_FACTOR, 0.25);

	precision taupiInv = 0.2 * T / etabar;
	precision taubulkInv = 0.0;
#else
	precision taupiInv = 1.0;	// fill in quasiparticle later
	precision taubulkInv = 1.0;
#endif


	// anisotropic transport coefficients
	transport_coefficients aniso;
	aniso.compute_transport_coefficients(e, pl, pt);

	// pl coefficients
	precision zeta_LL = aniso.zeta_LL;
	precision zeta_TL = aniso.zeta_TL;
#ifdef WTZMU
	precision lambda_WuL = aniso.lambda_WuL;
	precision lambda_WTL = aniso.lambda_WTL;
#endif
#ifdef PIMUNU
	precision lambda_piL = aniso.lambda_piL;
#endif
	
	// pt coefficients
#if (PT_MATCHING == 1)
	precision zeta_LT = aniso.zeta_LT;
	precision zeta_TT = aniso.zeta_TT;
#ifdef WTZMU
	precision lambda_WuT = aniso.lambda_WuT;
	precision lambda_WTT = aniso.lambda_WTT;
#endif
#ifdef PIMUNU
	precision lambda_piTW = aniso.lambda_piTW;
#endif
#endif

	// WTz coefficients
#ifdef WTZMU
	precision eta_uW = aniso.eta_uW;
	precision eta_TW = aniso.eta_TW;
	precision tau_zW = aniso.tau_zW;
	precision delta_WW = aniso.delta_WW;
	precision lambda_WuW = aniso.lambda_WuW;
	precision lambda_WTW = aniso.lambda_WTW;
#ifdef PIMUNU
	precision lambda_piuW = aniso.lambda_piuW;
	precision lambda_piTW = aniso.lambda_piTW;
#endif
#endif

	// piT coefficients
#ifdef PIMUNU
	precision eta_T = aniso.eta_T;
	precision delta_pipi = aniso.delta_pipi;
	precision tau_pipi = aniso.tau_pipi;
	precision lambda_pipi = aniso.lambda_pipi;
#ifdef WTZMU
	precision lambda_Wupi = aniso.lambda_Wupi;
	precision lambda_WTpi = aniso.lambda_WTpi;
#endif
#endif


	// L^munu components and derivatives
	precision dp  = pl - pt;

	precision Ltt = dp * zt2;
	precision Ltn = dp * ztzn;

	// L^munu components and derivatives
	precision dLtt_dx = (dpl_dx - dpt_dx) * zt2  +  2.0 * dp * zt * dzt_dx;
	precision dLtt_dy = (dpl_dy - dpt_dy) * zt2  +  2.0 * dp * zt * dzt_dy;
	precision dLtt_dn = (dpl_dn - dpt_dn) * zt2  +  2.0 * dp * zt * dzt_dn;

	precision dLtn_dx = (dpl_dx - dpt_dx) * ztzn  +  dp * (dzt_dx * zn  +  dzn_dx * zt);
	precision dLtn_dy = (dpl_dy - dpt_dy) * ztzn  +  dp * (dzt_dy * zn  +  dzn_dy * zt);
	precision dLtn_dn = (dpl_dn - dpt_dn) * ztzn  +  dp * (dzt_dn * zn  +  dzn_dn * zt);

	precision dLnn_dn = (dpl_dn - dpt_dn) * zn2  +  2.0 * dp * zn * dzn_dn;


#ifdef WTZMU
	// W^munu components
	precision Wtt = 2.0 * WtTz * zt;
	precision Wtx = WxTz * zt;
	precision Wty = WyTz * zt;
	precision Wtn = WtTz * zn  +  WnTz * zt;
	precision Wxn = WxTz * zn;
	precision Wyn = WyTz * zn;
	precision Wnn = 2.0 * WnTz * zn;

	// W^munu derivatives in conservation laws
	precision dWtt_dx = 2.0 * (dWtTz_dx * zt  +  WtTz * dzt_dx);
	precision dWtt_dy = 2.0 * (dWtTz_dy * zt  +  WtTz * dzt_dy);
	precision dWtt_dn = 2.0 * (dWtTz_dn * zt  +  WtTz * dzt_dn);

	precision dWtx_dx = dWxTz_dx * zt  +  WxTz * dzt_dx;
	precision dWtx_dy = dWxTz_dy * zt  +  WxTz * dzt_dy;
	precision dWtx_dn = dWxTz_dn * zt  +  WxTz * dzt_dn;

	precision dWty_dx = dWyTz_dx * zt  +  WyTz * dzt_dx;
	precision dWty_dy = dWyTz_dy * zt  +  WyTz * dzt_dy;
	precision dWty_dn = dWyTz_dn * zt  +  WyTz * dzt_dn;

	precision dWtn_dx = dWtTz_dx * zn  +  WtTz * dzn_dx  +  dWnTz_dx * zt  +  WnTz * dzt_dx;
	precision dWtn_dy = dWtTz_dy * zn  +  WtTz * dzn_dy  +  dWnTz_dy * zt  +  WnTz * dzt_dy;
	precision dWtn_dn = dWtTz_dn * zn  +  WtTz * dzn_dn  +  dWnTz_dn * zt  +  WnTz * dzt_dn;

	precision dWxn_dx = dWxTz_dx * zn  +  WxTz * dzn_dx;
	precision dWxn_dn = dWxTz_dn * zn  +  WxTz * dzn_dn;

	precision dWyn_dy = dWyTz_dy * zn  +  WyTz * dzn_dy;
	precision dWyn_dn = dWyTz_dn * zn  +  WyTz * dzn_dn;

	precision dWnn_dn = 2.0 * (dWnTz_dn * zn  +  WnTz * dzn_dn);


	// 2nd order gradient terms in dpl and dpt
	precision WTz_D_z = WtTz * D_zt  -  t2 * WnTz * D_zn;
	precision WTz_Dz_u = WtTz * Dz_ut  -  WxTz * Dz_ux  -  WyTz * Dz_uy  -  t2 * WnTz * Dz_un;
	precision WTz_z_NabT_u = WtTz * z_NabTt_u  -  WxTz * z_NabTx_u  -  WyTz * z_NabTy_u  -  t2 * WnTz * z_NabTn_u;

	precision IplW = - 2.0 * WTz_D_z  +  lambda_WuL * WTz_Dz_u  +  lambda_WTL * WTz_z_NabT_u;

#if (PT_MATCHING == 1)
	precision IptW = WTz_D_z  +  lambda_WuT * WTz_Dz_u  -  lambda_WTT * WTz_z_NabT_u;
#endif

	// gradient terms in dWTz (some are unprojected)
	precision It_W1 = 2.0 * eta_uW * Dz_ut;
	precision Ix_W1 = 2.0 * eta_uW * Dz_ux;
	precision Iy_W1 = 2.0 * eta_uW * Dz_uy;
	precision In_W1 = 2.0 * eta_uW * Dz_un;

	precision It_W2 = 2.0 * eta_TW * z_NabTt_u;
	precision Ix_W2 = 2.0 * eta_TW * z_NabTx_u;
	precision Iy_W2 = 2.0 * eta_TW * z_NabTy_u;
	precision In_W2 = 2.0 * eta_TW * z_NabTn_u;

	precision It_W3 = tau_zW * D_zt;
	precision In_W3 = tau_zW * D_zn;

	precision It_W4 = pitt * D_zt  -  t2 * pitn * D_zn;
	precision In_W4 = pitn * D_zt  -  t2 * pinn * D_zn;

	precision It_W5 = delta_WW * WtTz * thetaT;
	precision Ix_W5 = delta_WW * WxTz * thetaT;
	precision Iy_W5 = delta_WW * WyTz * thetaT;
	precision In_W5 = delta_WW * WnTz * thetaT;

	precision It_W6 = lambda_WuW * WtTz * thetaL;
	precision Ix_W6 = lambda_WuW * WxTz * thetaL;
	precision Iy_W6 = lambda_WuW * WyTz * thetaL;
	precision In_W6 = lambda_WuW * WnTz * thetaL;

	precision It_W7 = lambda_WTW * (sTtt * WtTz  -  sTtx * WxTz  -  sTty * WyTz  - t2 * sTtn * WnTz);
	precision Ix_W7 = lambda_WTW * (sTtx * WtTz  -  sTxx * WxTz  -  sTxy * WyTz  - t2 * sTxn * WnTz);
	precision Iy_W7 = lambda_WTW * (sTty * WtTz  -  sTxy * WxTz  -  sTyy * WyTz  - t2 * sTyn * WnTz);
	precision In_W7 = lambda_WTW * (sTtn * WtTz  -  sTxn * WxTz  -  sTyn * WyTz  - t2 * sTnn * WnTz);

	precision It_W8 = 0.0;		// keep vorticity term zero for now..
	precision Ix_W8 = 0.0;
	precision Iy_W8 = 0.0;
	precision In_W8 = 0.0;

#ifdef PIMUNU
	precision It_W9 = lambda_piuW * (pitt * Dz_ut  -  pitx * Dz_ux  -  pity * Dz_uy  -  t2 * pitn * Dz_un);
	precision Ix_W9 = lambda_piuW * (pitx * Dz_ut  -  pixx * Dz_ux  -  pixy * Dz_uy  -  t2 * pixn * Dz_un);
	precision Iy_W9 = lambda_piuW * (pity * Dz_ut  -  pixy * Dz_ux  -  piyy * Dz_uy  -  t2 * piyn * Dz_un);
	precision In_W9 = lambda_piuW * (pitn * Dz_ut  -  pixn * Dz_ux  -  piyn * Dz_uy  -  t2 * pinn * Dz_un);

	precision It_W10 = lambda_piTW * (pitt * z_NabTt_u  -  pitx * z_NabTx_u  -  pity * z_NabTy_u  -  t2 * pitn * z_NabTn_u);
	precision It_W10 = lambda_piTW * (pitx * z_NabTt_u  -  pixx * z_NabTx_u  -  pixy * z_NabTy_u  -  t2 * pixn * z_NabTn_u);
	precision It_W10 = lambda_piTW * (pity * z_NabTt_u  -  pixy * z_NabTx_u  -  piyy * z_NabTy_u  -  t2 * piyn * z_NabTn_u);
	precision It_W10 = lambda_piTW * (pitn * z_NabTt_u  -  pixn * z_NabTx_u  -  piyn * z_NabTy_u  -  t2 * pinn * z_NabTn_u);
#else
	precision It_W9 = 0.0;
	precision Ix_W9 = 0.0;
	precision Iy_W9 = 0.0;
	precision In_W9 = 0.0;

	precision It_W10 = 0.0;
	precision Ix_W10 = 0.0;
	precision Iy_W10 = 0.0;
	precision In_W10 = 0.0;
#endif

	// add source terms together
	precision It_W = It_W1  -  It_W2  -  It_W3  -  It_W4  +  It_W5  -  It_W6  +  It_W7  +  It_W8  +  It_W9  -  It_W10;
	precision Ix_W = Ix_W1  -  Ix_W2  			   		  +  Ix_W5  -  Ix_W6  +  Ix_W7  +  Ix_W8  +  Ix_W9  -  Ix_W10;
	precision Iy_W = Iy_W1  -  Iy_W2  			 		  +  Iy_W5  -  Iy_W6  +  Iy_W7  +  Iy_W8  +  Iy_W9  -  Iy_W10;
	precision In_W = In_W1  -  In_W2  -  In_W3  -  In_W4  +  In_W5  -  In_W6  +  In_W7  +  In_W8  +  In_W9  -  In_W10;

	Xi.transverse_project_vector(It_W, Ix_W, Iy_W, In_W);


	// Christofel terms (G_W^\mu  =  u^\alpha . \Gamma^\mu_{\alpha\beta} . WTz^\mu)
	precision Gt_W = t * un * WnTz;
	precision Gn_W = (ut * WnTz  +  un * WtTz) / t;


	// product rule terms (P_W^\mu  =  - u^\mu . Wtz^\nu . a_\nu  +  z^\mu . Wtz^\nu . Dz_\nu)
	precision Wtz_a  = WtTz * at  -  WxTz * ax  -  WyTz * ay  -  t2 * WnTz * an;
	precision Wtz_Dz = WtTz * D_zt  -  t2 * WnTz * D_zn;

	precision Pt_W = - ut * Wtz_a  +  zt * Wtz_Dz;
	precision Px_W = - ux * Wtz_a;
	precision Py_W = - uy * Wtz_a;
	precision Pn_W = - un * Wtz_a  +  zn * Wtz_Dz;

#else
	precision Wtt = 0.0;
	precision Wtx = 0.0;
	precision Wty = 0.0;
	precision Wtn = 0.0;

	precision Wxn = 0.0;
	precision Wyn = 0.0;
	precision Wnn = 0.0;

	precision dWtt_dx = 0.0;
	precision dWtt_dy = 0.0;
	precision dWtt_dn = 0.0;

	precision dWtx_dx = 0.0;
	precision dWtx_dy = 0.0;
	precision dWtx_dn = 0.0;

	precision dWty_dx = 0.0;
	precision dWty_dy = 0.0;
	precision dWty_dn = 0.0;

	precision dWtn_dx = 0.0;
	precision dWtn_dy = 0.0;
	precision dWtn_dn = 0.0;

	precision dWxn_dx = 0.0;
	precision dWxn_dn = 0.0;

	precision dWyn_dy = 0.0;
	precision dWyn_dn = 0.0;

	precision dWnn_dn = 0.0;

	precision IplW = 0.0;

#if (PT_MATCHING == 1)
	precision IptW = 0.0
#endif

#endif


#ifdef PIMUNU
	precision pisT = pitt * sTtt  +  pixx * sTxx  +  piyy * sTyy  +  t2 * t2 * pinn * sTnn  +  2.0 * (pixy * sTxy  -  pitx * sTtx  -  pity * sTty  +  t2 * (pixn * sTxn  +  piyn * sTyn  -  pitn * sTtn));
#else
	precision pisT = 0.0;

#endif






	// conservation laws
	precision tnn = (e + pt) * un2  +  pt / t2  +  dp * zn2  +  Wnn  +  pinn;


	S[0] =	- (ttt / t  +  t * tnn)  +  div_v * (Ltt  +  Wtt  +  pitt  -  pt)
			+  vx * (dLtt_dx  +  dWtt_dx  +  dpitt_dx  -  dpt_dx)  -  dWtx_dx  -  dpitx_dx
			+  vy * (dLtt_dy  +  dWtt_dy  +  dpitt_dy  -  dpt_dy)  -  dWty_dy  -  dpity_dy
			+  vn * (dLtt_dn  +  dWtt_dn  +  dpitt_dn  -  dpt_dn)  -  dWtn_dn  -  dpitn_dn  -  dLtn_dn;


	S[1] =	- ttx / t  -  dpt_dx  +  div_v * (Wtx  +  pitx)
			+  vx * (dWtx_dx  +  dpitx_dx)  -  dpixx_dx
			+  vy * (dWtx_dy  +  dpitx_dy)  -  dpixy_dy
			+  vn * (dWtx_dn  +  dpitx_dn)  -  dpixn_dn  -  dWxn_dn
			-  0.5 * (vx * dLtn_dn  -  Ltn * dvx_dn);	// go over this line again


	S[2] =	- tty / t  -  dpt_dy  +  div_v * (Wty  +  pity)
			+  vx * (dWty_dx  +  dpity_dx)  -  dpixy_dx
			+  vy * (dWty_dy  +  dpity_dy)  -  dpiyy_dy
			+  vn * (dWty_dn  +  dpity_dn)  -  dpiyn_dn  -  dWyn_dn
			-  0.5 * (vy * dLtn_dn  -  Ltn * dvy_dn);  // go over this line again


	S[3] =	- 3.0 * ttn / t  -  dpt_dn / t2  +  div_v * (Ltn  +  Wtn  +  pitn)
			+  vx * (dLtn_dx  +  dWtn_dx  +  dpitn_dx)  -  dWxn_dx  -  dpixn_dx
			+  vy * (dLtn_dy  +  dWtn_dy  +  dpitn_dy)  -  dWyn_dy  -  dpiyn_dy
			+  vn * (dLtn_dn  +  dWtn_dn  +  dpitn_dn)  -  dWnn_dn  -  dpinn_dn  -  dLnn_dn;




	// pl relaxation equation
	precision dpl = - dp * taupiInv / 1.5  +  zeta_LL * thetaL  +  zeta_TL * thetaT  +  IplW;
	S[4] =	dpl / ut  +  div_v * pl;


	a = 5;		// reset index

#if (PT_MATCHING == 1)
	precision dpt =	dp * taupiInv / 3.0  +  zeta_LT * thetaL  +  zeta_TT * thetaT  +  IptW;
	S[a] = dpt / ut  +  div_v * pt;		a++;
#endif

#ifdef PIMUNU
	precision dpitt = - pitt * taupiInv  +  Itt_pi_pro  +  Ptt_pi  -  Gtt_pi;	// imagine it look like this 
	precision dpitx = - pitx * taupiInv  +  Itx_pi_pro  +  Ptx_pi;
	precision dpity = - pity * taupiInv  +  Ity_pi_pro  +  Pty_pi;
	precision dpitn = - pitn * taupiInv  +  Itn_pi_pro  +  Ptn_pi  -  Gtn_pi;
	precision dpixx = - pixx * taupiInv  +  Ixx_pi_pro  +  Pxx_pi;
	precision dpixy = - pixy * taupiInv  +  Ixy_pi_pro  +  Pxy_pi;
	precision dpixn = - pixn * taupiInv  +  Ixn_pi_pro  +  Pxn_pi;
	precision dpiyy = - piyy * taupiInv  +  Iyy_pi_pro  +  Pyy_pi;
	precision dpiyn = - piyn * taupiInv  +  Iyn_pi_pro  +  Pyn_pi;
	precision dpinn = - pinn * taupiInv  +  Inn_pi_pro  +  Pnn_pi  -  Gnn_pi; 

	S[a] = dpitt / ut  +  div_v * pitt;		a++;
	S[a] = dpitx / ut  +  div_v * pitx;		a++;
	S[a] = dpity / ut  +  div_v * pity;		a++;
	S[a] = dpitn / ut  +  div_v * pitn;		a++;
	S[a] = dpixx / ut  +  div_v * pixx;		a++;
	S[a] = dpixy / ut  +  div_v * pixy;		a++;
	S[a] = dpixn / ut  +  div_v * pixn;		a++;
	S[a] = dpiyy / ut  +  div_v * piyy;		a++;
	S[a] = dpiyn / ut  +  div_v * piyn;		a++;
	S[a] = dpinn / ut  +  div_v * pinn;		a++;
#endif

#ifdef WTZMU
	precision dWtTz = - WtTz * taupiInv  +  It_W  +  Pt_W  -  Gt_W;	
	precision dWxTz = - WxTz * taupiInv  +  Ix_W  +  Px_W;
	precision dWyTz = - WyTz * taupiInv  +  Iy_W  +  Py_W;
	precision dWnTz = - WnTz * taupiInv  +  In_W  +  Pn_W  -  Gn_W;

	S[a] = dWtTz / ut  +  div_v * WtTz;		a++;
	S[a] = dWxTz / ut  +  div_v * WxTz;		a++;
	S[a] = dWyTz / ut  +  div_v * WyTz;		a++;
	S[a] = dWnTz / ut  +  div_v * WnTz;	
#endif

}




