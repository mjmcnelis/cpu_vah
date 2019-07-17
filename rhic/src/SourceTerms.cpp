/*
 * SourceTerms.cpp
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */
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
	precision t2 = t * t;		// useful expressions
	precision t3 = t2 * t;

	precision tut  = t * ut;
	precision ut2  = ut * ut;
	precision un2  = un * un;
	precision utun = ut * un;


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

	// also need to check for missing conservation law terms (was in x, y)

	


	// fluid velocity 
	//:::::::::::::::::::::::::::::::::::::::::::::
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

	precision utperp  = sqrt(1.0  +  ux * ux  +  uy * uy);
	precision utperp2 = utperp * utperp;

	// scalar expansion rate
	precision theta = dut_dt  +  dux_dx  +  duy_dy  +  dun_dn  +  ut / t;
	//:::::::::::::::::::::::::::::::::::::::::::::



	// spatial velocity 
	//:::::::::::::::::::::::::::::::::::::::::::::
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



	// longitudinal basis vector
	//:::::::::::::::::::::::::::::::::::::::::::::
	precision zt = t * un / utperp;
	precision zn = ut / t / utperp;

	precision zt2  = zt * zt;
	precision ztzn = zt * zn;
	precision zn2  = zn * zn;

	precision t2zn = t2 * zn;

	precision unzn = un * zn;
	precision unzt_minus_utzn = un * zt  -  ut * zn;

	// spatial derivatives
	precision dzt_dt = t * (dun_dt  -  un * duT2_dt / utperp2) / utperp  +  zt / t;
	precision dzt_dx = t * (dun_dx  -  un * duT2_dx / utperp2) / utperp;
	precision dzt_dy = t * (dun_dy  -  un * duT2_dy / utperp2) / utperp;
	precision dzt_dn = t * (dun_dn  -  un * duT2_dn / utperp2) / utperp;

	precision dzn_dt = (dut_dt  -  ut * duT2_dt / utperp2) / (t * utperp)  -  zn / t;
	precision dzn_dx = (dut_dx  -  ut * duT2_dx / utperp2) / (t * utperp);
	precision dzn_dy = (dut_dy  -  ut * duT2_dy / utperp2) / (t * utperp);
	precision dzn_dn = (dut_dn  -  ut * duT2_dn / utperp2) / (t * utperp);



	// longitudinal and transverse expansion rates: thetaL = z_\mu Dz u^\mu, thetaT = NablaT_\mu u^\mu
	//:::::::::::::::::::::::::::::::::::::::::::::
	precision thetaL = - zt2 * dut_dt  +  t2 * zn2 * dun_dn  +  ztzn * (t2 * dun_dt  -  dut_dn)  +  t * zn2 * ut;
	precision thetaT = theta - thetaL;
	precision thetaT_over2 = thetaT / 2.0;
	//:::::::::::::::::::::::::::::::::::::::::::::

	

#if (NUMBER_OF_RESIDUAL_CURRENTS != 0)
	// transverse projection tensor
	//:::::::::::::::::::::::::::::::::::::::::::::
	precision Xitt = 1.0  -  ut2  +  zt2;
	precision Xitx = - ut * ux;
	precision Xity = - ut * uy;
	precision Xitn = - utun  +  ztzn;
	precision Xixx = - 1.0  -  ux * ux;
	precision Xixy = - ux * uy;
	precision Xixn = - ux * un;
	precision Xiyy = - 1.0  -  uy * uy;
	precision Xiyn = - uy * un;
	precision Xinn = - 1.0 / t2  -  un2  + zn2;  
	//:::::::::::::::::::::::::::::::::::::::::::::



	// covariant time derivatives
	precision D_zt = ut * dzt_dt  +  ux * dzt_dx  +  uy * dzt_dy  +  un * dzt_dn  +  t * un * zn;
	precision D_zn = ut * dzn_dt  +  ux * dzn_dx  +  uy * dzn_dy  +  un * dzn_dn  +  (ut * zn  +  un * zt) / t;
	//:::::::::::::::::::::::::::::::::::::::::::::



	// longitudinal covariant derivative of u: Dz u^\mu
	//:::::::::::::::::::::::::::::::::::::::::::::
	precision Dz_ut = - zt * dut_dt  -  zn * dut_dn  -  t * zn * un; 
	precision Dz_ux = - zt * dux_dt  -  zn * dux_dn; 
	precision Dz_uy = - zt * duy_dt  -  zn * duy_dn; 
	precision Dz_un = - zt * dun_dt  -  zn * dun_dn  -  (ut * zn  +  un * zt) / t; 
	//:::::::::::::::::::::::::::::::::::::::::::::



	// transverse derivatives of u contracted in longitudinal direction: z_\nu NablaT^\mu u^\nu
	//:::::::::::::::::::::::::::::::::::::::::::::
	precision znu_dt_unu = zt * dut_dt  -  t2zn * dun_dt;
	precision znu_dx_unu = zt * dut_dx  -  t2zn * dun_dx;
	precision znu_dy_unu = zt * dut_dy  -  t2zn * dun_dy;
	precision znu_dn_unu = zt * dut_dn  -  t2zn * dun_dn;

	precision znu_NablaTt_unu = Xitt * (znu_dt_unu  - t * unzn)  +  Xitx * znu_dx_unu  +  Xity * znu_dy_unu  +  Xitn * (znu_dn_unu  +  t * unzt_minus_utzn);
	precision znu_NablaTx_unu = Xitx * (znu_dt_unu  - t * unzn)  +  Xixx * znu_dx_unu  +  Xixy * znu_dy_unu  +  Xixn * (znu_dn_unu  +  t * unzt_minus_utzn);
	precision znu_NablaTy_unu = Xity * (znu_dt_unu  - t * unzn)  +  Xixy * znu_dx_unu  +  Xiyy * znu_dy_unu  +  Xiyn * (znu_dn_unu  +  t * unzt_minus_utzn);
	precision znu_NablaTn_unu = Xitn * (znu_dt_unu  - t * unzn)  +  Xixn * znu_dx_unu  +  Xiyn * znu_dy_unu  +  Xinn * (znu_dn_unu  +  t * unzt_minus_utzn);
	//:::::::::::::::::::::::::::::::::::::::::::::



	// transverse shear stress tensor: \sigmaT^{\mu\nu} (explicit expressions very cumbersome and bug prone...is there a better way?)
	//:::::::::::::::::::::::::::::::::::::::::::::
	precision dut_dx_minus_dux_dt 	= dut_dx  -  dux_dt;
	precision dut_dy_minus_duy_dt 	= dut_dy  -  duy_dt;
	precision dut_dn_minus_t2dun_dt = dut_dn  -  t2 * dun_dt;
	precision dux_dy_plus_duy_dx 	= dux_dy  +  duy_dx;
	precision dux_dn_plus_t2dun_dx 	= dux_dn  +  t2 * dun_dx;
	precision duy_dn_plus_t2dun_dy 	= duy_dn  +  t2 * dun_dy;


	precision sTtt =	Xitt * (Xitt * dut_dt  +  Xitx * dut_dx_minus_dux_dt  +  Xity * dut_dy_minus_duy_dt  +  Xitn * dut_dn_minus_t2dun_dt - thetaT_over2 + Xinn*tut/2.)
					-  	Xitx * (Xitx * dux_dx  +  Xity * dux_dy_plus_duy_dx	  +  Xitn * dux_dn_plus_t2dun_dx)  
					-  	Xity * (Xity * duy_dy  +  Xitn * duy_dn_plus_t2dun_dy)
					-  	Xitn * (Xitn * (t2 * dun_dn  +  tut));
	

	precision sTtx = 	Xitx * (Xitt * dut_dt  -  Xixx * dux_dx  -  thetaT_over2  +  Xinn * tut / 2.)  -  Xity * Xixy * duy_dy  -  Xitn * Xixn * (t2 * dun_dn  +  tut)
					+	0.5 * (		(Xitx * Xitx  +  Xixx * Xitt) * dut_dx_minus_dux_dt
								+	(Xity * Xitx  +  Xixy * Xitt) * dut_dy_minus_duy_dt
								+	(Xitn * Xitx  +  Xixn * Xitt) * dut_dn_minus_t2dun_dt	
								-	(Xity * Xixx  +  Xixy * Xitx) * dux_dy_plus_duy_dx
								-	(Xitn * Xixx  +  Xixn * Xitx) * dux_dn_plus_t2dun_dx
								-	(Xitn * Xixy  +  Xixn * Xity) * duy_dn_plus_t2dun_dy	);


	precision sTty = 	Xity * (Xitt * dut_dt  -  Xixy * dux_dx  -  thetaT_over2  +  Xinn * tut / 2.)  -  Xity * Xiyy * duy_dy  -  Xitn * Xiyn * (t2 * dun_dn  +  tut)
					+	0.5 * (		(Xitx * Xity  +  Xixy * Xitt) * dut_dx_minus_dux_dt
								+	(Xity * Xity  +  Xiyy * Xitt) * dut_dy_minus_duy_dt
								+	(Xitn * Xity  +  Xiyn * Xitt) * dut_dn_minus_t2dun_dt	
								-	(Xity * Xixy  +  Xiyy * Xitx) * dux_dy_plus_duy_dx
								-	(Xitn * Xixy  +  Xiyn * Xitx) * dux_dn_plus_t2dun_dx
								-	(Xitn * Xiyy  +  Xiyn * Xity) * duy_dn_plus_t2dun_dy	);


	precision sTtn = 	Xitn * (Xitt * dut_dt  -  Xixn * dux_dx  -  thetaT_over2  +  Xinn * tut / 2.)  -  Xity * Xiyn * duy_dy  -  Xitn * Xinn * (t2 * dun_dn  +  tut)
					+	0.5 * (		(Xitx * Xitn  +  Xixn * Xitt) * dut_dx_minus_dux_dt
								+	(Xity * Xitn  +  Xiyn * Xitt) * dut_dy_minus_duy_dt
								+	(Xitn * Xitn  +  Xinn * Xitt) * dut_dn_minus_t2dun_dt	
								-	(Xity * Xixn  +  Xiyn * Xitx) * dux_dy_plus_duy_dx
								-	(Xitn * Xixn  +  Xinn * Xitx) * dux_dn_plus_t2dun_dx
								-	(Xitn * Xiyn  +  Xinn * Xity) * duy_dn_plus_t2dun_dy	);

				
	precision sTxx =	Xitx * (Xitx * dut_dt  +  Xixx * dut_dx_minus_dux_dt  +  Xixy * dut_dy_minus_duy_dt   +  Xixn * dut_dn_minus_t2dun_dt)  
					-  	Xixx * (Xixx * dux_dx  +  Xixy * dux_dy_plus_duy_dx	  +  Xixn * dux_dn_plus_t2dun_dx  +  thetaT_over2  -  Xinn * tut / 2.)  
					-  	Xixy * (Xixy * duy_dy  +  Xixn * duy_dn_plus_t2dun_dy)
					-  	Xixn * (Xixn * (t2 * dun_dn  +  tut));



	precision sTxy =	- Xixy * (Xixx * dux_dx  +  Xiyy * duy_dy  +  thetaT_over2  -  Xinn * tut / 2.)  +  Xitx * Xity * dut_dt  -  Xixn * Xiyn * (t2 * dun_dn  +  tut)

					+	0.5 * (		(Xixx * Xity  +  Xixy * Xitx) * dut_dx_minus_dux_dt
								+	(Xixy * Xity  +  Xiyy * Xitx) * dut_dy_minus_duy_dt
								+	(Xixn * Xity  +  Xiyn * Xitx) * dut_dn_minus_t2dun_dt	
								-	(Xixy * Xixy  +  Xiyy * Xixx) * dux_dy_plus_duy_dx
								-	(Xixn * Xixy  +  Xiyn * Xixx) * dux_dn_plus_t2dun_dx
								-	(Xixn * Xiyy  +  Xiyn * Xixy) * duy_dn_plus_t2dun_dy	);


	precision sTxn =	- Xixn * (Xixx * dux_dx  +  Xiyn * duy_dy  +  thetaT_over2  -  Xinn * tut / 2.)  +  Xitx * Xitn * dut_dt  -  Xixn * Xinn * (t2 * dun_dn  +  tut)

					+	0.5 * (		(Xixx * Xitn  +  Xixn * Xitx) * dut_dx_minus_dux_dt
								+	(Xixy * Xitn  +  Xiyn * Xitx) * dut_dy_minus_duy_dt
								+	(Xixn * Xitn  +  Xinn * Xitx) * dut_dn_minus_t2dun_dt	
								-	(Xixy * Xixn  +  Xiyn * Xixx) * dux_dy_plus_duy_dx
								-	(Xixn * Xixn  +  Xinn * Xixx) * dux_dn_plus_t2dun_dx
								-	(Xixn * Xiyn  +  Xinn * Xixy) * duy_dn_plus_t2dun_dy	);


	precision sTyy =	Xity * (Xity * dut_dt  +  Xixy * dut_dx_minus_dux_dt  +  Xiyy * dut_dy_minus_duy_dt   +  Xiyn * dut_dn_minus_t2dun_dt)  
					-  	Xixy * (Xixy * dux_dx  +  Xiyy * dux_dy_plus_duy_dx	  +  Xiyn * dux_dn_plus_t2dun_dx)  
					-  	Xiyy * (Xiyy * duy_dy  +  Xiyn * duy_dn_plus_t2dun_dy +  thetaT_over2  -  Xinn * tut / 2.)
					-  	Xiyn * (Xiyn * (t2 * dun_dn  +  tut));


	precision sTyn =	- Xiyn * (Xixy * dux_dx  +  Xiyy * duy_dy  +  thetaT_over2  -  Xinn * tut / 2.)  +  Xity * Xitn * dut_dt  -  Xiyn * Xinn * (t2 * dun_dn  +  tut)

					+	0.5 * (		(Xixy * Xitn  +  Xixn * Xity) * dut_dx_minus_dux_dt
								+	(Xiyy * Xitn  +  Xiyn * Xity) * dut_dy_minus_duy_dt
								+	(Xiyn * Xitn  +  Xinn * Xity) * dut_dn_minus_t2dun_dt	

								-	(Xiyy * Xixn  +  Xiyn * Xixy) * dux_dy_plus_duy_dx
								-	(Xiyn * Xixn  +  Xinn * Xixy) * dux_dn_plus_t2dun_dx

								-	(Xiyn * Xiyn  +  Xinn * Xiyy) * duy_dn_plus_t2dun_dy	);


	precision sTnn =	Xitn * (Xitn * dut_dt  +  Xixn * dut_dx_minus_dux_dt  +  Xiyn * dut_dy_minus_duy_dt   +  Xinn * dut_dn_minus_t2dun_dt)  
					-  	Xixn * (Xixn * dux_dx  +  Xiyn * dux_dy_plus_duy_dx	  +  Xinn * dux_dn_plus_t2dun_dx)  
					-  	Xiyn * (Xiyn * duy_dy  +  Xinn * duy_dn_plus_t2dun_dy)
					-  	Xinn * (Xinn * (t2 * dun_dn  +  tut / 2.) +  thetaT_over2);
#else
	precision D_zt = 0.0;
	precision D_zn = 0.0;

	precision Dz_ut = 0.0; 
	precision Dz_ux = 0.0; 
	precision Dz_uy = 0.0; 
	precision Dz_un = 0.0; 

	precision znu_NablaTt_unu = 0.0;
	precision znu_NablaTx_unu = 0.0;
	precision znu_NablaTy_unu = 0.0;
	precision znu_NablaTn_unu = 0.0;		

	precision sTtt = 0.0;
	precision sTtx = 0.0;
	precision sTty = 0.0;
	precision sTtn = 0.0;
	precision sTxx = 0.0;
	precision sTxy = 0.0;
	precision sTxn = 0.0;
	precision sTyy = 0.0;
	precision sTyn = 0.0;
	precision sTnn = 0.0;
#endif
					


	// relaxation times
	precision p = equilibriumPressure(e);
	precision T = effectiveTemperature(e);

	precision taupiInv = 0.2 * T / etabar;

	// anisotropic transport coefficients
	transport_coefficients aniso;
	aniso.compute_transport_coefficients(e, pl, pt);

	// pl coefficients
	precision zeta_LL = aniso.I_240  -  3.0 * pl;
	precision zeta_TL = aniso.I_221  -  pl;

	precision lambda_WuL = 1.0;	
	precision lambda_WTL = 1.0;
	
	

	// pt coefficients
#if (PT_MATCHING == 1)
	precision zeta_LT = aniso.I_221  -  pt;
	precision zeta_TT = 2.0 * (aniso.I_202 - pt);
	precision lambda_WuT = 1.0;
	precision lambda_WTT = 1.0;
#endif



	// L^munu components
	//:::::::::::::::::::::::::::::::::::::::::::::
	precision dp  = pl - pt;

	precision Ltt = dp * zt2;
	precision Ltn = dp * ztzn;
	precision tnn = (e + pt) * un2  +  pt / t2  +  dp * zn2;
	//:::::::::::::::::::::::::::::::::::::::::::::



	// L^munu derivatives in conservation laws
	//:::::::::::::::::::::::::::::::::::::::::::::
	precision dLtt_dx = (dpl_dx - dpt_dx) * zt2  +  2.0 * dp * zt * dzt_dx;
	precision dLtt_dy = (dpl_dy - dpt_dy) * zt2  +  2.0 * dp * zt * dzt_dy;
	precision dLtt_dn = (dpl_dn - dpt_dn) * zt2  +  2.0 * dp * zt * dzt_dn;

	precision dLtn_dx = (dpl_dx - dpt_dx) * ztzn  +  dp * (dzt_dx * zn  +  dzn_dx * zt);
	precision dLtn_dy = (dpl_dy - dpt_dy) * ztzn  +  dp * (dzt_dy * zn  +  dzn_dy * zt);
	precision dLtn_dn = (dpl_dn - dpt_dn) * ztzn  +  dp * (dzt_dn * zn  +  dzn_dn * zt);

	precision dLnn_dn = (dpl_dn - dpt_dn) * zn2  +  2.0 * dp * zn * dzn_dn;
	//:::::::::::::::::::::::::::::::::::::::::::::



	 
#ifdef WTZMU
	// W^munu components
	//:::::::::::::::::::::::::::::::::::::::::::::
	precision Wtt = 2.0 * WtTz * zt;
	precision Wtx = WxTz * zt;
	precision Wty = WyTz * zt;
	precision Wtn = WtTz * zn  +  WnTz * zt;
	
	precision Wxn = WxTz * zn;
	precision Wyn = WyTz * zn;
	precision Wnn = 2.0 * WnTz * zn;
	//:::::::::::::::::::::::::::::::::::::::::::::



	// W^munu derivatives in conservation laws
	//:::::::::::::::::::::::::::::::::::::::::::::
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
	//:::::::::::::::::::::::::::::::::::::::::::::



	// 2nd order gradient terms in relaxation equations (pl, pt)
	//:::::::::::::::::::::::::::::::::::::::::::::
	precision WTzmu_D_zmu  = WtTz * D_zt  -  t2 * WnTz * D_zn;

	precision WTzmu_Dz_umu = WtTz * Dz_ut  -  WxTz * Dz_ux  -  WyTz * Dz_uy  -  t2 * WnTz * Dz_un;
	
	precision WTzmu_znu_NablaTmu_unu = WtTz * znu_NablaTt_unu  -  WxTz * znu_NablaTx_unu  -  WyTz * znu_NablaTy_unu  -  t2 * WnTz * znu_NablaTt_unu;
	//:::::::::::::::::::::::::::::::::::::::::::::
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


	precision WTzmu_D_zmu  = 0.0;
	precision WTzmu_Dz_umu = 0.0;
	precision WTzmu_znu_NablaTmu_unu = 0.0;
#endif


#ifdef PIMUNU
	// compute 2nd order terms
	precision pimunu_sTmunu =	pitt * sTtt  +  pixx * sTxx  +  piyy * sTyy  +  t2 * t2 * pinn * sTnn  
							 +  2.0 * (pixy * sTxy  -  pitx * sTtx  -  pity * sTty  +  t2 * (pixn * sTxn  +  piyn * sTyn  -  pitn * sTtn));
#else
	// set 2nd order terms to zero 
#endif


	
	//:::::::::::::::::::::::::::::::::::::::::::::



	// conservation laws


	//S[0] = - (ttt / t  +  t * tnn)  +  div_v * (Ltt - pt)  +  vx * (dLtt_dx - dpt_dx)  +  vy * (dLtt_dy - dpt_dy)  +  vn * (dLtt_dn - dpt_dn)  - dLtn_dn;
	//S[1] = - ttx / t  -  dpt_dx  -  0.5 * (vx * dLtn_dn  -  Ltn * dvx_dn);
	//S[2] = - tty / t  -  dpt_dy  -  0.5 * (vy * dLtn_dn  -  Ltn * dvy_dn);


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
	precision dpl = - dp * taupiInv / 1.5  +  zeta_LL * thetaL  +  zeta_TL * thetaT  -  2.0 * WTzmu_D_zmu;
	S[4] =	dpl / ut  +  div_v * pl;


	a = 5;		// reset index

#if (PT_MATCHING == 1)
	precision dpt =	dp * taupiInv / 3.0  +  zeta_LT * thetaL  +  zeta_TT * thetaT  +  WTzmu_D_zmu;
	S[a] = dpt / ut  +  div_v * pt;
	a++;
#endif

}




