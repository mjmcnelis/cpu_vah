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


inline precision sign(precision x)
{
	if(x < 0.0) return -1.0;
	else return 1.0;
}


inline precision minmod(precision x, precision y)
{
	return (sign(x) + sign(y)) * fmin(fabs(x), fabs(y)) / 2.0;
}


inline precision minmod3(precision x, precision y, precision z)
{
   return minmod(x, minmod(y, z));
}


inline precision spatial_derivative(precision qm, precision q, precision qp)
{
	//return minmod3(THETA * (q - qm), (qp - qm) / 2.0, THETA * (qp - q));
	return (qp - qm) / 2.0;
}




void source_terms(precision * const __restrict__ S, const precision * const __restrict__ q, precision e, precision t, const precision * const __restrict__ qi1, const precision * const __restrict__ qj1, const precision * const __restrict__ qk1, const precision * const __restrict__ e1, const precision * const __restrict__ ui1, const precision * const __restrict__ uj1, const precision * const __restrict__ uk1, precision ut, precision ux, precision uy, precision un, precision ut_p, precision ux_p, precision uy_p, precision un_p, precision dt, precision dx, precision dy, precision dn, precision etabar)
{
	precision t2 = t * t;		// useful expressions
	precision t3 = t2 * t;

	precision ut2 = ut * ut;
	precision un2 = un * un;


	// conserved variables
	precision ttt = q[0];
	precision ttx = q[1];
	precision tty = q[2];
	precision ttn = q[3];

	precision pl  = q[4];

#if (PT_MATCHING == 1)
	precision pt  = q[5];
#else
	precision pt  = 0.5 * (e - pl);
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

	// fluid velocity derivatives
	precision ut_sim = ui1[0];		// ut [i-1, i+1]
	precision ut_sip = ui1[1];
	precision ux_sim = ui1[2];		// ux [i-1, i+1]
	precision ux_sip = ui1[3];
	precision uy_sim = ui1[4];		// uy [i-1, i+1]
	precision uy_sip = ui1[5];
	precision un_sim = ui1[6];		// un [i-1, i+1]
	precision un_sip = ui1[7];

	precision dut_dx = spatial_derivative(ut_sim, ut, ut_sip) / dx;
	precision dux_dx = spatial_derivative(ux_sim, ux, ux_sip) / dx;
	precision duy_dx = spatial_derivative(uy_sim, uy, uy_sip) / dx;
	precision dun_dx = spatial_derivative(un_sim, un, un_sip) / dx;

	precision ut_sjm = uj1[0];		// ut [j-1, j+1]
	precision ut_sjp = uj1[1];
	precision ux_sjm = uj1[2];		// ux [j-1, j+1]
	precision ux_sjp = uj1[3];
	precision uy_sjm = uj1[4];		// uy [j-1, j+1]
	precision uy_sjp = uj1[5];
	precision un_sjm = uj1[6];		// un [j-1, j+1]
	precision un_sjp = uj1[7];

	precision dut_dy = spatial_derivative(ut_sjm, ut, ut_sjp) / dy;
	precision dux_dy = spatial_derivative(ux_sjm, ux, ux_sjp) / dy;
	precision duy_dy = spatial_derivative(uy_sjm, uy, uy_sjp) / dy;
	precision dun_dy = spatial_derivative(un_sjm, un, un_sjp) / dy;

	precision ut_skm = uk1[0];		// ut [k-1, k+1]
	precision ut_skp = uk1[1];
	precision ux_skm = uk1[2];		// ux [k-1, k+1]
	precision ux_skp = uk1[3];
	precision uy_skm = uk1[4];		// uy [k-1, k+1]
	precision uy_skp = uk1[5];
	precision un_skm = uk1[6];		// un [k-1, k+1]
	precision un_skp = uk1[7];

	precision dut_dn = spatial_derivative(ut_skm, ut, ut_skp) / dn;
	precision dux_dn = spatial_derivative(ux_skm, ux, ux_skp) / dn;
	precision duy_dn = spatial_derivative(uy_skm, uy, uy_skp) / dn;
	precision dun_dn = spatial_derivative(un_skm, un, un_skp) / dn;


	// spatial velocity components
	precision vx = ux / ut;
	precision vy = uy / ut;
	precision vn = un / ut;

	// spatial velocity derivatives
	precision dvx_dn = (dux_dn  -  vx * dut_dn) / ut;
	precision dvy_dn = (duy_dn  -  vy * dut_dn) / ut;

	// divergence of v
	precision div_v = (dux_dx  +  duy_dy  +  dun_dn  -  vx * dut_dx  -  vy * dut_dy  -  vn * dut_dn) / ut;


	// time derivatives of u
	precision dut_dt = (ut - ut_p) / dt;
	precision dux_dt = (ux - ux_p) / dt;
	precision duy_dt = (uy - uy_p) / dt;
	precision dun_dt = (un - un_p) / dt;

	// radial velocity derivatives
	precision duT2_dx = ux * dux_dx  +  uy * duy_dx;
	precision duT2_dy = ux * dux_dy  +  uy * duy_dy;
	precision duT2_dn = ux * dux_dn  +  uy * duy_dn;

	precision utperp  = sqrt(1.0  +  ux * ux  +  uy * uy);
	precision utperp2 = utperp * utperp;

	// longitudinal basis vector
	precision zt = t * un / utperp;
	precision zn = ut / t / utperp;

	precision zt2  = zt * zt;
	precision ztzn = zt * zn;
	precision zn2  = zn * zn;

	// longitudinal vector derivatives
	precision dzt_dx = t * (dun_dx  -  un * duT2_dx / utperp2) / utperp;
	precision dzt_dy = t * (dun_dy  -  un * duT2_dy / utperp2) / utperp;
	precision dzt_dn = t * (dun_dn  -  un * duT2_dn / utperp2) / utperp;

	precision dzn_dx = (dut_dx  -  ut * duT2_dx / utperp2) / (t * utperp);
	precision dzn_dy = (dut_dy  -  ut * duT2_dy / utperp2) / (t * utperp);
	precision dzn_dn = (dut_dn  -  ut * duT2_dn / utperp2) / (t * utperp);

	// scalar expansion rate
	precision theta = dut_dt  +  dux_dx  +  duy_dy  +  dun_dn  +  ut / t;

	// longitudinal and transverse expansion rates
	precision thetaL = - zt2 * dut_dt  +  t2 * zn2 * dun_dn  +  ztzn * (t2 * dun_dt  -  dut_dn)  +  t * zn2 * ut;
	precision thetaT = theta - thetaL;


	// relaxation times
	precision peq = equilibriumPressure(e);
	precision T   = effectiveTemperature(e);

	precision taupiInv = 0.2 * T / etabar;

	// anisotropic transport coefficients
	transport_coefficients aniso;
	aniso.compute_transport_coefficients(e, pl, pt);

	// pl coefficients
	precision zeta_LL = aniso.I_240  -  3.0 * pl;
	precision zeta_TL = aniso.I_221  -  pl;

	// pt coefficients
#if (PT_MATCHING == 1)
	precision zeta_LT = aniso.I_221  -  pt;
	precision zeta_TT = 2.0 * (aniso.I_202 - pt);
#endif

	// L^munu components and derivatives
	precision dp  = pl - pt;

	precision Ltt = dp * zt2;
	precision Ltn = dp * ztzn;
	//precision Lnn = dp * zn2;

	precision tnn = (e + pt) * un2  +  pt / t2  +  dp * zn2;

	precision dLtt_dx = (dpl_dx - dpt_dx) * zt2  +  2.0 * dp * zt * dzt_dx;
	precision dLtt_dy = (dpl_dy - dpt_dy) * zt2  +  2.0 * dp * zt * dzt_dy;
	precision dLtt_dn = (dpl_dn - dpt_dn) * zt2  +  2.0 * dp * zt * dzt_dn;

	precision dLtn_dx = (dpl_dx - dpt_dx) * ztzn  +  dp * (dzt_dx * zn  +  dzn_dx * zt);
	precision dLtn_dy = (dpl_dy - dpt_dy) * ztzn  +  dp * (dzt_dy * zn  +  dzn_dy * zt);
	precision dLtn_dn = (dpl_dn - dpt_dn) * ztzn  +  dp * (dzt_dn * zn  +  dzn_dn * zt);

	precision dLnn_dn = (dpl_dn - dpt_dn) * zn2  +  2.0 * dp * zn * dzn_dn;


	// conservation laws
	S[0] = - (ttt / t  +  t * tnn)  +  div_v * (Ltt - pt)  +  vx * (dLtt_dx - dpt_dx)  +  vy * (dLtt_dy - dpt_dy)  +  vn * (dLtt_dn - dpt_dn)  - dLtn_dn;
	S[1] = - ttx / t  -  dpt_dx  -  0 * 0.5 * (vx * dLtn_dn  -  Ltn * dvx_dn);
	S[2] = - tty / t  -  dpt_dy  -  0 * 0.5 * (vy * dLtn_dn  -  Ltn * dvy_dn);
	S[3] = - 3.0 * ttn / t  -  dpt_dn / t2  +  div_v * Ltn  +  vx * dLtn_dx  +  vy * dLtn_dy  +  vn * dLtn_dn  -  dLnn_dn;


	// pl relaxation equation
	precision dpl = - dp * taupiInv / 1.5  +  zeta_LL * thetaL  +  zeta_TL * thetaT;
	S[4] = dpl / ut  +  div_v * pl;


#if (PT_MATCHING == 1)
	precision dpt =  dp * taupiInv / 3.0  +  zeta_LT * thetaL  +  zeta_TT * thetaT;
	S[5] = dpt / ut  +  div_v * pt;
#endif

}




