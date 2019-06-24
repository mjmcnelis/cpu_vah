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


inline precision approximate_derivative(precision qm, precision q, precision qp)
{
	//return minmod3(THETA * (q - qm), (qp - qm) / 2.0, THETA * (qp - q));
	return (qp - qm) / 2.0;
}


void source_terms(precision * const __restrict__ S, const precision * const __restrict__ ql, precision e_s, precision t, const precision * const __restrict__ qi1, const precision * const __restrict__ qj1, const precision * const __restrict__ qk1, const precision * const __restrict__ e1, const precision * const __restrict__ ui1, const precision * const __restrict__ uj1, const precision * const __restrict__ uk1, precision ut, precision ux, precision uy, precision un, precision ut_p, precision ux_p, precision uy_p, precision un_p, precision dt, precision dx, precision dy, precision dn, precision etabar)
{
	precision t2 = t * t;
	precision t3 = t2 * t;

	// conserved variables
	precision ttt = ql[0];
	precision ttx = ql[1];
	precision tty = ql[2];
	precision ttn = ql[3];

	precision pl = ql[4];
	precision pt = ql[5];


	pt = 0.5 * (e_s - pl);	// temporary


	// fluid velocity
	precision ut_sim = ui1[0];		// ut [i-1, i+1]
	precision ut_sip = ui1[1];

	precision ut_sjm = uj1[0];		// ut [j-1, j+1]
	precision ut_sjp = uj1[1];

	precision ut_skm = uk1[0];		// ut [k-1, k+1]
	precision ut_skp = uk1[1];


	precision ux_sim = ui1[2];		// ux [i-1, i+1]
	precision ux_sip = ui1[3];

	precision ux_sjm = uj1[2];		// ux [j-1, j+1]
	precision ux_sjp = uj1[3];

	precision ux_skm = uk1[2];		// ux [k-1, k+1]
	precision ux_skp = uk1[3];


	precision uy_sim = ui1[4];		// uy [i-1, i+1]
	precision uy_sip = ui1[5];

	precision uy_sjm = uj1[4];		// uy [j-1, j+1]
	precision uy_sjp = uj1[5];

	precision uy_skm = uk1[4];		// uy [k-1, k+1]
	precision uy_skp = uk1[5];


	precision un_sim = ui1[6];		// un [i-1, i+1]
	precision un_sip = ui1[7];

	precision un_sjm = uj1[6];		// un [j-1, j+1]
	precision un_sjp = uj1[7];

	precision un_skm = uk1[6];		// un [k-1, k+1]
	precision un_skp = uk1[7];



	precision dut_dx = approximate_derivative(ut_sim, ut, ut_sip) / dx;
	precision dut_dy = approximate_derivative(ut_sjm, ut, ut_sjp) / dy;
	precision dut_dn = approximate_derivative(ut_skm, ut, ut_skp) / dn;

	precision dux_dx = approximate_derivative(ux_sim, ux, ux_sip) / dx;
	precision dux_dy = approximate_derivative(ux_sjm, ux, ux_sjp) / dy;
	precision dux_dn = approximate_derivative(ux_skm, ux, ux_skp) / dn;

	precision duy_dx = approximate_derivative(uy_sim, uy, uy_sip) / dx;
	precision duy_dy = approximate_derivative(uy_sjm, uy, uy_sjp) / dy;
	precision duy_dn = approximate_derivative(uy_skm, uy, uy_skp) / dn;

	precision dun_dx = approximate_derivative(un_sim, un, un_sip) / dx;
	precision dun_dy = approximate_derivative(un_sjm, un, un_sjp) / dy;
	precision dun_dn = approximate_derivative(un_skm, un, un_skp) / dn;




	// primary variables
	precision e_sim = e1[0];		// e [i-1, i+1]
	precision e_sip = e1[1];

	precision e_sjm = e1[2];		// e [j-1, j+1]
	precision e_sjp = e1[3];

	precision e_skm = e1[4];		// e [k-1, k+1]
	precision e_skp = e1[5];


	// primary variable spatial derivatives
	precision de_dx = approximate_derivative(e_sim, e_s, e_sip) / dx;		// after put in pt matching won't need this either
	precision de_dy = approximate_derivative(e_sjm, e_s, e_sjp) / dy;
	precision de_dn = approximate_derivative(e_skm, e_s, e_skp) / dn;




	// longitudinal pressure
	// conserved variables
	int n = 8;						// start at pl since neighbors of ttt, ttx, tty, ttn aren't needed here

	precision pl_sim = qi1[n];		// pl [i-1, i+1]	(sim = i-1)
	precision pl_sip = qi1[n + 1];	//				  	(sip = i+1)

	precision pl_sjm = qj1[n];		// pl [j-1, j+1]	(sjm = j-1)
	precision pl_sjp = qj1[n + 1];  //				 	(sjp = j+1)

	precision pl_skm = qk1[n];		// pl [k-1, k+1]  	(skm = k-1)
	precision pl_skp = qk1[n + 1];	//					(skp = k+1)

	precision dpl_dx = approximate_derivative(pl_sim, pl, pl_sip) / dx;
	precision dpl_dy = approximate_derivative(pl_sjm, pl, pl_sjp) / dy;
	precision dpl_dn = approximate_derivative(pl_skm, pl, pl_skp) / dn;

	precision dpt_dx = 0.5 * (de_dx - dpl_dx);	// temporary
	precision dpt_dy = 0.5 * (de_dy - dpl_dy);
	precision dpt_dn = 0.5 * (de_dn - dpl_dn);


	// spatial velocity components
	precision vx = ux / ut;
	precision vy = uy / ut;
	precision vn = un / ut;

	// divergence of v
	precision div_v = (dux_dx  +  duy_dy  +  dun_dn  -  vx * dut_dx  -  vy * dut_dy  -  vn * dut_dn) / ut;


	precision un2 = un * un;
	precision ut2 = ut * ut;


	// time derivatives of u
	precision dtut = (ut - ut_p) / dt;
	precision dtux = (ux - ux_p) / dt;
	precision dtuy = (uy - uy_p) / dt;
	precision dtun = (un - un_p) / dt;


	precision dut = ut * dtut  +  ux * dut_dx  +  uy * dut_dy  +  un * dut_dn;
	precision dux = ut * dtux  +  ux * dux_dx  +  uy * dux_dy  +  un * dux_dn;
	precision duy = ut * dtuy  +  ux * duy_dx  +  uy * duy_dy  +  un * duy_dn;
	precision dun = ut * dtun  +  ux * dun_dx  +  uy * dun_dy  +  un * dun_dn;






	// expansion rate
	precision theta = dtut  +  dux_dx  +  duy_dy  +  dun_dn  +  ut / t;
	precision theta_over_3 = theta / 3.0;

	// shear tensor
	precision stt = - t * ut * un2  +  (dtut  -  ut * dut)  +  (ut2 - 1.0) * theta_over_3;
	// precision stx = -(t * un2 * ux) / 2 + (dtux - dut_dx) / 2 - (ux * dut + ut * dux) / 2 + ut * ux * theta / 3;
	// precision sty = -(t * un2 * uy) / 2 + (dtuy - dut_dy) / 2 - (uy * dut + ut * duy) / 2 + ut * uy * theta / 3;
	precision stn = - un * (2.0 * ut2  +  t2 * un2) / (2.0 * t)  +  (dtun  -  dut_dn / t2) / 2.0  -  (un * dut  +  ut * dun) / 2.0  +  ut * un * theta_over_3;
	// precision sxx = -(dux_dx + ux * dux) + (1 + ux*ux) * theta / 3;
	// precision sxy = -(duy_dx + dux_dy) / 2 - (uy * dux + ux * duy) / 2	+ ux * uy * theta / 3;
	// precision sxn = -ut * ux * un / t - (dun_dx + dux_dn / t2) / 2 - (un * dux + ux * dun) / 2 + ux * un * theta / 3;
	// precision syy = -(duy_dy + uy * duy) + (1 + uy*uy) * theta / 3;
	// precision syn = -ut * uy * un / t - (dun_dy + duy_dn / t2) / 2 - (un * duy + uy * dun) / 2 + uy * un * theta / 3;
	precision snn = - ut * (1.0  +  2.0 * t2 * un2) / t3  -  dun_dn / t2  -  un * dun + (1.0 / t2  +  un2) * theta_over_3;





	// transverse flow velocity
	precision uT2 = ux * ux  +  uy * uy;
	precision uT = sqrt(uT2);

	precision uTdxuT = ux * dux_dx  +  uy * duy_dx;
	precision uTdyuT = ux * dux_dy  +  uy * duy_dy;
	precision uTdnuT = ux * dux_dn  +  uy * duy_dn;


	precision F = 1.0 + uT2;
	precision F2 = F * F;

	precision utperp  = sqrt(1.0  +  ux * ux  +  uy * uy);
	precision utperp2 = utperp * utperp;


	// longitudinal basis vector
	precision zt = t * un / utperp;
	precision zn = ut / t / utperp;

	precision zt2  = zt * zt;
	precision ztzn = zt * zn;
	precision zn2  = zn * zn;

	precision A = zt2 * stt  -  2 * t2 * ztzn * stn  +  t2 * t2 * zn2 * snn;

	// longitudinal and transverse expansion rates
	precision thetaT = 2.0 * theta_over_3  +  A;
	precision thetaL = theta - thetaT;


	// anisotropic transport coefficients
	if(e_s == 0) e_s = 1.e-7;
	precision p = equilibriumPressure(e_s);
	precision T = effectiveTemperature(e_s);

	precision taupiInv = 0.2 * T / etabar;


	transport_coefficients aniso;

	aniso.compute_transport_coefficients(e_s, pl, pt);

	precision zeta_zz = aniso.I_240  -  3.0 * pl;
	precision zeta_zT = aniso.I_221  -  pl;				// okay now it's working

	
	// precision a  = pl / e_s;
	// precision a2 = a * a;
	// precision a3 = a2 * a;
	// precision a4 = a3 * a;
	// precision a5 = a4 * a;
	// precision a6 = a5 * a;
	// precision a7 = a6 * a;

	// precision Rtilde = (-6.674731906076046e-6 + 0.004617789933500251*a + 0.7207562721999754*a2 + 9.097427250602184*a3 - 4.475814747302824*a4 - 36.37501529319408*a5 +
 //     46.868405146729316*a6 - 15.833867583743228*a7)/
 //   (0.06856675185266 + 2.9181587012768597*a + 11.951184087839218*a2 - 29.708257843442173*a3 - 2.618233802059826*a4 + 34.646239784689065*a5 -
 //     19.62596366454439*a6 + 2.374808442453899*a7);

	// precision zeta_zz = Rtilde * e_s  -  3.0 * pl;
	// precision zeta_zT = (Rtilde * e_s +  pl) / 2.0;
	



	precision dp = pl - pt;

	precision Ltt = dp * zt2;
	precision Ltn = dp * ztzn;
	precision Lnn = dp * zn2;


	precision dzt_dx = t * (dun_dx  -  un * uTdxuT / utperp2) / utperp;
	precision dzt_dy = t * (dun_dy  -  un * uTdyuT / utperp2) / utperp;
	precision dzt_dn = t * (dun_dn  -  un * uTdnuT / utperp2) / utperp;

	precision dzn_dx = (dut_dx  -  ut * uTdxuT / utperp2) / (t * utperp);
	precision dzn_dy = (dut_dy  -  ut * uTdyuT / utperp2) / (t * utperp);
	precision dzn_dn = (dut_dn  -  ut * uTdnuT / utperp2) / (t * utperp);

	precision dLtt_dx = (dpl_dx - dpt_dx) * zt2  +  2.0 * dp * zt * dzt_dx;
	precision dLtt_dy = (dpl_dy - dpt_dy) * zt2  +  2.0 * dp * zt * dzt_dy;
	precision dLtt_dn = (dpl_dn - dpt_dn) * zt2  +  2.0 * dp * zt * dzt_dn;

	precision dLtn_dx = (dpl_dx - dpt_dx) * ztzn  +  dp * (dzt_dx * zn  +  dzn_dx * zt);
	precision dLtn_dy = (dpl_dy - dpt_dy) * ztzn  +  dp * (dzt_dy * zn  +  dzn_dy * zt);
	precision dLtn_dn = (dpl_dn - dpt_dn) * ztzn  +  dp * (dzt_dn * zn  +  dzn_dn * zt);

	precision dLnn_dn = (dpl_dn - dpt_dn) * zn2  +  2.0 * dp * zn * dzn_dn;




	// source terms
	precision tnn = (e_s + pt) * un2  +  pt / t2  +  Lnn;
	//precision dpl = - (pl - p) * taupiInv  +  zeta_zz * thetaL  -  zeta_zT * thetaT;
	precision dpl = - (pl - p) * taupiInv  +  zeta_zz * thetaL  +  zeta_zT * thetaT;

	// conservation laws
	S[0] = - (ttt / t  +  t * tnn)  +  div_v * (Ltt - pt)  -  vx * dpt_dx  -  vy * dpt_dy  -  vn * dpt_dn  +  vx * dLtt_dx  +  vy * dLtt_dy  +  vn * dLtt_dn - dLtn_dn;
	S[1] = - ttx / t  -  dpt_dx;
	S[2] = - tty / t  -  dpt_dy;
	S[3] = - 3.0 * ttn / t  -  dpt_dn / t2  +  div_v * Ltn  +  vx * dLtn_dx  +  vy * dLtn_dy  +  2.0 * vn * dLtn_dn  -  dLnn_dn;

	// relaxation equations
	S[4] = dpl / ut  +  div_v * pl;

}




