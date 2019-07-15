
#include <stdlib.h>
#include <math.h>
#include "../include/FluxTerms.h"
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"



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


precision approx_derivative(precision qm, precision q, precision qp)
{
	return minmod3(THETA * (q - qm), (qp - qm) / 2.0, THETA * (qp - q));
}


precision right_half_cell_extrapolation_forward(precision qmm, precision qm, precision q, precision qp, precision qpp)
{

	return qp - approx_derivative(q, qp, qpp) / 2.0;	// Eq. (63)
}


precision left_half_cell_extrapolation_forward(precision qmm, precision qm, precision q, precision qp, precision qpp)
{
	return q + approx_derivative(qm, q, qp) / 2.0;		// Eq. (64)
}


precision right_half_cell_extrapolation_backward(precision qmm, precision qm, precision q, precision qp, precision qpp)
{
	return q - approx_derivative(qm, q, qp) / 2.0;		// Eq. (65)
}


precision left_half_cell_extrapolation_backward(precision qmm, precision qm, precision q, precision qp, precision qpp)
{
	return qm + approx_derivative(qmm, qm, q) / 2.0;	// Eq. (66)
}


void flux_terms_new(precision * const __restrict__ Hp, precision * const __restrict__ Hm, const precision * const __restrict__ q_data, const precision * const __restrict__ q1_data, const precision * const __restrict__ q2_data, const precision * const __restrict__ ui_data, const precision * const __restrict__ ut_data, precision ui, precision ut)
{
	// neighbor fluid velocities
	precision uimm = ui_data[0];	// spatial component i
	precision uim  = ui_data[1];
	precision uip  = ui_data[2];
	precision uipp = ui_data[3];

	precision utmm = ut_data[0];	// time component
	precision utm  = ut_data[1];
	precision utp  = ut_data[2];
	precision utpp = ut_data[3];

	// neigbhor spatial velocities
	precision vmm = uimm / utmm;
	precision vm  = uim  / utm;
	precision v   = ui   / ut;
	precision vp  = uip  / utp;
	precision vpp = uipp / utpp;


	// precision dvp, dv, dvm;
	// spatial velocity derivatives (not as good)
	// if(ut > 1.21061)
	// {
	// 	precision dui  = approx_derivative(uim, ui, uip);
	// 	precision dut  = approx_derivative(utm, ut, utp);
	// 	dv  = (dui   -   v * dut)  / ut;
	// }
	// else
	// {
	// 	dv  = approx_derivative(vm, v, vp);
	// }
	// if(utp > 1.21061)
	// {
	// 	precision duip = approx_derivative(ui, uip, uipp);
	// 	precision dutp = approx_derivative(ut, utp, utpp);
	// 	dvp = (duip  -  vp * dutp) / utp;
	// }
	// else
	// {
	// 	dvp = approx_derivative(v, vp, vpp);
	// }
	// if(utm > 1.21061)
	// {
	// 	precision duim = approx_derivative(uimm, uim, ui);
	// 	precision dutm = approx_derivative(utmm, utm, ut);	
	// 	dvm = (duim  -  vm * dutm) / utm;
	// }
	// else
	// {
	// 	dvm = approx_derivative(vmm, vm, v);
	// }
	
	
	precision dvp = approx_derivative(v, vp, vpp);	
	precision dv  = approx_derivative(vm, v, vp);
	precision dvm = approx_derivative(vmm, vm, v);

	// extrapolated spatial velocities
	precision vRp = vp  -  0.5 * dvp;	// q^{+}_{i+1/2}	
	precision vLp = v   +  0.5 * dv;	// q^{-}_{i+1/2}	
	precision vRm = v   -  0.5 * dv;	// q^{+}_{i-1/2}
	precision vLm = vm  +  0.5 * dvm;	// q^{-}_{i-1/2}

	// local propagation speeds
	precision ap = fmax(fabs(vLp), fabs(vRp));
	precision am = fmax(fabs(vLm), fabs(vRm));

	int p = 0;

	// compute the flux terms Hp, Hm
	for(int n = 0; n < NUMBER_CONSERVED_VARIABLES; n++)
	{
		// neighbor conserved variables
		precision q   =  q_data[n];
		precision qm  = q1_data[p];
		precision qp  = q1_data[p + 1];
		precision qmm = q2_data[p];
		precision qpp = q2_data[p + 1];

		// conserved variable derivatives
		precision dqp = approx_derivative(q,   qp, qpp);
		precision dq  = approx_derivative(qm,  q,  qp);
		precision dqm = approx_derivative(qmm, qm, q);

		// extrapolated conserved variables
		precision qRp = qp  -  0.5 * dqp;	// q^{+}_{i+1/2}	
		precision qLp = q   +  0.5 * dq;	// q^{-}_{i+1/2}	
		precision qRm = q   -  0.5 * dq;	// q^{+}_{i-1/2}
		precision qLm = qm  +  0.5 * dqm;	// q^{-}_{i-1/2}

		// neighbor fluxes
		precision Fm  = qm * vm;
		precision F   = q  * v;      
		precision Fp  = qp * vp;

		// extrapolated fluxes
		precision FRp = Fp  -  0.5 * (qp * dvp  +  vp * dqp);	// F^{+}_{i+1/2}	(R = +, p = i + 1/2)	Eq. (63)
		precision FLp = F   +  0.5 * (q  * dv   +  v  * dq);	// F^{-}_{i+1/2}	(L = -, p = i + 1/2)	Eq. (64)
		precision FRm = F   -  0.5 * (q  * dv   +  v  * dq);	// F^{+}_{i-1/2}	(R = -, m = i - 1/2)	Eq. (65)
		precision FLm = Fm  +  0.5 * (qm * dvm  +  vm * dqm);	// F^{-}_{i-1/2}	(L = -, m = i - 1/2)	Eq. (66)

		// precision FRp = qRp * vRp;
		// precision FLp = qLp * vLp;
		// precision FRm = qRm * vRm;
		// precision FLm = qLm * vLm;

		// Hp, Hm from Eq.(61)
		Hp[n] = (FRp  +  FLp  -  ap * (qRp - qLp)) / 2.0;
		Hm[n] = (FRm  +  FLm  -  am * (qRm - qLm)) / 2.0;

		p += 2;
	}
}


void flux_terms(precision * const __restrict__ H, const precision * const __restrict__ q_data, const precision * const __restrict__ q1_data, const precision * const __restrict__ q2_data, const precision * const __restrict__ ui_data, const precision * const __restrict__ ut_data, precision ui, precision ut,
	precision (* const right_half_cell_extrapolation)(precision qmm, precision qm, precision q, precision qp, precision qpp),
	precision (* const left_half_cell_extrapolation)(precision qmm, precision qm, precision q, precision qp, precision qpp))
{
	// neighbor fluid velocities
	precision uimm = ui_data[0];	// spatial component i
	precision uim  = ui_data[1];
	precision uip  = ui_data[2];
	precision uipp = ui_data[3];

	precision utmm = ut_data[0];	// time component
	precision utm  = ut_data[1];
	precision utp  = ut_data[2];
	precision utpp = ut_data[3];


	// left / right extrapolated velocities
	// precision uiL = left_half_cell_extrapolation(uimm, uim, ui, uip, uipp);
	// precision uiR = right_half_cell_extrapolation(uimm, uim, ui, uip, uipp);

	// precision utL = left_half_cell_extrapolation(utmm, utm, ut, utp, utpp);
	// precision utR = right_half_cell_extrapolation(utmm, utm, ut, utp, utpp);

	// precision vL = uiL / utL;
	// precision vR = uiR / utR;

	precision vL =  left_half_cell_extrapolation(uimm / utmm, uim / utm, ui / ut, uip / utp, uipp / utpp);
	precision vR = right_half_cell_extrapolation(uimm / utmm, uim / utm, ui / ut, uip / utp, uipp / utpp);


	// local propagation speed
	precision a = fmax(fabs(vL), fabs(vR));

	// left / right extrapolated values of q and F = vq
	precision qmm, qm, q, qp, qpp;
	precision Fmm, Fm, F, Fp, Fpp;
	precision qL, qR, FL, FR;
	int p = 0;

	for(int n = 0; n < NUMBER_CONSERVED_VARIABLES; n++)
	{
		q 	=  q_data[n];

		qm 	= q1_data[p];
		qp 	= q1_data[p + 1];
		
		qmm = q2_data[p];
		qpp = q2_data[p + 1];

		Fmm = qmm * uimm / utmm;
		Fm 	= qm  * uim  / utm;
		F 	= q   * ui   / ut;
		Fp 	= qp  * uip  / utp;
		Fpp = qpp * uipp / utpp;

		// left / right extrapolated velocities
		qL = left_half_cell_extrapolation(qmm, qm, q, qp, qpp);
		qR = right_half_cell_extrapolation(qmm, qm, q, qp, qpp);

		FL = left_half_cell_extrapolation(Fmm, Fm, F, Fp, Fpp);
		FR = right_half_cell_extrapolation(Fmm, Fm, F, Fp, Fpp);

		// FL = qL * vL;
		// FR = qR * vR;

		// H from Eq.(61)
		H[n] = (FR  +  FL  -  a * (qR - qL)) / 2.0;

		p += 2;
	}
}


/* // old version 
void flux_terms(precision * const __restrict__ H, const precision * const __restrict__ Q_data, const precision * const __restrict__ Q1_data, const precision * const __restrict__ Q2_data, const precision * const __restrict__ V_data, precision u_i, precision ut,
	precision (* const right_half_cell_extrapolation)(precision qmm, precision qm, precision q, precision qp, precision qpp),
	precision (* const left_half_cell_extrapolation)(precision qmm, precision qm, precision q, precision qp, precision qpp))
{
	// neighbor fluid velocities
	precision vmm = V_data[0];
	precision vm  = V_data[1];
	precision v   = u_i / ut;
	precision vp  = V_data[2];
	precision vpp = V_data[3];

	// left / right extrapolated speeds
	precision vL = fabs(left_half_cell_extrapolation(vmm, vm, v, vp, vpp));
	precision vR = fabs(right_half_cell_extrapolation(vmm, vm, v, vp, vpp));

	// local propagation speed
	precision a = fmax(vL, vR);

	// left / right extrapolated values of q and F = vq
	precision qmm, qm, q, qp, qpp;
	precision Fmm, Fm, F, Fp, Fpp;
	precision qL, qR, FL, FR;
	int p = 0;

	for(int n = 0; n < NUMBER_CONSERVED_VARIABLES; n++)
	{
		q 	= Q_data[n];

		qm 	= Q1_data[p];
		qp 	= Q1_data[p + 1];

		qmm = Q2_data[p];
		qpp = Q2_data[p + 1];

		Fmm = qmm * vmm;
		Fm 	= qm  * vm;
		F 	= q   * v;
		Fp 	= qp  * vp;
		Fpp = qpp * vpp;

		qL = left_half_cell_extrapolation(qmm, qm, q, qp, qpp);
		qR = right_half_cell_extrapolation(qmm, qm, q, qp, qpp);

		FL = left_half_cell_extrapolation(Fmm, Fm, F, Fp, Fpp);
		FR = right_half_cell_extrapolation(Fmm, Fm, F, Fp, Fpp);

		// H from Eq.(61)
		H[n] = (FR + FL  -  a * (qR - qL)) / 2.0;

		p += 2;
	}
}
*/

