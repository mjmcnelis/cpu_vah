
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

