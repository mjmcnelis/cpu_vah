/*
 * SemiDiscreteKurganovTadmorScheme.cpp
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#include <stdlib.h>
#include <math.h>
#include "../include/FluxTerms.h"
#include "../include/DynamicalVariables.h"


inline PRECISION sign(PRECISION x)
{
	if(x < 0.0) return -1.0;
	else return 1.0;
}

inline PRECISION minmod(PRECISION x, PRECISION y)
{
	return (sign(x) + sign(y)) * fmin(fabs(x), fabs(y)) / 2.0;
}


inline PRECISION minmod3(PRECISION x, PRECISION y, PRECISION z)
{
   return minmod(x, minmod(y, z));
}


PRECISION approx_derivative(PRECISION qm, PRECISION q, PRECISION qp)
{
	return minmod3(THETA * (q - qm), (qp - qm) / 2.0, THETA * (qp - q));
}


PRECISION right_half_cell_extrapolation_forward(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp)
{

	return qp - approx_derivative(q, qp, qpp) / 2.0;	// Eq. (63)
}


PRECISION left_half_cell_extrapolation_forward(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp)
{
	return q + approx_derivative(qm, q, qp) / 2.0;		// Eq. (64)
}


PRECISION right_half_cell_extrapolation_backward(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp)
{
	return q - approx_derivative(qm, q, qp) / 2.0;		// Eq. (65)
}


PRECISION left_half_cell_extrapolation_backward(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp)
{
	return qm + approx_derivative(qmm, qm, q) / 2.0;	// Eq. (66)
}


void flux_terms(PRECISION * const __restrict__ H, const PRECISION * const __restrict__ Q_data, const PRECISION * const __restrict__ Q1_data, const PRECISION * const __restrict__ Q2_data, const PRECISION * const __restrict__ V_data, PRECISION u_i, PRECISION ut,
	PRECISION (* const right_half_cell_extrapolation)(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp),
	PRECISION (* const left_half_cell_extrapolation)(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp))
{
	// neighbor fluid velocities
	PRECISION vmm = V_data[0];
	PRECISION vm  = V_data[1];
	PRECISION v   = u_i / ut;
	PRECISION vp  = V_data[2];
	PRECISION vpp = V_data[3];

	// left / right extrapolated speeds
	PRECISION vL = fabs(left_half_cell_extrapolation(vmm, vm, v, vp, vpp));
	PRECISION vR = fabs(right_half_cell_extrapolation(vmm, vm, v, vp, vpp));

	// local propagation speed
	PRECISION a = fmax(vL, vR);

	// left / right extrapolated values of q and F = vq
	PRECISION qmm, qm, q, qp, qpp;
	PRECISION Fmm, Fm, F, Fp, Fpp;
	PRECISION qL, qR, FL, FR;
	int p = 0;

	for(short n = 0; n < NUMBER_CONSERVED_VARIABLES; n++)
	{
		q 	= Q_data[n];

		qm 	= Q1_data[p];
		qp 	= Q1_data[p + 1];

		qmm = Q2_data[p];
		qpp = Q2_data[p + 1];

		Fmm = qmm * vmm;
		Fm 	= qm * vm;
		F 	= q * v;
		Fp 	= qp * vp;
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

