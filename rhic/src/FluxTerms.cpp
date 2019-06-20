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
#include "../include/EnergyMomentumTensor.h"


// 6/12 moving all the KT-algorithm stuff in one file



// the purpose of inline:
//  - overhead can occur if the function call takes longer than the function's execution time
// 	- this is likely to occur only for simple functions that execute a small amount of code and are called frequently
//	- inline inserts the function code at the point of the function call (a sort of substitution)


inline double sign(PRECISION x)
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


PRECISION approximateDerivative(PRECISION qm, PRECISION q, PRECISION qp)
{
	// qm = q_{i-1}		|	q = q_i		|	qp = q_{i+1}

	PRECISION backward_derivative = THETA * (q - qm);	// with THETA flux factor
	PRECISION central_derivative = (qp - qm) / 2.0;
	PRECISION forward_derivative = THETA * (qp - q);	// with THETA flux factor

	return minmod3(backward_derivative, central_derivative, forward_derivative);
}


PRECISION rightHalfCellExtrapolationForward(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp)
{

	return qp - approximateDerivative(q, qp, qpp) / 2.0;	// Eq. (63)
}


PRECISION leftHalfCellExtrapolationForward(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp)
{
	return q + approximateDerivative(qm, q, qp) / 2.0;		// Eq. (64)
}


PRECISION rightHalfCellExtrapolationBackwards(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp)
{
	return q - approximateDerivative(qm, q, qp) / 2.0;		// Eq. (65)
}


PRECISION leftHalfCellExtrapolationBackwards(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp)
{
	return qm + approximateDerivative(qmm, qm, q) / 2.0;	// Eq. (66)
}


PRECISION spectralRadiusX(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un)
{
	return fabs(ux / ut);	// can I inline this?
}

PRECISION spectralRadiusY(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un)
{
	return fabs(uy / ut);
}

PRECISION spectralRadiusZ(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un)
{
	return fabs(un / ut);
}


inline PRECISION localPropagationSpeed(PRECISION utr, PRECISION uxr, PRECISION uyr, PRECISION unr, PRECISION utl, PRECISION uxl, PRECISION uyl, PRECISION unl, PRECISION (*spectralRadius)(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un))
{
	// the local propagation speed in the KT algorithm

	// maximal local speed at the cell boundaries x_{j +/- 1/2}
	// spectralRadius is one of the inputs, so either X, Y, or Z is passed

	PRECISION rhoLeftMovingWave = spectralRadius(utl, uxl, uyl, unl);
	PRECISION rhoRightMovingWave = spectralRadius(utr, uxr, uyr, unr);

	return fmax(rhoLeftMovingWave, rhoRightMovingWave);
}


PRECISION Fx(PRECISION q, PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un)
{
	return ux * q / ut;
}


PRECISION Fy(PRECISION q, PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un)
{
	return uy * q / ut;
}


PRECISION Fz(PRECISION q, PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un)
{
	return un * q / ut;
}


void flux_terms(const PRECISION * const __restrict__ Q_data, const PRECISION * const __restrict__ Q1_data, const PRECISION * const __restrict__ Q2_data, const PRECISION * const __restrict__ V_data, PRECISION u_i, PRECISION ut, PRECISION * const __restrict__ H,
	PRECISION (* const rightHalfCellExtrapolation)(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp),
	PRECISION (* const leftHalfCellExtrapolation)(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp))
{
	// neighbor fluid velocities
	PRECISION vmm = V_data[0];
	PRECISION vm  = V_data[1];
	PRECISION v   = u_i / ut;
	PRECISION vp  = V_data[2];
	PRECISION vpp = V_data[3];

	// left / right extrapolated speeds
	PRECISION vL = fabs(leftHalfCellExtrapolation(vmm, vm, v, vp, vpp));
	PRECISION vR = fabs(rightHalfCellExtrapolation(vmm, vm, v, vp, vpp));

	// local propagation speed
	PRECISION a = fmax(vL, vR);

	// compute left / right extrapolated values of F = vq, q and v
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

		qL = leftHalfCellExtrapolation(qmm, qm, q, qp, qpp);
		qR = rightHalfCellExtrapolation(qmm, qm, q, qp, qpp);

		FL = leftHalfCellExtrapolation(Fmm, Fm, F, Fp, Fpp);
		FR = rightHalfCellExtrapolation(Fmm, Fm, F, Fp, Fpp);

		// H from Eq.(61)
		H[n] = (FR + FL  -  a * (qR - qL)) / 2.0;

		p += 2;
	}
}



void flux(const PRECISION * const __restrict__ Q_data, PRECISION * const __restrict__ H,
	PRECISION (* const rightHalfCellExtrapolation)(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp),
	PRECISION (* const leftHalfCellExtrapolation)(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp),
	PRECISION (* const spectralRadius)(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un),
	PRECISION (* const fluxFunction)(PRECISION q, PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un),
	PRECISION t, PRECISION ePrev, PRECISION uT, int i, int j, int k)
{
	PRECISION t2 = t * t;

	PRECISION qmm, qm, q, qp, qpp;

	PRECISION qL[NUMBER_CONSERVED_VARIABLES];
	PRECISION qR[NUMBER_CONSERVED_VARIABLES];

	int p = 0;

	for(short n = 0; n < NUMBER_CONSERVED_VARIABLES; n++)
	{
		qmm = Q_data[p];
		qm 	= Q_data[p + 1];
		q 	= Q_data[p + 2];
		qp 	= Q_data[p + 3];
		qpp = Q_data[p + 4];

		qL[n] = leftHalfCellExtrapolation(qmm, qm, q, qp, qpp);
		qR[n] = rightHalfCellExtrapolation(qmm, qm, q, qp, qpp);

		p += 5;
	}


	// left and right extrapolated values of the primary variables (e, p, u^mu)
	// specifically need the left/right extrapolated values of u^mu to compute the local propagation speed
	PRECISION eL, pL, utL, uxL, uyL, unL;
	PRECISION eR, pR, utR, uxR, uyR, unR;

	// this function needs to be replaced
	int statusL = get_inferred_variables(t, qL, ePrev, uT, &eL, &pL, &utL, &uxL, &uyL, &unL, i, j, k, 0, 0, 0, 0);
	int statusR = get_inferred_variables(t, qR, ePrev, uT, &eR, &pR, &utR, &uxR, &uyR, &unR, i, j, k, 0, 0, 0, 0);


	//utL = sqrt(1.0  +  uxL * uxL  +  uyL * uyL  +  t2 * unL * unL);		// renormalize fluid velocity
	//utR = sqrt(1.0  +  uxR * uxR  +  uyR * uyR  +  t2 * unR * unR);


	PRECISION a = localPropagationSpeed(utR, uxR, uyR, unR, utL, uxL, uyL, unL, spectralRadius);
	PRECISION qR_n, qL_n, F_R, F_L;

	for(short n = 0; n < NUMBER_CONSERVED_VARIABLES; n++)
	{
		qL_n = qL[n];
		qR_n = qR[n];

		F_L = fluxFunction(qL_n, utL, uxL, uyL, unL);
		F_R = fluxFunction(qR_n, utR, uxR, uyR, unR);

		// H from Eq.(61)
		H[n] = (F_R + F_L  -  a * (qR_n - qL_n)) / 2.0;
	}
}




