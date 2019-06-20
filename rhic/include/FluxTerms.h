/*
 * SemiDiscreteKurganovTadmorScheme.h
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#ifndef FLUXTERMS_H_
#define FLUXTERMS_H_

#include "DynamicalVariables.h"


#define THETA 1.8	// flux limiter parameter (should be set near 1 for fluctuating initial conditions)

// what __restrict__ means:
// https://en.wikipedia.org/wiki/Restrict


// computes the finite spatial derivative with a minmod limiter
PRECISION approximateDerivative(PRECISION qm, PRECISION q, PRECISION qp);


// qmm = q_{i-2}	|	qm = q_{i-1}	|	q = q_i		|	qp = q_{i+1}	|	qpp = q_{i+2}
//
PRECISION rightHalfCellExtrapolationForward(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp);
PRECISION leftHalfCellExtrapolationForward(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp);
PRECISION rightHalfCellExtrapolationBackwards(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp);
PRECISION leftHalfCellExtrapolationBackwards(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp);



// spectral radius components u_i / ut for computing the local propagation speeds
//
PRECISION spectralRadiusX(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un);
PRECISION spectralRadiusY(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un);
PRECISION spectralRadiusZ(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un);



// flux vector components F_i = v_i * q of variable q
//
PRECISION Fx(PRECISION q, PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un);
PRECISION Fy(PRECISION q, PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un);
PRECISION Fz(PRECISION q, PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un);



// semi-discrete KT algorithm: computes the H-flux term Eq.(61) for a given spatial grid point
// updated version on 6/13
void flux_terms(const PRECISION * const __restrict__ Q_data, const PRECISION * const __restrict__ Q1_data, const PRECISION * const __restrict__ Q2_data, const PRECISION * const __restrict__ V_data, PRECISION u_i, PRECISION ut, PRECISION * const __restrict__ H,
	PRECISION (* const rightHalfCellExtrapolation)(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp),
	PRECISION (* const leftHalfCellExtrapolation)(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp));


void flux(const PRECISION * const __restrict__ Q_data, PRECISION * const __restrict__ H,
		PRECISION (* const rightHalfCellExtrapolation)(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp),
		PRECISION (* const leftHalfCellExtrapolation)(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp),
		PRECISION (* const spectralRadius)(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un),
		PRECISION (* const fluxFunction)(PRECISION q, PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un),
		PRECISION t, PRECISION ePrev, PRECISION uT, int i, int j, int k);

#endif


