/*
 * SemiDiscreteKurganovTadmorScheme.h
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#ifndef FLUXTERMS_H_
#define FLUXTERMS_H_

#include "DynamicalVariables.h"

#define THETA 1.8	// flux limiter parameter
					// set near 1 for fluctuating initial conditions (e.g. 1.1)
					// set near 2 for smooth initial conditions (e.g. 1.8)

// qmm = q_{i-2}	|	qm = q_{i-1}	|	q = q_i		|	qp = q_{i+1}	|	qpp = q_{i+2}
//
PRECISION right_half_cell_extrapolation_forward(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp);
PRECISION left_half_cell_extrapolation_forward(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp);
PRECISION right_half_cell_extrapolation_backward(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp);
PRECISION left_half_cell_extrapolation_backward(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp);

// semi-discrete KT algorithm: computes the H-flux term Eq.(61) for a given spatial grid point
// updated version on 6/13
void flux_terms(PRECISION * const __restrict__ H, const PRECISION * const __restrict__ Q_data, const PRECISION * const __restrict__ Q1_data, const PRECISION * const __restrict__ Q2_data, const PRECISION * const __restrict__ V_data, PRECISION u_i, PRECISION ut,
	PRECISION (* const right_half_cell_extrapolation)(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp),
	PRECISION (* const left_half_cell_extrapolation)(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp));

#endif


