
#ifndef FLUXTERMS_H_
#define FLUXTERMS_H_

#include "Precision.h"

#define THETA 1.8	// flux limiter parameter
					// set near 1 for fluctuating initial conditions (e.g. 1.1)
					// set near 2 for smooth initial conditions (e.g. 1.8)

// qmm = q_{i-2}	|	qm = q_{i-1}	|	q = q_i		|	qp = q_{i+1}	|	qpp = q_{i+2}
//
precision right_half_cell_extrapolation_forward(precision qmm, precision qm, precision q, precision qp, precision qpp);
precision left_half_cell_extrapolation_forward(precision qmm, precision qm, precision q, precision qp, precision qpp);
precision right_half_cell_extrapolation_backward(precision qmm, precision qm, precision q, precision qp, precision qpp);
precision left_half_cell_extrapolation_backward(precision qmm, precision qm, precision q, precision qp, precision qpp);

// semi-discrete KT algorithm: computes the H-flux term Eq.(61) for a given spatial grid point
// updated version on 6/13
void flux_terms(precision * const __restrict__ H, const precision * const __restrict__ q_data, const precision * const __restrict__ q1_data, const precision * const __restrict__ q2_data, const precision * const __restrict__ ui_data, const precision * const __restrict__ ut_data, precision ui, precision ut,
	precision (* const right_half_cell_extrapolation)(precision qmm, precision qm, precision q, precision qp, precision qpp),
	precision (* const left_half_cell_extrapolation)(precision qmm, precision qm, precision q, precision qp, precision qpp));

#endif


