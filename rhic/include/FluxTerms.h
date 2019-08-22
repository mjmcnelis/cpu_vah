
#ifndef FLUXTERMS_H_
#define FLUXTERMS_H_

#include "Precision.h"

#define THETA 1.8	// flux limiter parameter
					// set near 1 for fluctuating initial conditions (e.g. 1.1)
					// set near 2 for smooth initial conditions (e.g. 1.8)

precision approximate_derivative(precision qm, precision q, precision qp);

precision compute_max_local_propagation_speed(const precision * const __restrict__ v_data, precision v);

void flux_terms(precision * const __restrict__ Hp, precision * const __restrict__ Hm, const precision * const __restrict__ q_data, const precision * const __restrict__ q1_data, const precision * const __restrict__ q2_data, const precision * const __restrict__ v_data, precision v);


#endif


