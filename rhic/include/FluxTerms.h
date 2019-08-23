
#ifndef FLUXTERMS_H_
#define FLUXTERMS_H_

#include "Precision.h"

precision approximate_derivative(precision qm, precision q, precision qp);

precision compute_max_local_propagation_speed(const precision * const __restrict__ v_data, precision v, precision Theta);

void flux_terms(precision * const __restrict__ Hp, precision * const __restrict__ Hm, const precision * const __restrict__ q_data, const precision * const __restrict__ q1_data, const precision * const __restrict__ q2_data, const precision * const __restrict__ v_data, precision v, precision Theta);


#endif


