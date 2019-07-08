
#ifndef REGULATION_H_
#define REGULATION_H_

#include "Precision.h"
#include "DynamicalVariables.h"


void regulate_dissipative_currents(precision t, CONSERVED_VARIABLES * const __restrict__ q, precision * const __restrict__ e, const FLUID_VELOCITY * const __restrict__ u, int nx, int ny, int nz);


#endif


