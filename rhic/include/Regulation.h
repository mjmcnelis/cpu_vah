
#ifndef REGULATION_H_
#define REGULATION_H_

#include "DynamicalVariables.h"

void regulate_dissipative_currents(PRECISION t, CONSERVED_VARIABLES * const __restrict__ Q_current, PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, const FLUID_VELOCITY * const __restrict__ u, int nx, int ny, int nz, int ncx, int ncy);

#endif


