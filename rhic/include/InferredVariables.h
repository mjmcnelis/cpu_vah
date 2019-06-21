
#ifndef INFERREDVARIABLES_H_
#define INFERREDVARIABLES_H_

#include "Precision.h"
#include "DynamicalVariables.h"


// compute e, p, u^mu
void set_inferred_variables(const CONSERVED_VARIABLES * const __restrict__ q, precision * const __restrict__ e, FLUID_VELOCITY * const __restrict__ u, precision t, int nx, int ny, int nz);


#endif



