
#ifndef GHOSTCELLS_H_
#define GHOSTCELLS_H_

#include "Precision.h"
#include "DynamicalVariables.h"

// ghost cell boundary conditions
void set_ghost_cells(CONSERVED_VARIABLES * const __restrict__ q, precision * const __restrict__ e, FLUID_VELOCITY * const __restrict__ u, int nx, int ny, int nz);

#endif


