
#ifndef GHOSTCELLS_H_
#define GHOSTCELLS_H_

#include "DynamicalVariables.h"

// ghost cell boundary conditions
void set_ghost_cells(CONSERVED_VARIABLES * const __restrict__ q, PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, FLUID_VELOCITY * const __restrict__ u, int nx, int ny, int nz, int ncx, int ncy);

#endif


