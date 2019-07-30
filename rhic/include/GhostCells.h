
#ifndef GHOSTCELLS_H_
#define GHOSTCELLS_H_

#include "Precision.h"
#include "DynamicalVariables.h"

// ghost cell boundary conditions
void set_ghost_cells(hydro_variables * const __restrict__ q, precision * const __restrict__ e, fluid_velocity * const __restrict__ u, int nx, int ny, int nz);

#endif


