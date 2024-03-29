
#ifndef INFERREDVARIABLES_H_
#define INFERREDVARIABLES_H_

#include "Precision.h"
#include "DynamicalVariables.h"
#include "Parameters.h"

void set_inferred_variables_aniso_hydro(precision t, const hydro_variables * const __restrict__ q, precision * const __restrict__ e, fluid_velocity * const __restrict__ u, lattice_parameters lattice, hydro_parameters hydro);

void set_inferred_variables_viscous_hydro(precision t, const hydro_variables * const __restrict__ q, precision * const __restrict__ e, fluid_velocity * const __restrict__ u, lattice_parameters lattice, hydro_parameters hydro);

#endif



