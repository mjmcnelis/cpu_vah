
#ifndef INFERREDVARIABLES_H_
#define INFERREDVARIABLES_H_

#include "Precision.h"
#include "DynamicalVariables.h"
#include "Parameters.h"

typedef struct
{
	precision energy;
	precision ux;
	precision uy;
	precision un;

} inferred_variables;


inferred_variables solve_inferred_variables_aniso_hydro(hydro_variables q_aniso, precision t, precision t2, precision t4, hydro_parameters hydro);
inferred_variables solve_inferred_variables_viscous_hydro(hydro_variables q_viscous, precision t, precision e_prev, hydro_parameters hydro);

void set_inferred_variables_aniso_hydro(const hydro_variables * const __restrict__ q, precision * const __restrict__ e, fluid_velocity * const __restrict__ u, precision t, lattice_parameters lattice, hydro_parameters hydro);

void set_inferred_variables_viscous_hydro(const hydro_variables * const __restrict__ q, precision * const __restrict__ e, fluid_velocity * const __restrict__ u, precision t, lattice_parameters lattice, hydro_parameters hydro);

#endif



