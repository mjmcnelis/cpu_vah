
#ifndef REGULATION_H_
#define REGULATION_H_

#include "Precision.h"
#include "DynamicalVariables.h"
#include "InferredVariables.h"
#include "Parameters.h"

hydro_variables regulate_viscous_currents_aniso(hydro_variables q_aniso, inferred_variables root, precision t, precision t2, precision t4, hydro_parameters hydro);
hydro_variables regulate_viscous_currents_viscous(hydro_variables q_viscous, inferred_variables root, precision t, precision t2, precision t4, hydro_parameters hydro);

void regulate_residual_currents(precision t, hydro_variables * const __restrict__ q, precision * const __restrict__ e, const fluid_velocity * const __restrict__ u, lattice_parameters lattice, hydro_parameters hydro);

void regulate_viscous_currents(precision t, hydro_variables * const __restrict__ q, precision * const __restrict__ e, const fluid_velocity * const __restrict__ u, lattice_parameters lattice, hydro_parameters hydro);

#endif


