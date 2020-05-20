
#ifndef HYDRODYNAMICS_H_
#define HYDRODYNAMICS_H_

#include "Parameters.h"

const double hbarc = 0.197326938;

void run_hydro(lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro, int sample);

#endif