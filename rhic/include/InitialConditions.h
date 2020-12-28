
#ifndef INITIALCONDITIONS_H_
#define INITIALCONDITIONS_H_

#include "Parameters.h"
#include "Precision.h"
#include <vector>

void set_initial_conditions(precision t, lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro, std::vector<double> trento);

#endif
