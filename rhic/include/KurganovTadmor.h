
#ifndef KURGANOVTADMOR_H_
#define KURGANOVTADMOR_H_

#include "Precision.h"
#include "Parameters.h"

void evolve_hydro_one_time_step(precision t, precision dt, precision dt_prev, lattice_parameters lattice, hydro_parameters hydro);

#endif