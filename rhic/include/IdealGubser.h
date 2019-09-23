#ifndef IDEALGUBSER_H_
#define IDEALGUBSER_H_

#include "Parameters.h"

void set_ideal_gubser_initial_conditions(lattice_parameters lattice, precision dt, initial_condition_parameters initial, hydro_parameters hydro);

void run_analytic_ideal_gubser(lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro);

#endif