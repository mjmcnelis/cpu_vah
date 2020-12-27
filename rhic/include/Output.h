
#ifndef OUTPUT_H_
#define OUTPUT_H_

#include "Parameters.h"

int central_index(lattice_parameters lattice);

void output_freezeout_slice_x(double t, lattice_parameters lattice, hydro_parameters hydro);
void output_freezeout_slice_z(double t, lattice_parameters lattice, hydro_parameters hydro);

void output_hydro_simulation(double t, double dt_prev, lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro);

#endif
