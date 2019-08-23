
#ifndef FILEIO_H_
#define FILEIO_H_

#include "Parameters.h"

int central_index(lattice_parameters lattice);	

void output_dynamical_variables(double t, double dt, lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro);


#endif
