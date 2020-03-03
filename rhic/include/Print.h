#ifndef PRINT_H_
#define PRINT_H_

#include "Parameters.h"

void print_hydro_mode(hydro_parameters hydro);

void print_run_time(double duration, double steps, lattice_parameters lattice);

void print_hydro_center(int n, double t, lattice_parameters lattice, hydro_parameters hydro, long cells_above_Tswitch);

void print_parameters(lattice_parameters lattice, hydro_parameters hydro);

#endif