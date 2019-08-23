#ifndef PRINT_H_
#define PRINT_H_

#include "Parameters.h"

void print_hydro_mode();
void print_hydro_center(double t, double e, int s);

void print_parameters(lattice_parameters lattice, double t0, double T_switch, double etabar);

#endif