
#ifndef ANISOBJORKEN
#define ANISOBJORKEN

#include "Parameters.h"

void run_semi_analytic_aniso_bjorken(lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro);

void set_aniso_bjorken_initial_condition(int nx, int ny, int nz, initial_condition_parameters initial, hydro_parameters hydro);

#endif