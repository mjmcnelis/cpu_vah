#ifndef HYDRODYNAMICS_H_
#define HYDRODYNAMICS_H_

#include "Parameters.h"

const double hbarc = 0.197326938;

void run_hydro(void * latticeParams, void * initCondParams, void * hydroParams, hydro_parameters hydro);

#endif