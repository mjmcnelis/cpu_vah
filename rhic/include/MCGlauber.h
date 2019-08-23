#ifndef MCGLAUBER_H_
#define MCGLAUBER_H_

#include "Parameters.h"

void MC_Glauber_energy_density_transverse_profile(double * const __restrict__ energyDensityTransverse, int nx, int ny, double dx, double dy, initial_condition_parameters initial);

#endif