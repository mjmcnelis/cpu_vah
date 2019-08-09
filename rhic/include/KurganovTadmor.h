#ifndef KURGANOVTADMOR_H_
#define KURGANOVTADMOR_H_

#include "Precision.h"

void evolve_hydro_one_time_step(precision t, precision dt, precision dt_prev, int nx, int ny, int nz, precision dx, precision dy, precision dz, precision etabar_const);

#endif