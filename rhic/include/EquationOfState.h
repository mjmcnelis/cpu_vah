
#ifndef EQUATIONOFSTATE_H_
#define EQUATIONOFSTATE_H_

#include "Precision.h"

#define CONFORMAL_EOS

precision energy_density_cutoff(precision e_min, precision e);

precision equilibriumPressure(precision e);

precision speedOfSoundSquared(precision e);

precision effectiveTemperature(precision e, precision conformal_eos_prefactor);

precision equilibriumEnergyDensity(precision T, precision conformal_eos_prefactor);


#endif
