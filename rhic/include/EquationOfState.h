
#ifndef EQUATIONOFSTATE_H_
#define EQUATIONOFSTATE_H_

#include "Precision.h"

#define CONFORMAL_EOS

// ideal gas of massless quarks and gluons
#define EOS_FACTOR 15.6269 // Nc=3, Nf=3
//#define EOS_FACTOR 13.8997 // Nc=3, Nf=2.5

#define E_MIN 1.e-7		// minimum energy density

precision equilibriumPressure(precision e);

precision speedOfSoundSquared(precision e);

precision effectiveTemperature(precision e);

precision equilibriumEnergyDensity(precision T);


#endif
