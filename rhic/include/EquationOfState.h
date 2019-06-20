/*
 * EquationOfState.h
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#ifndef EQUATIONOFSTATE_H_
#define EQUATIONOFSTATE_H_

#include "DynamicalVariables.h"

#define CONFORMAL_EOS

// ideal gas of massless quarks and gluons
//#define EOS_FACTOR 15.6269 // Nc=3, Nf=3
#define EOS_FACTOR 13.8997 // Nc=3, Nf=2.5

PRECISION equilibriumPressure(PRECISION e);

PRECISION speedOfSoundSquared(PRECISION e);

PRECISION effectiveTemperature(PRECISION e);

PRECISION equilibriumEnergyDensity(PRECISION T);

PRECISION derivativeEnergyDensityWithRespectToTemperature(PRECISION T);

#endif /* EQUATIONOFSTATE_H_ */
