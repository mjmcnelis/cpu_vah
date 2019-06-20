/*
 * FullyDiscreteKurganovTadmorScheme.h
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#ifndef KURGANOVTADMOR_H_
#define KURGANOVTADMOR_H_

#include "DynamicalVariables.h"

void rungeKutta2(PRECISION t, PRECISION dt, int nx, int ny, int nz, PRECISION dx, PRECISION dy, PRECISION dz, PRECISION etabar);

#endif
