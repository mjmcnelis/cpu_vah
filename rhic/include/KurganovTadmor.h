/*
 * FullyDiscreteKurganovTadmorScheme.h
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#ifndef KURGANOVTADMOR_H_
#define KURGANOVTADMOR_H_

#include "DynamicalVariables.h"

void rungeKutta2(PRECISION t, PRECISION dt, CONSERVED_VARIABLES * __restrict__ q, CONSERVED_VARIABLES * __restrict__ Q, int nx, int ny, int nz, int ncx, int ncy, PRECISION dx, PRECISION dy, PRECISION dz, PRECISION etabar);

#endif
