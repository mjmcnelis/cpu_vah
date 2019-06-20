/*
 * EnergyMomentumTensor.h
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#ifndef INFERREDVARIABLES_H_
#define INFERREDVARIABLES_H_

#include "DynamicalVariables.h"

// compute e, p, u^mu
void set_inferred_variables(const CONSERVED_VARIABLES * const __restrict__ q, PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, FLUID_VELOCITY * const __restrict__ u, PRECISION t, int nx, int ny, int nz, int ncx, int ncy);

#endif



