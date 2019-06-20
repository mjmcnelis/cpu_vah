/*
 * EnergyMomentumTensor.h
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#ifndef ENERGYMOMENTUMTENSOR_H_
#define ENERGYMOMENTUMTENSOR_H_

#include "DynamicalVariables.h"

// this will need to be removed
PRECISION getTransverseFluidVelocityMagnitude(const FLUID_VELOCITY * const __restrict__ u, int s);



// get e, p, u^mu
// my version
void get_inferred_variables_new(PRECISION t, const PRECISION * const __restrict__ q, PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, PRECISION * const __restrict__ ut, PRECISION * const __restrict__ ux, PRECISION * const __restrict__ uy, PRECISION * const __restrict__ un);

void set_inferred_variables_new(const CONSERVED_VARIABLES * const __restrict__ q, PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, FLUID_VELOCITY * const __restrict__ u, PRECISION t, int nx, int ny, int nz, int ncx, int ncy);


// get e, p, u^mu
// for my case I don't need ePrev, uT_0, indices, fullTimeStepInversion
// here q[n] is an array that holds the conserved variables
int get_inferred_variables(PRECISION t, const PRECISION * const __restrict__ q, PRECISION ePrev, PRECISION uT_0, PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, PRECISION * const __restrict__ ut, PRECISION * const __restrict__ ux, PRECISION * const __restrict__ uy, PRECISION * const __restrict__ un, int i, int j, int k, double xi, double yj, double zk, int fullTimeStepInversion);

// here we're setting u, up (and q,e,p?)
// this calls getInferredVariables right?
void set_inferred_variables(const CONSERVED_VARIABLES * const __restrict__ q, PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, FLUID_VELOCITY * const __restrict__ up, FLUID_VELOCITY * const __restrict__ u, PRECISION t, int nx, int ny, int nz, int ncx, int ncy, PRECISION *fTSolution);

#endif



