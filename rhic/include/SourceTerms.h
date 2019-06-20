/*
 * SourceTerms.h
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#ifndef SOURCETERMS_H_
#define SOURCETERMS_H_

#include "DynamicalVariables.h"




void source_terms(const PRECISION * const __restrict__ Q, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u,
PRECISION utp, PRECISION uxp, PRECISION uyp, PRECISION unp,
PRECISION t, const PRECISION * const __restrict__ evec, const PRECISION * const __restrict__ pvec,
int s, int d_ncx, int d_ncy, int d_ncz, PRECISION etabar, PRECISION dt, PRECISION dx, PRECISION dy, PRECISION dn,
int i, int j, int k, double x, double y, double z,
const CONSERVED_VARIABLES * const __restrict__ currentVars, const PRECISION * const __restrict__ qi1, const PRECISION * const __restrict__ qj1, const PRECISION * const __restrict__ qk1, const PRECISION * const __restrict__ e1, const PRECISION * const __restrict__ p1, PRECISION e_s, PRECISION p_s, const PRECISION * const __restrict__ ui1, const PRECISION * const __restrict__ uj1, const PRECISION * const __restrict__ uk1, PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un, PRECISION ut_p, PRECISION ux_p, PRECISION uy_p, PRECISION un_p);


void loadSourceTerms(const PRECISION * const __restrict__ Q, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u,
PRECISION utp, PRECISION uxp, PRECISION uyp, PRECISION unp,
PRECISION t, const PRECISION * const __restrict__ evec, const PRECISION * const __restrict__ pvec,
int s, int d_ncx, int d_ncy, int d_ncz, PRECISION d_etabar, PRECISION d_dt, PRECISION d_dx, PRECISION d_dy, PRECISION d_dz,
int i, int j, int k, double x, double y, double z,
const CONSERVED_VARIABLES * const __restrict__ currentVars);

// can I move this into source terms?
void loadSourceTermsX(const PRECISION * const __restrict__ I, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u, int s, PRECISION d_dx);
void loadSourceTermsY(const PRECISION * const __restrict__ J, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u, int s, PRECISION d_dy);
void loadSourceTermsZ(const PRECISION * const __restrict__ K, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u, int s, PRECISION t,PRECISION d_dz);

#endif
