/*
 * SourceTerms.h
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#ifndef SOURCETERMS_H_
#define SOURCETERMS_H_

#include "DynamicalVariables.h"


void source_terms(PRECISION * const __restrict__ S, const PRECISION * const __restrict__ Q, PRECISION e_s, PRECISION t, const PRECISION * const __restrict__ qi1, const PRECISION * const __restrict__ qj1, const PRECISION * const __restrict__ qk1, const PRECISION * const __restrict__ e1, const PRECISION * const __restrict__ ui1, const PRECISION * const __restrict__ uj1, const PRECISION * const __restrict__ uk1, PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un, PRECISION ut_p, PRECISION ux_p, PRECISION uy_p, PRECISION un_p, PRECISION dt, PRECISION dx, PRECISION dy, PRECISION dn, PRECISION etabar);


#endif
