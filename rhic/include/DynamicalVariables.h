/*
 * DynamicalVariables.h
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#ifndef DYNAMICALVARIABLES_H_
#define DYNAMICALVARIABLES_H_

#define PRECISION double

#define NUMBER_CONSERVATION_LAWS 4

// keep these undefined for now
//#define PIMUNU
//#define W_TZ_MU		// there's a segfault when turn on  6/13 (because PIMUNU needs to be turned on too)
// force both of them to turn on or none

#ifndef PIMUNU
#define PIMUNU_COMPONENTS 0
#else
#define PIMUNU_COMPONENTS 10
#endif

#ifndef W_TZ_MU
#define WTZMU_COMPONENTS 0
#else
#define WTZMU_COMPONENTS 4
#endif

#define NUMBER_RESIDUAL_CURRENTS (PIMUNU_COMPONENTS + WTZMU_COMPONENTS)

// the +1 is for pl
#define NUMBER_CONSERVED_VARIABLES (NUMBER_CONSERVATION_LAWS + 1 + NUMBER_RESIDUAL_CURRENTS)

typedef struct
{
	PRECISION * ttt;
	PRECISION * ttx;
	PRECISION * tty;
	PRECISION * ttn;
	PRECISION * pl;		// need to add pt
	//PRECISION *pt;
#ifdef PIMUNU
	PRECISION * pitt;
	PRECISION * pitx;
	PRECISION * pity;
	PRECISION * pitn;
	PRECISION * pixx;
	PRECISION * pixy;
	PRECISION * pixn;
	PRECISION * piyy;
	PRECISION * piyn;
	PRECISION * pinn;
#endif
#ifdef W_TZ_MU
	PRECISION * WtTz;
	PRECISION * WxTz;
	PRECISION * WyTz;
	PRECISION * WnTz;
#endif
} CONSERVED_VARIABLES;

typedef struct
{
	PRECISION * ut;
	PRECISION * ux;
	PRECISION * uy;
	PRECISION * un;
} FLUID_VELOCITY;

// q, u = current variables
// up = previous fluid velocity
// Q = updated variables
// qS, uS = intermediate variables
// extern means the variables are declared but defined elsewhere to allow other source files to use it

extern CONSERVED_VARIABLES *q, *Q, *qS;
extern FLUID_VELOCITY *u, *up, *uS;
extern PRECISION *e, *p;

// swap q <-> Q
void set_current_conserved_variables();

// swap u <-> up
void swap_fluid_velocity(FLUID_VELOCITY ** arr1, FLUID_VELOCITY ** arr2);

// ghost cell boundary conditions
void set_ghost_cells(CONSERVED_VARIABLES * const __restrict__ q, PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, FLUID_VELOCITY * const __restrict__ u, int nx, int ny, int nz, int ncx, int ncy);

// memory
void allocate_memory(int len);
void free_memory();

#endif


