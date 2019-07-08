
#ifndef DYNAMICALVARIABLES_H_
#define DYNAMICALVARIABLES_H_

#include "EquationofState.h"
#include "Precision.h"

#define NUMBER_CONSERVATION_LAWS 4

// keep these undefined for now
//#define PIMUNU
//#define WTZMU		// there's a segfault when turn on  6/13 (because PIMUNU needs to be turned on too)
// force both of them to turn on or none

#ifndef PIMUNU
	#define PIMUNU_COMPONENTS 0
#else
	#define PIMUNU_COMPONENTS 10
#endif


#ifndef WTZMU
	#define WTZMU_COMPONENTS 0
#else
	#define WTZMU_COMPONENTS 4
#endif


#define NUMBER_RESIDUAL_CURRENTS (PIMUNU_COMPONENTS + WTZMU_COMPONENTS)

#define PL_MATCHING 1

#define PT_MATCHING 1
// #ifndef CONFORMAL_EOS
// 	#define PT_MATCHING 1
// #else
// 	#define PT_MATCHING 0
// #endif


#define NUMBER_CONSERVED_VARIABLES (NUMBER_CONSERVATION_LAWS + PL_MATCHING + PT_MATCHING + NUMBER_RESIDUAL_CURRENTS)

typedef struct
{
	precision * ttt;
	precision * ttx;
	precision * tty;
	precision * ttn;
	precision * pl;

#if (PT_MATCHING == 1)
	precision * pt;
#endif

#ifdef PIMUNU
	precision * pitt;
	precision * pitx;
	precision * pity;
	precision * pitn;
	precision * pixx;
	precision * pixy;
	precision * pixn;
	precision * piyy;
	precision * piyn;
	precision * pinn;
#endif
#ifdef W_TZ_MU
	precision * WtTz;
	precision * WxTz;
	precision * WyTz;
	precision * WnTz;
#endif
} CONSERVED_VARIABLES;

typedef struct
{
	precision * ut;
	precision * ux;
	precision * uy;
	precision * un;
} FLUID_VELOCITY;

// q, u = current variables
// up = previous fluid velocity
// Q = updated variables
// qS, uS = intermediate variables
// extern means the variables are declared but defined elsewhere to allow other source files to use it

extern CONSERVED_VARIABLES *q, *Q, *qS;	// should call it uI
extern FLUID_VELOCITY *u, *up, *uS;
extern precision *e;

// swap q <-> Q
void set_current_conserved_variables();

// swap u <-> up
void swap_fluid_velocity(FLUID_VELOCITY ** arr1, FLUID_VELOCITY ** arr2);

void allocate_memory(int len);	// memory
void free_memory();

#endif


