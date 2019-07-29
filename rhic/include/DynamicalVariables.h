
#ifndef DYNAMICALVARIABLES_H_
#define DYNAMICALVARIABLES_H_

#include "EquationofState.h"
#include "Precision.h"

#define NUMBER_OF_CONSERVATION_LAWS 4

#define PIMUNU
//#define WTZMU

#ifdef PIMUNU
	#define PIMUNU_COMPONENTS 10
#else
	#define PIMUNU_COMPONENTS 0
#endif

#ifdef WTZMU
	#define WTZMU_COMPONENTS 4
#else
	#define WTZMU_COMPONENTS 0
#endif

#define NUMBER_OF_RESIDUAL_CURRENTS (PIMUNU_COMPONENTS + WTZMU_COMPONENTS)

#define PL_MATCHING 1

#ifdef CONFORMAL_EOS
	#define PT_MATCHING 0
#else
	#define PT_MATCHING 1
#endif

#define NUMBER_CONSERVED_VARIABLES (NUMBER_OF_CONSERVATION_LAWS + PL_MATCHING + PT_MATCHING + NUMBER_OF_RESIDUAL_CURRENTS)

typedef struct
{
	precision ttt;
	precision ttx;
	precision tty;
	precision ttn;
	precision pl;
#if (PT_MATCHING == 1)
	precision pt;
#endif
#ifdef PIMUNU
	precision pitt;
	precision pitx;
	precision pity;
	precision pitn;
	precision pixx;
	precision pixy;
	precision pixn;
	precision piyy;
	precision piyn;
	precision pinn;
#endif
#ifdef WTZMU
	precision WtTz;
	precision WxTz;
	precision WyTz;
	precision WnTz;
#endif
} conserved_variables;

typedef struct
{
	precision ux;
	precision uy;
	precision un;
} fluid_velocity;


// q, u = current variables
// up = previous fluid velocity
// Q = updated variables
// qS, uS = intermediate variables
// extern means the variables are declared but defined elsewhere to allow other source files to use it

extern conserved_variables *q, *Q, *qS;	// should call it qI, uI
extern fluid_velocity *u, *up, *uS;
extern precision *e;

void test_memory_time(int nt, int nx, int ny, int nz);

// swap q <-> Q
void set_current_conserved_variables();

// swap u <-> up
void swap_fluid_velocity(fluid_velocity ** arr1, fluid_velocity ** arr2);

void allocate_memory(int len);	// memory
void free_memory();

#endif


