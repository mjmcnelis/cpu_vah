
#ifndef DYNAMICALVARIABLES_H_
#define DYNAMICALVARIABLES_H_

#include "EquationOfState.h"
#include "Precision.h"

#define NUMBER_OF_CONSERVATION_LAWS 4

//#define PIMUNU
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
} hydro_variables;

typedef struct
{
	precision ux;
	precision uy;
	precision un;
} fluid_velocity;


// q, u = current variables
// up = previous fluid velocity
// Q = updated variables
// qI, uI = intermediate variables
// extern means the variables are declared but defined elsewhere to allow other source files to use it

extern hydro_variables *q, *Q, *qI;
extern fluid_velocity *u, *up, *uI;
extern precision *e;

// swap q <-> Q
void set_current_hydro_variables();

// swap u <-> up
void swap_fluid_velocity(fluid_velocity ** arr1, fluid_velocity ** arr2);

void allocate_memory(int len);	// memory
void free_memory();

#endif


