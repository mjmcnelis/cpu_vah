
#ifndef DYNAMICALVARIABLES_H_
#define DYNAMICALVARIABLES_H_

#include "EquationOfState.h"
#include "Precision.h"
#include "Parameters.h"

// should have a macros.h (only for tunable macro parameters?)

#define ANISO_HYDRO				// comment to run 2nd order viscous hydrodynamics
								// (for cross-checking vahydro in limit of small viscosity)

//#define PIMUNU 					// name shared but different shear stresses (maybe isolate them)

//#define VORTICITY				// include the vorticity source terms in the relaxation equations


#ifdef ANISO_HYDRO
	#ifdef CONFORMAL_EOS
		#define PT_MATCHING 0
	#else
		#define PT_MATCHING 1
	#endif
	//#define WTZMU 			// uncomment to turn on WTz (comment for 2+1d simulations b/c it's zero)
#else
	#ifndef CONFORMAL_EOS
		#define PI
	#endif
#endif


#ifdef PIMUNU
	#define PIMUNU_SCALAR 1
	#define PIMUNU_COMPONENTS 10
#else
	#define PIMUNU_SCALAR 0
	#define PIMUNU_COMPONENTS 0
#endif


#ifdef WTZMU
	#define WTZMU_SCALAR 1
	#define WTZMU_COMPONENTS 4
#else
	#define WTZMU_SCALAR 0
	#define WTZMU_COMPONENTS 0
#endif


#ifdef PI
	#define PI_COMPONENTS 1
#else
	#define PI_COMPONENTS 0
#endif


#define NUMBER_OF_RESIDUAL_CURRENTS (PIMUNU_COMPONENTS + WTZMU_COMPONENTS)
#define NUMBER_OF_VISCOUS_CURRENTS (PIMUNU_COMPONENTS + PI_COMPONENTS)

#ifdef ANISO_HYDRO
	#define NUMBER_CONSERVED_VARIABLES (5 + PT_MATCHING + NUMBER_OF_RESIDUAL_CURRENTS)
	#define NUMBER_HYDRO_SCALARS (3 + PT_MATCHING + PIMUNU_SCALAR + WTZMU_SCALAR)
#else
	#define NUMBER_CONSERVED_VARIABLES (4 + NUMBER_OF_VISCOUS_CURRENTS)
#endif


typedef struct
{
	precision ttt;
	precision ttx;
	precision tty;
	precision ttn;

#ifdef ANISO_HYDRO
	precision pl;

#if (PT_MATCHING == 1)
	precision pt;
#endif
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

#ifdef PI
	precision Pi;
#endif

} hydro_variables;

typedef struct
{
	precision ux;
	precision uy;
	precision un;
} fluid_velocity;


// q, u = current hydro, fluid velocity variables
// up = previous fluid velocity
// Q = updated hydro variables
// qI, uI = intermediate variables
// extern means the variables are declared but defined elsewhere to allow other source files to use it

extern hydro_variables *q, *Q, *qI;
extern fluid_velocity *u, *up, *uI;
extern precision *e;

// swap q <-> Q
void set_current_hydro_variables();

// swap u <-> up
void swap_fluid_velocity(fluid_velocity ** arr1, fluid_velocity ** arr2);

// memory
void allocate_memory(lattice_parameters lattice);
void free_memory();

#endif









