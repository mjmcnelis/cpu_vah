
#ifndef DYNAMICALVARIABLES_H_
#define DYNAMICALVARIABLES_H_

#include "Macros.h"
#include "EquationOfState.h"
#include "Precision.h"
#include "Parameters.h"

#ifdef BOOST_INVARIANT
	#define CONSERVATION_LAWS 3
#else
	#define CONSERVATION_LAWS 4
#endif


#ifdef PIMUNU
	#ifdef BOOST_INVARIANT
		#ifdef ANISO_HYDRO
			#define PIMUNU_COMPONENTS 6
		#else
			#define PIMUNU_COMPONENTS 7
		#endif
	#else
		#define PIMUNU_COMPONENTS 10
	#endif
#else
	#define PIMUNU_COMPONENTS 0
#endif


#ifdef WTZMU
	#define WTZMU_COMPONENTS 4
#else
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
	#define NUMBER_CONSERVED_VARIABLES (CONSERVATION_LAWS + 1 + PT_MATCHING + NUMBER_OF_RESIDUAL_CURRENTS)
#else
	#define NUMBER_CONSERVED_VARIABLES (CONSERVATION_LAWS + NUMBER_OF_VISCOUS_CURRENTS)
#endif


// hydro variables
typedef struct
{
	precision ttt;
	precision ttx;
	precision tty;
#ifndef BOOST_INVARIANT
	precision ttn;
#endif

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
	precision pixx;
	precision pixy;
	precision piyy;

#ifndef BOOST_INVARIANT
	precision pitn;
	precision pixn;
	precision piyn;
	precision pinn;
#else
	#ifndef ANISO_HYDRO
		precision pinn;		// piperp^\eta\eta = 0 in 2+1d
	#endif
#endif
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


// fluid velocity
typedef struct
{
	precision ux;
	precision uy;
#ifndef BOOST_INVARIANT
	precision un;
#endif

} fluid_velocity;


// external variables
extern hydro_variables *q, *Q, *qI;
extern fluid_velocity *u, *up, *uI;
extern precision *e;


// swap variables
void swap_hydro_variables(hydro_variables ** hydro_1, hydro_variables ** hydro_2);
void swap_fluid_velocity(fluid_velocity ** velocity_1, fluid_velocity ** velocity_2);


// memory
void allocate_memory(lattice_parameters lattice);
void free_memory();

#endif









