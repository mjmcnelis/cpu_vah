/*
 * DynamicalVariables.h
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#ifndef DYNAMICALVARIABLES_H_
#define DYNAMICALVARIABLES_H_

#define NUMBER_CONSERVATION_LAWS 4

// keep these undefined for now
//#define PIMUNU
//#define W_TZ_MU		// there's a segfault when turn on  6/13 (because PIMUNU needs to be turned on too)
// force both of them to turn on or none 

/*********************************************************/

#ifndef PIMUNU
#define NUMBER_PROPAGATED_PIMUNU_COMPONENTS 0
#else
#define NUMBER_PROPAGATED_PIMUNU_COMPONENTS 10
#endif

#ifndef W_TZ_MU
#define NUMBER_PROPAGATED_WTZMU_COMPONENTS 0
#else
#define NUMBER_PROPAGATED_WTZMU_COMPONENTS 4
#endif

#define NUMBER_DISSIPATIVE_CURRENTS (NUMBER_PROPAGATED_WTZMU_COMPONENTS + NUMBER_PROPAGATED_PIMUNU_COMPONENTS)

// what do I do with this? find out what this marco does
// Because of PL matching it's not exactly ideal hydro
#if NUMBER_DISSIPATIVE_CURRENTS == 0
#define IDEAL
#endif

// the +1 is for PL
#define NUMBER_CONSERVED_VARIABLES (NUMBER_CONSERVATION_LAWS + 1 + NUMBER_DISSIPATIVE_CURRENTS)
/*********************************************************/

#define PRECISION double

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


// this was an old debugging
// can replace it with something else
typedef struct
{
	PRECISION *knudsenNumberTaupiT;
	PRECISION *knudsenNumberTaupiL;
	PRECISION *knudsenNumberTaupi;
	PRECISION *knudsenNumberTauPi;
	PRECISION *Rpi;
	PRECISION *RPi;
	PRECISION *Rw;
	PRECISION *Rpi2;
	PRECISION *RPi2;
	PRECISION *Rw2;
	PRECISION *fTSolution;
	PRECISION *regulations;
	PRECISION *regMag;
	PRECISION *regTr;
	PRECISION *regU0;
	PRECISION *regU1;
	PRECISION *regU2;
	PRECISION *regU3;
	PRECISION *regZ0;
	PRECISION *regZ1;
	PRECISION *regZ2;
	PRECISION *regZ3;
	PRECISION *stt;
	PRECISION *sxx;
	PRECISION *syy;
	PRECISION *snn;
	PRECISION *taupi;
	PRECISION *dxux;
	PRECISION *dyuy;
	PRECISION *theta;
} VALIDITY_DOMAIN;

// q, u = current variables
// up = previous fluid velocity
// Q = updated variables
// qS, uS = intermediate variables
// extern means the variables are declared but defined elsewhere to allow other source files to use it

extern CONSERVED_VARIABLES *q, *Q, *qS;
extern FLUID_VELOCITY *u, *up, *uS;
extern PRECISION *e, *p;


// for debugging only
extern VALIDITY_DOMAIN *validityDomain;
extern double *fTSol_X1,*fTSol_Y1,*fTSol_1,*fTSol_X2,*fTSol_Y2,*fTSol_2;

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


