
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"
#include "../include/Parameters.h"

using namespace std;

CONSERVED_VARIABLES *q, *Q, *qS;	// the extern variables are defined here
FLUID_VELOCITY *u, *up, *uS;		// so what is this purpose then?
precision *e;


inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}

void allocate_memory(int len)
{
	size_t bytes = sizeof(precision);

	// primary variables
	e = (precision *)calloc(len, bytes);

	// fluid velocity at current time step
	u = (FLUID_VELOCITY *)calloc(1, sizeof(FLUID_VELOCITY));
	u->ut = (precision *)calloc(len, bytes);
	u->ux = (precision *)calloc(len, bytes);
	u->uy = (precision *)calloc(len, bytes);
	u->un = (precision *)calloc(len, bytes);

	// fluid velocity at previous time step
	up = (FLUID_VELOCITY *)calloc(1, sizeof(FLUID_VELOCITY));
	up->ut = (precision *)calloc(len, bytes);
	up->ux = (precision *)calloc(len, bytes);
	up->uy = (precision *)calloc(len, bytes);
	up->un = (precision *)calloc(len, bytes);

	// fluid velocity at intermediate time step
	uS = (FLUID_VELOCITY *)calloc(1, sizeof(FLUID_VELOCITY));
	uS->ut = (precision *)calloc(len, bytes);
	uS->ux = (precision *)calloc(len, bytes);
	uS->uy = (precision *)calloc(len, bytes);
	uS->un = (precision *)calloc(len, bytes);


	// conserved variables at current time step
	q = (CONSERVED_VARIABLES *)calloc(1, sizeof(CONSERVED_VARIABLES));
	q->ttt = (precision *)calloc(len, bytes);
	q->ttx = (precision *)calloc(len, bytes);
	q->tty = (precision *)calloc(len, bytes);
	q->ttn = (precision *)calloc(len, bytes);

	q->pl  = (precision *)calloc(len, bytes);

#if (PT_MATCHING == 1)
	q->pt  = (precision *)calloc(len, bytes);
#endif

#ifdef PIMUNU
	q->pitt = (precision *)calloc(len, bytes);
	q->pitx = (precision *)calloc(len, bytes);
	q->pity = (precision *)calloc(len, bytes);
	q->pitn = (precision *)calloc(len, bytes);
	q->pixx = (precision *)calloc(len, bytes);
	q->pixy = (precision *)calloc(len, bytes);
	q->pixn = (precision *)calloc(len, bytes);
	q->piyy = (precision *)calloc(len, bytes);
	q->piyn = (precision *)calloc(len, bytes);
	q->pinn = (precision *)calloc(len, bytes);
#endif
#ifdef WTZMU
	q->WtTz = (precision *)calloc(len, bytes);
	q->WxTz = (precision *)calloc(len, bytes);
	q->WyTz = (precision *)calloc(len, bytes);
	q->WnTz = (precision *)calloc(len, bytes);
#endif

	// conversed variables at intermediate time step
	qS = (CONSERVED_VARIABLES *)calloc(1, sizeof(CONSERVED_VARIABLES));
	qS->ttt = (precision *)calloc(len, bytes);
	qS->ttx = (precision *)calloc(len, bytes);
	qS->tty = (precision *)calloc(len, bytes);
	qS->ttn = (precision *)calloc(len, bytes);

	qS->pl  = (precision *)calloc(len, bytes);

#if (PT_MATCHING == 1)
	qS->pt  = (precision *)calloc(len, bytes);
#endif

#ifdef PIMUNU
	qS->pitt = (precision *)calloc(len, bytes);
	qS->pitx = (precision *)calloc(len, bytes);
	qS->pity = (precision *)calloc(len, bytes);
	qS->pitn = (precision *)calloc(len, bytes);
	qS->pixx = (precision *)calloc(len, bytes);
	qS->pixy = (precision *)calloc(len, bytes);
	qS->pixn = (precision *)calloc(len, bytes);
	qS->piyy = (precision *)calloc(len, bytes);
	qS->piyn = (precision *)calloc(len, bytes);
	qS->pinn = (precision *)calloc(len, bytes);
#endif
#ifdef WTZMU
	qS->WtTz = (precision *)calloc(len, bytes);
	qS->WxTz = (precision *)calloc(len, bytes);
	qS->WyTz = (precision *)calloc(len, bytes);
	qS->WnTz = (precision *)calloc(len, bytes);
#endif

	// conserved variables at next time step
	Q = (CONSERVED_VARIABLES *)calloc(1, sizeof(CONSERVED_VARIABLES));
	Q->ttt = (precision *)calloc(len, bytes);
	Q->ttx = (precision *)calloc(len, bytes);
	Q->tty = (precision *)calloc(len, bytes);
	Q->ttn = (precision *)calloc(len, bytes);

	Q->pl  = (precision *)calloc(len, bytes);

#if (PT_MATCHING == 1)
	Q->pt  = (precision *)calloc(len, bytes);
#endif

#ifdef PIMUNU
	Q->pitt = (precision *)calloc(len, bytes);
	Q->pitx = (precision *)calloc(len, bytes);
	Q->pity = (precision *)calloc(len, bytes);
	Q->pitn = (precision *)calloc(len, bytes);
	Q->pixx = (precision *)calloc(len, bytes);
	Q->pixy = (precision *)calloc(len, bytes);
	Q->pixn = (precision *)calloc(len, bytes);
	Q->piyy = (precision *)calloc(len, bytes);
	Q->piyn = (precision *)calloc(len, bytes);
	Q->pinn = (precision *)calloc(len, bytes);
#endif
#ifdef WTZMU
	Q->WtTz = (precision *)calloc(len, bytes);
	Q->WxTz = (precision *)calloc(len, bytes);
	Q->WyTz = (precision *)calloc(len, bytes);
	Q->WnTz = (precision *)calloc(len, bytes);
#endif
}


void swap_conserved_variables(CONSERVED_VARIABLES **arr1, CONSERVED_VARIABLES **arr2)
{
	// don't understand the double pointer thing
	CONSERVED_VARIABLES *tmp = *arr1;
	*arr1 = *arr2;
	*arr2 = tmp;
}


void set_current_conserved_variables()
{
	swap_conserved_variables(&q, &Q);
}


void swap_fluid_velocity(FLUID_VELOCITY **arr1, FLUID_VELOCITY **arr2)
{
	FLUID_VELOCITY *tmp = *arr1;
	*arr1 = *arr2;
	*arr2 = tmp;
}


void free_fluid_velocity(FLUID_VELOCITY * u)
{
	free(u->ut);
	free(u->ux);
	free(u->uy);
	free(u->un);
	free(u);
}


void free_conserved_variables(CONSERVED_VARIABLES * q)
{
	free(q->ttt);
	free(q->ttx);
	free(q->tty);
	free(q->ttn);
	free(q->pl);

#if (PT_MATCHING == 1)
	free(q->pt);
#endif

#ifdef PIMUNU
	free(q->pitt);
	free(q->pitx);
	free(q->pity);
	free(q->pitn);
	free(q->pixx);
	free(q->pixy);
	free(q->pixn);
	free(q->piyy);
	free(q->piyn);
	free(q->pinn);
#endif
#ifdef WTZMU
	free(q->WtTz);
	free(q->WxTz);
	free(q->WyTz);
	free(q->WnTz);
#endif
	free(q);
}


void free_memory()
{
	free(e);

	free_fluid_velocity(u);
	free_fluid_velocity(up);
	free_fluid_velocity(uS);

	free_conserved_variables(q);
	free_conserved_variables(qS);
	free_conserved_variables(Q);
}



