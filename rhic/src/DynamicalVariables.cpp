/*
 * DynamicalVariables.cpp
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "../include/DynamicalVariables.h"
#include "../include/Parameters.h"

using namespace std;

CONSERVED_VARIABLES *q, *Q, *qS;	// the extern variables are defined here
FLUID_VELOCITY *u, *up, *uS;		// so what is this purpose then?
PRECISION *e;


inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}

void allocate_memory(int len)
{
	size_t bytes = sizeof(PRECISION);

	// primary variables
	e = (PRECISION *)calloc(len, bytes);

	// fluid velocity at current time step
	u = (FLUID_VELOCITY *)calloc(1, sizeof(FLUID_VELOCITY));
	u->ut = (PRECISION *)calloc(len, bytes);
	u->ux = (PRECISION *)calloc(len, bytes);
	u->uy = (PRECISION *)calloc(len, bytes);
	u->un = (PRECISION *)calloc(len, bytes);

	// fluid velocity at previous time step
	up = (FLUID_VELOCITY *)calloc(1, sizeof(FLUID_VELOCITY));
	up->ut = (PRECISION *)calloc(len, bytes);
	up->ux = (PRECISION *)calloc(len, bytes);
	up->uy = (PRECISION *)calloc(len, bytes);
	up->un = (PRECISION *)calloc(len, bytes);

	// fluid velocity at intermediate time step
	uS = (FLUID_VELOCITY *)calloc(1, sizeof(FLUID_VELOCITY));
	uS->ut = (PRECISION *)calloc(len, bytes);
	uS->ux = (PRECISION *)calloc(len, bytes);
	uS->uy = (PRECISION *)calloc(len, bytes);
	uS->un = (PRECISION *)calloc(len, bytes);


	// conserved variables at current time step
	q = (CONSERVED_VARIABLES *)calloc(1, sizeof(CONSERVED_VARIABLES));
	q->ttt = (PRECISION *)calloc(len, bytes);
	q->ttx = (PRECISION *)calloc(len, bytes);
	q->tty = (PRECISION *)calloc(len, bytes);
	q->ttn = (PRECISION *)calloc(len, bytes);
	q->pl  = (PRECISION *)calloc(len, bytes);
	q->pt  = (PRECISION *)calloc(len, bytes);

#ifdef PIMUNU
	q->pitt = (PRECISION *)calloc(len, bytes);
	q->pitx = (PRECISION *)calloc(len, bytes);
	q->pity = (PRECISION *)calloc(len, bytes);
	q->pitn = (PRECISION *)calloc(len, bytes);
	q->pixx = (PRECISION *)calloc(len, bytes);
	q->pixy = (PRECISION *)calloc(len, bytes);
	q->pixn = (PRECISION *)calloc(len, bytes);
	q->piyy = (PRECISION *)calloc(len, bytes);
	q->piyn = (PRECISION *)calloc(len, bytes);
	q->pinn = (PRECISION *)calloc(len, bytes);
#endif
#ifdef WTZMU
	q->WtTz = (PRECISION *)calloc(len, bytes);
	q->WxTz = (PRECISION *)calloc(len, bytes);
	q->WyTz = (PRECISION *)calloc(len, bytes);
	q->WnTz = (PRECISION *)calloc(len, bytes);
#endif

	// conversed variables at intermediate time step
	qS = (CONSERVED_VARIABLES *)calloc(1, sizeof(CONSERVED_VARIABLES));
	qS->ttt = (PRECISION *)calloc(len, bytes);
	qS->ttx = (PRECISION *)calloc(len, bytes);
	qS->tty = (PRECISION *)calloc(len, bytes);
	qS->ttn = (PRECISION *)calloc(len, bytes);
	qS->pl  = (PRECISION *)calloc(len, bytes);
	qS->pt  = (PRECISION *)calloc(len, bytes);

#ifdef PIMUNU
	qS->pitt = (PRECISION *)calloc(len, bytes);
	qS->pitx = (PRECISION *)calloc(len, bytes);
	qS->pity = (PRECISION *)calloc(len, bytes);
	qS->pitn = (PRECISION *)calloc(len, bytes);
	qS->pixx = (PRECISION *)calloc(len, bytes);
	qS->pixy = (PRECISION *)calloc(len, bytes);
	qS->pixn = (PRECISION *)calloc(len, bytes);
	qS->piyy = (PRECISION *)calloc(len, bytes);
	qS->piyn = (PRECISION *)calloc(len, bytes);
	qS->pinn = (PRECISION *)calloc(len, bytes);
#endif
#ifdef WTZMU
	qS->WtTz = (PRECISION *)calloc(len, bytes);
	qS->WxTz = (PRECISION *)calloc(len, bytes);
	qS->WyTz = (PRECISION *)calloc(len, bytes);
	qS->WnTz = (PRECISION *)calloc(len, bytes);
#endif

	// conserved variables at next time step
	Q = (CONSERVED_VARIABLES *)calloc(1, sizeof(CONSERVED_VARIABLES));
	Q->ttt = (PRECISION *)calloc(len, bytes);
	Q->ttx = (PRECISION *)calloc(len, bytes);
	Q->tty = (PRECISION *)calloc(len, bytes);
	Q->ttn = (PRECISION *)calloc(len, bytes);
	Q->pl  = (PRECISION *)calloc(len, bytes);
	Q->pt  = (PRECISION *)calloc(len, bytes);

#ifdef PIMUNU
	Q->pitt = (PRECISION *)calloc(len, bytes);
	Q->pitx = (PRECISION *)calloc(len, bytes);
	Q->pity = (PRECISION *)calloc(len, bytes);
	Q->pitn = (PRECISION *)calloc(len, bytes);
	Q->pixx = (PRECISION *)calloc(len, bytes);
	Q->pixy = (PRECISION *)calloc(len, bytes);
	Q->pixn = (PRECISION *)calloc(len, bytes);
	Q->piyy = (PRECISION *)calloc(len, bytes);
	Q->piyn = (PRECISION *)calloc(len, bytes);
	Q->pinn = (PRECISION *)calloc(len, bytes);
#endif
#ifdef WTZMU
	Q->WtTz = (PRECISION *)calloc(len, bytes);
	Q->WxTz = (PRECISION *)calloc(len, bytes);
	Q->WyTz = (PRECISION *)calloc(len, bytes);
	Q->WnTz = (PRECISION *)calloc(len, bytes);
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
#ifdef W_TZ_MU
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



