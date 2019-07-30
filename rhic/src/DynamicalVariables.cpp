
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"
#include "../include/Parameters.h"
using namespace std;

hydro_variables *q, *Q, *qI;	
fluid_velocity *u, *up, *uI;	
precision *e;


inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}

void allocate_memory(int len)
{
	size_t bytes = sizeof(precision);

	e = (precision *)calloc(len, bytes);

	u  = (fluid_velocity *)calloc(len, sizeof(fluid_velocity));
	up = (fluid_velocity *)calloc(len, sizeof(fluid_velocity));
	uI = (fluid_velocity *)calloc(len, sizeof(fluid_velocity));

	q  = (hydro_variables *)calloc(len, sizeof(hydro_variables));
	Q  = (hydro_variables *)calloc(len, sizeof(hydro_variables));
	qI = (hydro_variables *)calloc(len, sizeof(hydro_variables));
}


void swap_hydro_variables(hydro_variables **arr1, hydro_variables **arr2)
{
	hydro_variables *tmp = *arr1;
	*arr1 = *arr2;
	*arr2 = tmp;
}


void set_current_hydro_variables()
{
	swap_hydro_variables(&q, &Q);
}


void swap_fluid_velocity(fluid_velocity **arr1, fluid_velocity **arr2)
{
	fluid_velocity *tmp = *arr1;
	*arr1 = *arr2;
	*arr2 = tmp;
}


void free_memory()
{
	free(e);
	free(u);
	free(up);
	free(uI);
	free(q);
	free(qI);
	free(Q);
}



