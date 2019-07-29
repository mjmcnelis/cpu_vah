
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"
#include "../include/Parameters.h"
using namespace std;

conserved_variables *q, *Q, *qS;	
fluid_velocity *u, *up, *uS;	
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
	uS = (fluid_velocity *)calloc(len, sizeof(fluid_velocity));

	q  = (conserved_variables *)calloc(len, sizeof(conserved_variables));
	Q  = (conserved_variables *)calloc(len, sizeof(conserved_variables));
	qS = (conserved_variables *)calloc(len, sizeof(conserved_variables));
}



void swap_conserved_variables(conserved_variables **arr1, conserved_variables **arr2)
{
	// don't understand the double pointer thing
	conserved_variables *tmp = *arr1;
	*arr1 = *arr2;
	*arr2 = tmp;
}


void set_current_conserved_variables()
{
	swap_conserved_variables(&q, &Q);
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
	free(uS);
	free(q);
	free(qS);
	free(Q);
}



