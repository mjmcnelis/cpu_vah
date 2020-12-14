
#ifndef ANISOVARIABLES_H_
#define ANISOVARIABLES_H_

#include "Precision.h"
#include "DynamicalVariables.h"
#include "Parameters.h"

const int N_max = 100;	      			// max number of iterations   (probably need to adjust)
const precision tol_dX = 1.e-4;    		// tolerance for dX
const precision tol_F = 1.e-4;    		// tolerance for F

typedef struct
{
	precision lambda;
	precision aT;
	precision aL;
	int did_not_find_solution;
	int number_of_iterations;

} aniso_variables;

aniso_variables find_anisotropic_variables(precision e, precision pl, precision pt, precision B, precision mass, precision lambda_0, precision aT_0, precision aL_0);

void set_anisotropic_variables(const hydro_variables * const __restrict__ q, const precision * const __restrict__ e, precision * const __restrict__ lambda, precision * const __restrict__ aT, precision * const __restrict__ aL, lattice_parameters lattice, hydro_parameters hydro);

#endif




