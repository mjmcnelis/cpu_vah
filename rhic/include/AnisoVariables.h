
#ifndef ANISOVARIABLES_H_
#define ANISOVARIABLES_H_

#include "Precision.h"

const int N_max = 200;	      			// max number of iterations   (probably need to adjust)
const precision tol_dX = 1.0e-6;    	// tolerance for dX
const precision tol_F = 1.0e-12;    	// tolerance for F 

typedef enum {newton, broyden} jacobian;

typedef struct 
{
	precision lambda;
	precision aT;
	precision aL;

} aniso_variables;

aniso_variables find_anisotropic_variables(precision e, precision pl, precision pt, precision B, precision mass, precision lambda_0, precision aT_0, precision aL_0);


#endif




