#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include <gsl/gsl_linalg.h>
#include "../include/TransportAnisoNonconformal.h"
#include "../include/TransportAniso.h"
#include "../include/AnisoVariables.h"
#include "../include/Precision.h"
#include "../include/Macros.h"
using namespace std;


void compute_F(precision Ea, precision PTa, precision PLa, precision mass, precision * X, precision * F)
{
	precision lambda = X[0];
	precision aT = X[1];
	precision aL = X[2];

	precision aT2 = aT * aT;
	precision aL2 = aL * aL;
	precision aT2_minus_aL2 = aT2 - aL2;

	precision mbar = mass / lambda;
	precision mbar2 = mbar * mbar;

	precision prefactor = g * aT2 * aL * lambda * lambda * lambda * lambda / (4. * M_PI * M_PI);

	precision t_200, t_220, t_201;

	precision I_200 = 0;
	precision I_220 = 0;
	precision I_201 = 0;

	for(int i = 0; i < pbar_pts; i++)														// gauss integration
	{
		precision pbar = pbar_root_a2[i];													// pbar roots / weights for a = 2 (a = n + s)
		precision total_weight = pbar * pbar_weight_a2[i] * exp(pbar - sqrt(pbar * pbar + mbar2));

		precision w = sqrt(aL2  +  mbar2 / (pbar * pbar));
		precision z = aT2_minus_aL2 / (w * w);

		if(z > delta)																		// compute hypergeometric functions
		{
			precision sqrtz = sqrt(z);
			precision t = atan(sqrtz) / sqrtz;

			t_200 = 1.  +  (1. + z) * t;
			t_220 = (-1.  +  (1. + z) * t) / z;
			t_201 = (1.  +  (z - 1.) * t) / z;
		}
		else if(z < -delta && z > -1.)
		{
			precision sqrtmz = sqrt(-z);
			precision t = atanh(sqrtmz) / sqrtmz;

			t_200 = 1.  +  (1. + z) * t;
			t_220 = (-1.  +  (1. + z) * t) / z;
			t_201 = (1.  +  (z - 1.) * t) / z;
		}
		else if(fabs(z) <= delta)
		{
			precision z2 =  z * z;
			precision z3 = z2 * z;
			precision z4 = z3 * z;
			precision z5 = z4 * z;
			precision z6 = z5 * z;

			t_200 = 2. + 0.6666666666666667*z - 0.1333333333333333*z2 + 0.05714285714285716*z3 - 0.031746031746031744*z4 + 0.020202020202020193*z5 - 0.013986013986013984*z6;

			t_220 = 0.6666666666666667 - 0.1333333333333333*z + 0.05714285714285716*z2 - 0.031746031746031744*z3 + 0.020202020202020193*z4 - 0.013986013986013984*z5 + 0.010256410256410262*z6;

			t_201 = 1.3333333333333333 - 0.5333333333333333*z + 0.34285714285714286*z2 - 0.25396825396825395*z3 + 0.20202020202020202*z4 - 0.16783216783216784*z5 + 0.14358974358974358*z6;
		}
		else
		{
			printf("compute_F error: z = %lf is out of bounds\n", z);
			exit(-1);
		}

		I_200 += total_weight * t_200 * w;
		I_220 += total_weight * t_220 / w;
		I_201 += total_weight * t_201 / w;
	}

	I_200 *= prefactor;
	I_220 *= prefactor * aL2;
	I_201 *= prefactor * aT2 / 2.;

	F[0] = I_200 - Ea;
	F[1] = I_201 - PTa;
	F[2] = I_220 - PLa;
}


void compute_J(precision Ea, precision PTa, precision PLa, precision mass, precision * X, precision * F, precision ** J)
{
	precision lambda = X[0];
	precision aT = X[1];
	precision aL = X[2];

	precision lambda2 =  lambda * lambda;
  	precision lambda3 = lambda2 * lambda;

	precision aT2 = aT * aT;
	precision aL2 = aL * aL;
	precision aT2_minus_aL2 = aT2 - aL2;

	precision mbar = mass / lambda;
	precision mbar2 = mbar * mbar;

	precision lambda_aT3 = lambda * aT2 * aT;
  	precision lambda_aL3 = lambda * aL2 * aL;

  	precision prefactor = g * aT2 * aL * lambda2 * lambda3 / (4. * M_PI * M_PI);

  	precision t_200, t_201, t_220;
  	precision t_402, t_421, t_440;

    precision I_2001 = 0;
    precision I_2011 = 0;
    precision I_2201 = 0;
    precision I_402m1 = 0;
    precision I_421m1 = 0;
    precision I_440m1 = 0;

	for(int i = 0; i < pbar_pts; i++)														// gauss integration loop
	{
		precision pbar = pbar_root_a3[i];													// pbar roots / weights for a = 3 (a = n + s)
		precision pbar2 = pbar * pbar;
		precision Ebar = sqrt(pbar2 + mbar2);
		precision common_weight = pbar_weight_a3[i] * exp(pbar - Ebar);						// common weight factor

		precision w = sqrt(aL2  +  mbar2 / pbar2);
		precision z = aT2_minus_aL2 / (w * w);
		precision z2 =  z * z;

		if(z > delta)																		// compute hypergeometric functions
		{
			precision sqrtz = sqrt(z);
			precision t = atan(sqrtz) / sqrtz;

			t_200 = 1.  +  (1. + z) * t;
			t_220 = (-1.  +  (1. + z) * t) / z;
			t_201 = (1.  +  (z - 1.) * t) / z;
			t_402 = (3. * (z - 1.)  +  (z * (3.*z - 2.) + 3.) * t) / (4. * z2);
			t_421 = (3.  +  z  +  (1. + z) * (z - 3.) * t) / (4. * z2);
			t_440 = (-(3. + 5.*z) + 3. * (z + 1.) * (z + 1.) * t) / (4. * z2);
		}
		else if(z < -delta && z > -1.)
		{
			precision sqrtmz = sqrt(-z);
			precision t = atanh(sqrtmz) / sqrtmz;

			t_200 = 1.  +  (1. + z) * t;
			t_220 = (-1.  +  (1. + z) * t) / z;
			t_201 = (1.  +  (z - 1.) * t) / z;
			t_402 = (3. * (z - 1.)  +  (z * (3.*z - 2.) + 3.) * t) / (4. * z2);
			t_421 = (3.  +  z  +  (1. + z) * (z - 3.) * t) / (4. * z2);
			t_440 = (-(3. + 5.*z) + 3. * (z + 1.) * (z + 1.) * t) / (4. * z2);
		}
		else if(fabs(z) <= delta)
		{
			precision z3 = z2 * z;
			precision z4 = z3 * z;
			precision z5 = z4 * z;
			precision z6 = z5 * z;

			t_200 = 2. + 0.6666666666666667*z - 0.1333333333333333*z2 + 0.05714285714285716*z3 - 0.031746031746031744*z4 + 0.020202020202020193*z5 - 0.013986013986013984*z6;

			t_220 = 0.6666666666666667 - 0.1333333333333333*z + 0.05714285714285716*z2 - 0.031746031746031744*z3 + 0.020202020202020193*z4 - 0.013986013986013984*z5 + 0.010256410256410262*z6;

			t_201 = 1.3333333333333333 - 0.5333333333333333*z + 0.34285714285714286*z2 - 0.25396825396825395*z3 + 0.20202020202020202*z4 - 0.16783216783216784*z5 + 0.14358974358974358*z6;

			t_402 = 1.0666666666666667 - 0.4571428571428572*z + 0.3047619047619048*z2 - 0.23088023088023088*z3 + 0.1864801864801865*z4 - 0.15664335664335666*z5 + 0.13514328808446457*z6;

			t_421 = 0.2666666666666666 - 0.0761904761904762*z + 0.0380952380952381*z2 - 0.023088023088023088*z3 + 0.015540015540015537*z4 - 0.011188811188811189*z5 + 0.00844645550527904*z6;

				t_440 = 0.4 - 0.057142857142857106*z + 0.019047619047619063*z2 - 0.008658008658008663*z3 + 0.004662004662004657*z4 - 0.002797202797202792*z5 + 0.0018099547511312257*z6;
		}
		else
		{
			printf("compute J error: z = %lf is out of bounds\n", z);
			exit(-1);
		}

		I_2001 += Ebar * common_weight * t_200 * w;
		I_2011 += Ebar * common_weight * t_201 / w;
		I_2201 += Ebar * common_weight * t_220 / w;

		I_402m1 += pbar2 / Ebar * common_weight * t_402 / w;
    	I_421m1 += pbar2 / Ebar * common_weight * t_421 / w;
    	I_440m1 += pbar2 / Ebar * common_weight * t_440 / w;
	}

	I_2001 *= prefactor;
    I_2011 *= prefactor * aT2 / 2.;
    I_2201 *= prefactor * aL2;
    I_402m1 *= prefactor * aT2 * aT2 / 8.;
    I_421m1 *= prefactor * aT2 * aL2 / 2.;
    I_440m1 *= prefactor * aL2 * aL2;

    precision Eai = F[0] + Ea;					// compute Eai, PTai, PLai from F
	precision PTai = F[1] + PTa;
	precision PLai = F[2] + PLa;

    J[0][0] = I_2001 / lambda2;	  		J[1][0] = I_2011 / lambda2;				J[2][0] = I_2201 / lambda2;
    J[0][1] = 2. * (Eai + PTai) / aT;  	J[1][1] = 4. * I_402m1 / lambda_aT3;	J[2][1] = 2. * I_421m1 / lambda_aT3;
    J[0][2] = (Eai + PLai) / aL;	  	J[1][2] = I_421m1 / lambda_aL3;      	J[2][2] = I_440m1 / lambda_aL3;
}




precision line_backtrack(precision Ea, precision PTa, precision PLa, precision mass, precision * Xcurrent, precision * dX, precision g0, precision * F, precision * gradf)
{
	// I got this algorithm from Numerical Recipes in C
	//line_backtracking(&l, Ea, PTa, PLa, mass, X, dX, f, F, gradf); (for reference)
	// g0 = f

	precision l = 1.0;         						// default value for l
	precision alpha = 0.0001;   					// descent parameter (default was 0.0001)

	precision dX_abs = sqrt(dX[0] * dX[0]  +  dX[1] * dX[1]  +  dX[2] * dX[2]);

	precision X[3];

	for(int i = 0; i < 3; i++)
	{
		X[i] = Xcurrent[i] + l * dX[i];
	}

	compute_F(Ea, PTa, PLa, mass, X, F);			// update F at least once

	if(l < (tol_dX / dX_abs))						// what does this mean?
	{
		return l;
	}

	precision f = (F[0] * F[0]  +  F[1] * F[1]  +  F[2] * F[2]) / 2.;
	precision g1 = f;

	precision gprime0 = 0;

	for(int j = 0; j < 3; j++)
	{
		gprime0 += gradf[j] * dX[j];				// can I just pass this?
	}

	if(gprime0 >= 0)
	{
		// roundoff issue with gradient descent
		printf("line_backtracking error: descent derivative is not negative\n");
		exit(-1);
	}

	precision lroot;
	precision lprev;
	precision fprev;

	for(int n = 0; n < 20; n++)						// line search iterations (max is 20)
	{
		if(f <= (g0  +  l * alpha * gprime0))		// check for sufficient decrease in f
		{
			return l;
		}
		else if(n == 0)								// start with quadratic formula
		{
			lroot = - gprime0 / (2. * (g1 - g0 - gprime0));
		}
		else 										// cubic formula for the rest of iterations
		{
			precision a = ((g1  -  gprime0 * l  -  g0) / (l * l)  -  (fprev  -  gprime0 * lprev  -  g0) / (lprev * lprev)) / (l - lprev);
			precision b = (-lprev * (g1  -  gprime0 * l  -  g0) / (l * l)  +  l * (fprev  -  gprime0 * lprev  -  g0)  /  (lprev * lprev)) / (l - lprev);

			if(a == 0) 								// solve dg/dl = 0 at a = 0
			{
				lroot = - gprime0 / (2. * b);
			}
			else
			{
				precision z = b * b  -  3. * a * gprime0;

				if(z < 0)
				{
					lroot = l / 2.;
				}
				else if(b <= 0)
				{
					lroot = (-b + sqrt(z)) / (3. * a);
				}
				else
				{
					lroot = - gprime0 / (b + sqrt(z));
				}
			}

			lroot = fmin(lroot, 0.5 * l);
		}

		lprev = l;									// store current values for next iteration
		fprev = f;

		l = fmax(lroot, 0.5 * l);					// update l and f

		for(int i = 0; i < 3; i++)
		{
			X[i] = Xcurrent[i]  +  l * dX[i];
		}

		compute_F(Ea, PTa, PLa, mass, X, F);

		f = (F[0] * F[0]  +  F[1] * F[1]  +  F[2] * F[2]) / 2.;
	}

	return l;
}


void free_2D(precision ** M, int n)
{
	for(int i = 0; i < n; i++)
	{
		free(M[i]);
	}

    free(M);
}


aniso_variables find_anisotropic_variables(precision e, precision pl, precision pt, precision B, precision mass, precision lambda_0, precision aT_0, precision aL_0)
{
	const precision Ea = e - B;											// kinetic energy density
	const precision PTa = pt + B;										// kinetic transverse pressure
	const precision PLa = pl + B;										// kinetic longitudinal pressure

	if(Ea < 0 || PTa < 0 || PLa < 0)
	{
		printf("find_anisotropic_variables error: (E_a, PT_a, PL_a) = (%lf, %lf, %lf) is negative\n", Ea, PTa, PLa);

		aniso_variables variables;
		variables.lambda = lambda_0;
		variables.aT = aT_0;
		variables.aL = aL_0;
		variables.did_not_find_solution = 1;
		variables.number_of_iterations = 0;
		return variables;
	}

	precision lambda = lambda_0;
	precision aT = aT_0;
	precision aL = aL_0;

	precision X[3] = {lambda, aT, aL};									// current solution
	precision dX[3];							 				 		// dX iteration
  	precision F[3];  											 		// F(X)
  	precision dX_abs;             										// magnitude of dX
	precision F_abs;		      										// magnitude of F
  	precision f;		                                  	     		// f = F(X).F(X) / 2
  	precision gradf[3];										 			// gradient f

	precision **J = (precision**)malloc(3 * sizeof(precision*));        // J(X)

	for(int i = 0; i < 3; i++)
	{
		J[i] = (precision*)malloc(3 * sizeof(precision));
 	}

 	// precision tolmin = 1.0e-6;   									// tolerance for spurious convergence to local min of f = F.F/2 (what does this mean?)

	precision stepmax = 100.0;   										// scaled maximum step length allowed in line searches (what does this mean?)
	stepmax *= fmax(sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]), 3. * sqrt(3.));		// why??? (never used anyway..)

	compute_F(Ea, PTa, PLa, mass, X, F);								// compute F

	gsl_vector * x = gsl_vector_alloc(3);								// holds dX
    gsl_permutation * p = gsl_permutation_alloc(3);

	for(int n = 0; n < N_max; n++)										// Newton iteration loop
	{
	    compute_J(Ea, PTa, PLa, mass, X, F, J);							// compute J and f at X

	    f = (F[0] * F[0]  +  F[1] * F[1]  +  F[2] * F[2]) / 2.;

	    for(int i = 0; i < 3; i++)
    	{
    		gradf[i] = 0; 								  				// compute gradf_j = F_k * J_kj

    		for(int j = 0; j < 3; j++)
    		{
    			gradf[i] += F[j] * J[j][i];
    		}
    	}

	    double J_gsl[] = {J[0][0], J[0][1], J[0][2],
                  	  	  J[1][0], J[1][1], J[1][2],
                	  	  J[2][0], J[2][1], J[2][2]};					// Jacobian matrix in gsl format

        for(int i = 0; i < 3; i++) 										// change sign of F
	    {
	    	F[i] *= -1.;
	    }

	    int s; 															// solve matrix equations J.dX = -F
   		gsl_matrix_view A = gsl_matrix_view_array(J_gsl, 3, 3);
    	gsl_vector_view b = gsl_vector_view_array(F, 3);
       	gsl_linalg_LU_decomp(&A.matrix, p, &s);
       	gsl_linalg_LU_solve(&A.matrix, p, &b.vector, x);

       	for(int i = 0; i < 3; i++)
       	{
       		dX[i] = gsl_vector_get(x, i);
       	}

	    dX_abs = sqrt(dX[0] * dX[0]  +  dX[1] * dX[1]  +  dX[2] * dX[2]);

		if(dX_abs > stepmax)											// rescale dX if too large
		{
			for(int i = 0; i < 3; i++)
			{
				dX[i] *= stepmax / dX_abs;
			}
		}

		precision l = line_backtrack(Ea, PTa, PLa, mass, X, dX, f, F, gradf);	// compute partial step l and F(X + l.dX)

		for(int i = 0; i < 3; i++)										// update solution
	    {
	    	dX[i] *= l;
	    	X[i] += dX[i];
	    }

	    F_abs = sqrt(F[0] * F[0]  +  F[1] * F[1]  +  F[2] * F[2]);		// update convergence values
	   	dX_abs = sqrt(dX[0] * dX[0]  +  dX[1] * dX[1]  +  dX[2] * dX[2]);

		if(X[0] < 0 || X[1] < 0 || X[2] < 0)							// check if any variable goes negative
		{
			aniso_variables variables;
			variables.lambda = lambda_0;
			variables.aT = aT_0;
			variables.aL = aL_0;
			variables.did_not_find_solution = 1;
			variables.number_of_iterations = n + 1;

			free_2D(J, 3);
			gsl_permutation_free(p);
       		gsl_vector_free(x);

			return variables;
		}

		if(dX_abs <= tol_dX && F_abs <= tol_F)							// check for convergence
		{
			printf("%d\n", n+1);
			exit(-1);
			aniso_variables variables;
			variables.lambda = X[0];
			variables.aT = X[1];
			variables.aL = X[2];
			variables.did_not_find_solution = 0;
			variables.number_of_iterations = n + 1;

			free_2D(J, 3);
			gsl_permutation_free(p);
   			gsl_vector_free(x);

			return variables;
		}

	}	// newton iteration (n)

	//printf("find_anisotropic_variables error: number of iterations on (lambda_0, aT_0, aL_0) = (%lf, %lf, %lf) exceed maximum (%lf, %lf, %lf)\n", lambda_0, aT_0, aL_0, X[0], X[1], X[2]);

	aniso_variables variables;
	variables.lambda = lambda_0;
	variables.aT = aT_0;
	variables.aL = aL_0;
	variables.did_not_find_solution = 1;
	variables.number_of_iterations = N_max;

	free_2D(J, 3);
	gsl_permutation_free(p);
  	gsl_vector_free(x);

	return variables;
}


inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


void set_anisotropic_variables(const hydro_variables * const __restrict__ q, const precision * const __restrict__ e, precision * const __restrict__ lambda, precision * const __restrict__ aT, precision * const __restrict__ aL, lattice_parameters lattice, hydro_parameters hydro)
{
#ifdef ANISO_HYDRO
#ifdef LATTICE_QCD
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	double dx = lattice.lattice_spacing_x;
	double dy = lattice.lattice_spacing_y;

	int grid_size = nx * ny * nz;

	precision conformal_prefactor = hydro.conformal_eos_prefactor;

	int number_of_failed_solutions = 0;

	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			double y = (j - 2. - (ny - 1.)/2.) * dy;

			for(int i = 2; i < nx + 2; i++)
			{
				double x = (i - 2. - (nx - 1.)/2.) * dx;

				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				precision e_s = e[s];
				precision pl = q[s].pl;
				precision pt = q[s].pt;
				precision b = q[s].b;

				equation_of_state_new eos(e_s, conformal_prefactor);
				precision T = eos.T;
				precision mass = T * eos.z_quasi();

				precision lambda_prev = lambda[s];
				precision aT_prev = aT[s];
				precision aL_prev = aL[s];

				aniso_variables X_s = find_anisotropic_variables(e_s, pl, pt, b, mass, lambda_prev, aT_prev, aL_prev);

				lambda[s] = X_s.lambda;			// update anisotropic variables
				aT[s] = X_s.aT;
				aL[s] = X_s.aL;
				number_of_failed_solutions += X_s.did_not_find_solution;

				if(X_s.did_not_find_solution)
				{
					aniso_regulation[s] += 1;
				}

			}
		}
	}

	// if(number_of_failed_solutions > 0)
	// {
	// 	printf("failed solutions = %d\n", number_of_failed_solutions);
	// }

#endif
#endif
}







