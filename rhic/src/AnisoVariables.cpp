#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include "../include/TransportAnisoNonconformal.h"
#include "../include/TransportAniso.h"
#include "../include/AnisoVariables.h"
#include "../include/Precision.h"


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
			printf("Error: z = %lf is out of bounds\n", z);
			exit(-1);
		}

		I_200 += total_weight * t_200 * w;
		I_220 += total_weight * t_220 / w;
		I_201 += total_weight * t_201 / w;
	}

	I_200 *= prefactor;
	I_220 *= prefactor * aL2;
	I_201 *= prefactor * aT2 / 2.;

	// printf("I_200 = %lf\n", I_200);
	// printf("I_220 = %lf\n", I_220);
	// printf("I_201 = %lf\n\n", I_201);
	//exit(-1);

	F[0] = I_200 - Ea;
	F[1] = I_201 - PTa;
	F[2] = I_220 - PLa;
}




void compute_J(precision Ea, precision PTa, precision PLa, precision mass, precision * X, precision * dX, precision * F, precision * Fcurrent, precision ** J, precision ** Jcurrent, jacobian method)
{
    switch(method)
    {
    	case newton:	// exact Jacobian
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
					printf("Error: z = %lf is out of bounds\n", z);
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

		    // printf("J_00 = %lf\n", J[0][0]);
		    // printf("J_10 = %lf\n", J[1][0]);
		    // printf("J_20 = %lf\n\n", J[2][0]);
		    // printf("J_01 = %lf\n", J[0][1]);
		    // printf("J_11 = %lf\n", J[1][1]);
		    // printf("J_21 = %lf\n\n", J[2][1]);
		    // printf("J_02 = %lf\n", J[0][2]);
		    // printf("J_12 = %lf\n", J[1][2]);
		    // printf("J_22 = %lf\n", J[2][2]);
		    // exit(-1);

    		break;
    	}
    	case broyden:		// approximate Jacobian
    	{
    		precision dX_abs = sqrt(dX[0] * dX[0]  +  dX[1] * dX[1]  +  dX[2] * dX[2]);
    		precision dX_unit[3] = {0, 0, 0};
    		precision Jcurrent_dX[3] = {0, 0, 0};

    		for(int i = 0; i < 3; i++)
			{
				dX_unit[i] = dX[i] / dX_abs;

				for(int j = 0; j < 3; j++)
				{
					Jcurrent_dX[i] += Jcurrent[i][j] * dX[j];
				} 
			}
    		for(int i = 0; i < 3; i++)
    		{
    			for(int j = 0; j < 3; j++)
    			{
    				J[i][j] = Jcurrent[i][j]  +  (F[i] - Fcurrent[i] - Jcurrent_dX[i]) * dX_unit[j] / dX_abs;	// do I get any numerical error from this?
    			}
    		}

    		break;
    	}
    	default:
    	{
    		printf("compute_J error: fix Jacobian method\n");
    		exit(-1);
    	}
    }
}




//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//                    LUP LINEAR EQUATION SOLVER                    ::
//                                                                  ::
//     Solves linear equation Ax = b using LU decomposition         ::
//	   with implicit partial pivoting (LUP). To directly solve      ::
//     for Ax = b, run these two functions concurrently:            ::
//                                                                  ::
//            LUP_decomposition              LUP_solve              ::
//                                                                  ::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


void LUP_decomposition(double ** A, int n, int * pvector)
{
	// takes A and decomposes it into LU (with row permutations P)
	// A = n x n matrix; function does A -> PA = LU; (L,U) of PA stored in same ** array
	// n = size of A
	// pvector = permutation vector; set initial pvector[i] = i (to track implicit partial pivoting)
	// 			 function updates pvector if there are any rows exchanges in A (to be used on b in LUP_solve)

	int i;	   // rows
	int j;	   // columns
	int k;     // dummy matrix index
	int imax;  // pivot row index
	double big;
	double sum;
	double temp;
	double implicit_scale[n];

	// Initialize permutation vector
	// to default no-pivot values
	for(i = 0; i < n; i++)
	{
		pvector[i] = i;
	}
	// Implicit scaling info. for A
	for(i = 0; i < n; i++)
	{
		big = 0.0;
		for(j = 0; j < n; j++)
		{
			temp = fabs(A[i][j]);
			if(temp > big)
			{
				big = temp;  // update biggest element in the ith row
			}
		}
		if(big == 0.0)
		{
			printf("Singular matrix in the routine");
			break;
		}
		implicit_scale[i] = 1.0 / big;  // store implicit scale of row i (will be used later)
	}
	// LU Decomposition
	for(j = 0; j < n; j++)
	{
		// loop rows i = 1 to j - 1
		for(i = 0; i < j; i++)
		{
			sum = A[i][j];
			for(k = 0; k < i; k++)
			{
				sum -= A[i][k] * A[k][j];
			}
			A[i][j] = sum;  // update U[i][j] elements
		}

		big = 0.0;          // initialize search for the largest normalized pivot in j column

		// loop through rows i = j to n-1
		for(i = j; i < n; i++)
		{
			sum = A[i][j];
			for(k = 0; k < j; k++)
			{
				sum -= A[i][k] * A[k][j];
			}
			A[i][j] = sum;   // update U[j][j] and L[i][j] elements (no division in L yet until the pivot determined)

			temp = implicit_scale[i] * fabs(sum);
			if(temp >= big)
			{
				big = temp;  // searchs for the biggest normalized member (imax) in column j
				imax = i;	 // implicit scale * A[i][j] normalizes each row entry i in column j before comparing
			}			     // implicit scale different for each i, that's why it's important
		}
		if(j != imax)
		{
			// then exchange rows j and imax of A
			for(k = 0; k < n; k++)
				{
					temp = A[imax][k];
					A[imax][k] = A[j][k];
					A[j][k] = temp;
				}
			implicit_scale[imax] = implicit_scale[j];   // interchange scaling
		}

		pvector[j] = imax;   		  // update permutation vector keeps track of pivot indices

		if(A[j][j] == 0.0)
		{
			A[j][j] = 1.e-16;         // matrix is singular
		}
		if(j != n-1)                  // there is no L[n,n] element
		{
			temp = 1.0 / A[j][j];     // divide L[i,j] elements by the pivot
			for(i = j+1; i < n; i++)
			{
				A[i][j] *= temp;
			}
		}
	}
}


void LUP_solve(double ** PA, int n, int * pvector, double b[])
{
	// input vector b is transformed to the solution x  of Ax = b
	// PA is the input permutated matrix from LUP_decomposition (will not be updated here)
	// input pvector comes from LUP_decomposition (generally not default); used to switch rows of b
	int i;       // rows
	int j;       // columns
	int m = -1;  // used to skip first few for loops (j) involving a zero (priorly permutated) b[j] element (pointless to multiply by 0)
	int ip;      // assigned permutation row index pvector[i]
	double sum;
	// Forward substitution routine for Ly = b
	for(i = 0; i < n; i++)
	{
		ip = pvector[i];         // permute b accordingly
		sum = b[ip];             // starting value given right b[ip]
		b[ip] = b[i];            // switch value of b[ip] to b[i]
		if(m != -1)			     // if m stays -1, skip for loop until it's been assigned a j_min >= 0
		{
			for(j = m; j <= i-1; j++)
			{
				sum -= PA[i][j] * b[j];    // forward iteration
			}
		}
		else if(sum)
		{
			m = i;               // once encounter nonzero b[ip], m = j_min >= 0 from now on
		}
		b[i] = sum;              // update y[i] and store in b
	}
	// Backward substitution routine for Ux = y
	for(i = n-1; i >= 0; i--)
	{
		sum = b[i];
		for(j = i+1; j < n; j++)
		{
			sum -= PA[i][j] * b[j];       // backward iteration
		}
		b[i] = sum / PA[i][i];            // update x[i] and store in b (FINAL SOLUTION)
	}
}


void free_2D(precision ** M, int n)
{
	for (int i = 0; i < n; i++) free(M[i]);
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
		exit(-1);
	}

	// let's keep the matrix solver for now

	precision lambda = lambda_0;
	precision aT = aT_0;
	precision aL = aL_0;

	precision X[3] = {lambda, aT, aL};									// updated solution 
	precision Xcurrent[3] = {lambda, aT, aL};			     			// current solution 
	precision dX[3];							 				 		// dX iteration
  	precision F[3];  											 		// F(X) 
  	precision Fcurrent[3];  									 		// F(Xcurrent)
  	precision dX_abs;             										// magnitude of dX
	precision F_abs;		      										// magnitude of F
  	precision f;		                                  	     		// f = F(X).F(X) / 2
  	precision fcurrent;										 			// f = F(Xcurrent).F(Xcurrent) / 2
  	precision gradf[3];										 			// gradient f at Xcurrent

	precision **J = (precision**)malloc(3*sizeof(precision*));        	// J(X)
	precision **Jcurrent = (precision**)malloc(3*sizeof(precision*)); 	// J(Xcurrent)			(maybe call them prev instead of current)
	for(int k = 0; k < 3; k++)
	{
		J[k] = (precision*)malloc(3*sizeof(precision));
 		Jcurrent[k] = (precision*)malloc(3*sizeof(precision));
 	}
	
 	int pvector[3];									  		 			// permutation vector for LUP solver

 	jacobian method = newton;                               			// jacobian method (start with newton)
 	precision l = 1; 		  											// default partial step parameter

 	// not sure
 	precision tolmin = 1.0e-6;   										// tolerance for spurious convergence to local min of f = F.F/2 (what does this mean?)
	precision stepmax = 100;   											// scaled maximum step length allowed in line searches (what does this mean?)

	stepmax *= fmax(sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]), 3. * sqrt(3.0));		// why??? (never used anyway..)


	compute_F(Ea, PTa, PLa, mass, Xcurrent, F);							// first compute F for the first iteration					

	// printf("F_0 = %lf\n", F[0]);
	// printf("F_1 = %lf\n", F[1]);
	// printf("F_2 = %lf\n\n", F[2]);
	

	printf("Start search for anisotropic variables\n");					// there's a lot about this finder that I don't remember

	for(int i = 0; i < N_max; i++)										// Newton-Broyden iteration loop
	{
		if(i > 0)														// compute at least one Newton iteration before checking for convergence
		{
			//if(dX_abs <= tol_dX & F_abs <= tol_F)						// was this a bug???
			if(dX_abs <= tol_dX && F_abs <= tol_F)						// why do I need both, don't I just need F?
			{
				//cout << i << " ";

				aniso_variables variables;
				variables.lambda = X[0];
				variables.aT = X[1];
				variables.aL = X[2];

				free_2D(J, 3);											// free allocated memory 
				free_2D(Jcurrent, 3);									// can I do something else? 

				return variables; 
			}

			method = broyden;											// default method is Broyden for i > 0
		}

		// default J and f at Xcurrent
	    compute_J(Ea, PTa, PLa, mass, Xcurrent, dX, F, Fcurrent, J, Jcurrent, method);		
	    fcurrent = (F[0] * F[0]  +  F[1] * F[1]  +  F[2] * F[2]) / 2.;

	    for(int j = 0; j < 3; j++)
    	{
    		Fcurrent[j] = F[j];											// store current F,  default J and compute
    		gradf[j] = 0; 												// gradient of f at Xcurrent (gradf_j = F_k * J_kj)		

    		for(int k = 0; k < 3; k++)
    		{
    			Jcurrent[j][k] = J[j][k];
    			gradf[j] += F[k] * J[k][j];
    		}
    	}

  		// printf("gradf_0 = %lf\n", gradf[0]);
		// printf("gradf_1 = %lf\n", gradf[1]);
		// printf("gradf_2 = %lf\n", gradf[2]);

		// solve matrix equation: J * dX = - F (stick with LUP routine for now...)
		// would it be more effective if solved equation in dimensionless units?
		// I could try normalizing this equation later
		// for(int j = 0; j < 3; j++) 									// normalizing equation doesn't appear to make a difference
		// {	
		// 	F[j] /= pow(lambda, 4);

		// 	for(int k = 0; k < 3; k++)
  //   		{
  //   			J[j][k] /= pow(lambda, 4);
  //   		}
		// }

	    for(int j = 0; j < 3; j++) 										// change sign of F first
	    {
	    	F[j] *= -1;  
	    }

	    LUP_decomposition(J, 3, pvector);          						// LUP of J now stored in J (pvector also updated)
	    LUP_solve(J, 3, pvector, F);               						// F now stores dX

	    for(int j = 0; j < 3; j++) 
	    {
	    	 dX[j] = F[j];   											// load full Newton/Broyden iteration dX
	    }
	   
		// printf("dX_0 = %lf\n", dX[0]);
		// printf("dX_1 = %lf\n", dX[1]);
		// printf("dX_2 = %lf\n", dX[2]);

	    // rescale dX if too large
	    dX_abs = sqrt(dX[0] * dX[0]  +  dX[1] * dX[1]  +  dX[2] * dX[2]);	

		if(dX_abs > stepmax)
		{
			printf("dX is too large, rescale it by stepmax\n");			// I don't think this ever happens...

			for(int j = 0; j < 3; j++)
			{
				dX[j] *= stepmax / dX_abs;
			}
		}

		// compute gradient descent derivative
		precision gprime0 = 0;		

		for(int j = 0; j < 3; j++)										
		{
			gprime0 += gradf[j] * dX[j];
		}

		switch(method)
		{
			case newton:
			{
				if(gprime0 >= 0) 
				{
					printf("find_anisotropic_variables error: gradient descent = %lf has wrong sign\n", gprime0);
    				exit(-1);
				}
			}
			case broyden:
			{
				break;
			}
			default:
			{
				printf("find_anisotropic_variables error: gradient descent has no Jacobian method\n");
    			exit(-1);
			}
		}

		printf("gprime_0 = %lf\n", gprime0);
		exit(-1);


		switch(method)													// update solution X or redo iteration with newton
		{
			case newton:												// compute line backtracking step
			{
				// line backtracking algorithm updates l and F(Xcurrent + l.dX)
				// 
				linebacktrack(&l, Ea, PTa, PLa, mass, Xcurrent, dX, fcurrent, F, gradf, tol_dX);

			    for(int j = 0; j < 3; j++)
			    {
			    	X[j] = Xcurrent[j]  +  l * dX[j];
			    }

				break;
			}
			case broyden:
			{
				for(int j = 0; j < 3; j++)								// default Broyden iteration (l = 1)
				{
					X[j] = Xcurrent[j] + dX[j];
				}

				compute_F(Ea, PTa, PLa, mass, X, F);					// compute default values for F and f at updated X
				f = (F[0] * F[0]  +  F[1] * F[1]  +  F[2] * F[2]) / 2.;


				if(f > (fcurrent  -  0.2 * fabs(gprime0)))				// if Broyden step doesn't decrease f sufficiently
				{														// compute exact Jacobian and redo iteration with newton
					//cout << "Hey" << endl;
					// reset F
			    	for(int k = 0; k < n; k++)
			    	{
			    		F[k] = Fcurrent[k];
			    	}

			    	jtype = newton;

			    	// reinitialize J as exact
			    	calcJ(Ea, PTa, PLa, thermal_mass, Xcurrent, dX, F, Fcurrent, J, Jcurrent, jtype, n);


			    	// store J and calculate gradf
			    	for(int k = 0; k < n; k++)
			    	{
			    		for(int j = 0; j < n; j++)
			    		{
			    			Jcurrent[k][j] = J[k][j];
			    		}
				    	gradf_sum = 0.0; // clear sum
				    	for(int m = 0; m < n; m++) gradf_sum += F[m] * J[m][k];
				    	gradf[k] = gradf_sum;
				    }

				    // resolve matrix equation
				    for(int k = 0; k < n; k++) F[k] = - F[k];
				    LUP_decomposition(J, n, pvector);
				    LUP_solve(J, n, pvector, F);
				    for(int k = 0; k < n; k++) dX[k] = F[k];


				    dXabs = sqrt(dX[0]*dX[0]+dX[1]*dX[1]+dX[2]*dX[2]);
					if(dXabs > stepmax)
					{
						for(int k = 0; k < n; k++) dX[k] *= (stepmax/dXabs);
					}


					// default Newton iteration
					for(int k = 0; k < n; k++) X[k] = Xcurrent[k] + dX[k];


					// line backtracking
					linebacktrack(&l, Ea, PTa, PLa, thermal_mass, Xcurrent, dX, fcurrent, F, gradf, toldX, n);


					// redo update for X
				    for(int k = 0; k < n; k++)
				    {
					    X[k] = Xcurrent[k] + l*dX[k];
				    }
				}
				break;
			}
			default:
				throw "correct line backtracking routine\n";
		}








	}	// newton iteration (i)

	// I should make a counter that tells me how many times I use newton-redo or broyden (or make a string)

	printf("find_anisotropic_variables error: number of iterations exceed maximum\n");

	aniso_variables variables;			// does algorithm ever get stuck searching (e.g. hit local min)? 
	variables.lambda = X[0];
	variables.aT = X[1];
	variables.aL = X[2];

	free_2D(J, 3);										
	free_2D(Jcurrent, 3);	

	return variables;
}











