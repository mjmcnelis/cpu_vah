/*
 * InitialConditions.c
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#include <math.h> // for math functions
#include <stdio.h> // for printf
#include <stdlib.h> //TEMP
#include <iostream>
#include <algorithm>    // for max

#include "../include/InitialConditions.h"
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"
#include "../include/Parameters.h"
#include "../include/OpticalGlauber.h"
#include "../include/MCGlauber.h"
#include "../include/EquationOfState.h"

using namespace std;

#define THETA_FUNCTION(X) ((double)X < (double)0 ? (double)0 : (double)1)

const double hbarc = 0.197326938;


inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}

// initialize conserved variables
void set_initial_conserved_variables(double t, int nx, int ny, int nz)
{
	// loop over the physical grid points
	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				// get the variables (computed from initialConditions(...))
				precision e_s = e[s];

				precision ux = u->ux[s];
				precision uy = u->uy[s];
				precision un = u->un[s];
				precision ut = u->ut[s];

				precision pl = q->pl[s];
				precision pt = q->pt[s];

				precision utperp = sqrt(1.0  +  ux * ux  +  uy * uy);

				// z^mu components
				precision zt = t * un / utperp;
				precision zn = ut / t / utperp;

				// only need the time components
				precision pitt = 0.0;
				precision pitx = 0.0;
				precision pity = 0.0;
				precision pitn = 0.0;

				precision Wt = 0.0;
				precision Wx = 0.0;
				precision Wy = 0.0;
				precision Wn = 0.0;

			#ifdef PIMUNU
				pitt = q->pitt[s];
				pitx = q->pitx[s];
				pity = q->pity[s];
				pitn = q->pitn[s];
			#endif
			#ifdef WTZMU
				Wt = q->WtTz[s];
				Wx = q->WxTz[s];
				Wy = q->WyTz[s];
				Wn = q->WnTz[s];
			#endif

				// initialize the time components of Tmunu
				q->ttt[s] = (e_s + pt) * ut * ut  -   pt  +  (pl - pt) * zt * zt  +  2.0 * Wt * zt  +  pitt;
				q->ttx[s] = (e_s + pt) * ut * ux  +  Wx * zt  +  pitx;
				q->tty[s] = (e_s + pt) * ut * uy  +  Wy * zt  +  pity;
				q->ttn[s] = (e_s + pt) * ut * un  +  (pl - pt) * zt * zn  +  (Wt * zn  +  Wn * zt)  +  pitn;
			}
		}
	}
}

// initialize viscous part of Tmunu to zero
void set_equilibrium_initial_condition(int nx, int ny, int nz)
{
	// loop over physical grid points
	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				precision e_s = e[s];
				precision p_s = equilibriumPressure(e_s);

				q->pl[s] = p_s;		// set (pl, pt) to equilibium pressure
				q->pt[s] = p_s;

			#ifdef PIMUNU
		  		q->pitt[s] = 0.0;
		  		q->pitx[s] = 0.0;
		  		q->pity[s] = 0.0;
		  		q->pitn[s] = 0.0;
		  		q->pixx[s] = 0.0;
		  		q->pixy[s] = 0.0;
		  		q->pixn[s] = 0.0;
		  		q->piyy[s] = 0.0;
		  		q->piyn[s] = 0.0;
		  		q->pinn[s] = 0.0;
			#endif
			#ifdef W_TZ_MU
		  		q->WtTz[s] = 0.0;
		  		q->WxTz[s] = 0.0;
		  		q->WyTz[s] = 0.0;
		  		q->WnTz[s] = 0.0;
			#endif
			}
		}
	}
}


// Constant initial energy density profile (Bjorken)
void set_Bjorken_energy_density_and_flow_profile(int nx, int ny, int nz, void * initCondParams)
{
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	double T0 = initCond->initialCentralTemperatureGeV;		// central temperature (GeV)
	double e0 = equilibriumEnergyDensity(T0 / hbarc);		// energy density

	// loop over physical grid points
	for(int i = 2; i < nx + 2; i++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int k = 2; k < nz + 2; k++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				e[s] = (precision) e0;

				u->ut[s] = 1.0;
				u->ux[s] = 0.0;
				u->uy[s] = 0.0;
				u->un[s] = 0.0;

				// also initialize up = u
				up->ut[s] = u->ut[s];
				up->ux[s] = u->ux[s];
				up->uy[s] = u->uy[s];
				up->un[s] = u->un[s];
			}
		}
	}
}


// Longitudinal energy density profile function eL(eta)
void longitudinal_Energy_Density_Profile(double * const __restrict__ eL, int nz, double dz, void * initCondParams)
{
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	double etaFlat = initCond->rapidityMean;
	double etaVariance = initCond->rapidityVariance;

	// profile along eta direction is a smooth plateu that exponentially decays when |eta| > etaFlat
	for(int k = 0; k < nz; k++)
	{
		double eta = (k - (nz - 1.0)/2.0) * dz;				// physical eta points
		double etaScaled = fabs(eta)  -  0.5 * etaFlat;

		eL[k] = exp(- 0.5 * etaScaled * etaScaled / etaVariance * THETA_FUNCTION(etaScaled));
	}
}



// Optical or MC glauber initial energy density profile
void set_Glauber_energy_density_and_flow_profile(int nx, int ny, int nz, double dx, double dy, double dz, void * initCondParams)
{
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int initialConditionType = initCond->initialConditionType;

	double T0 = initCond->initialCentralTemperatureGeV;		// central temperature (GeV)
	double e0 = equilibriumEnergyDensity(T0 / hbarc);		// energy density scale factor

	double eT[nx * ny];		// normalized transverse profile
	double eL[nz];			// normalized longitudinal profile

	// compute normalized transverse and longitudinal profiles
	if(initialConditionType == 3)
	{
		Optical_Glauber_energy_density_transverse_profile(eT, nx, ny, dx, dy, initCondParams);
	}
	else if(initialConditionType == 4)
	{
		MC_Glauber_energy_density_transverse_profile(eT, nx, ny, dx, dy, initCondParams);
	}
	longitudinal_Energy_Density_Profile(eL, nz, dz, initCondParams);


	double e_min = 1.e-3;	// minimum energy density

	// loop over physical grid points
	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				double e_s = max(e_min, e0 * eT[i - 2 + (j - 2) * nx] * eL[k - 2]);

				e[s] = (precision) e_s;

				u->ut[s] = 1.0;
				u->ux[s] = 0.0;
				u->uy[s] = 0.0;
				u->un[s] = 0.0;

				// also initialize up = u
				up->ut[s] = u->ut[s];
				up->ux[s] = u->ux[s];
				up->uy[s] = u->uy[s];
				up->un[s] = u->un[s];
			}
		}
	}
}



// Initial conditions for the Gubser flow test
void set_Gubser_energy_density_and_flow_profile(int nx, int ny, int nz, double dx, double dy, double dz, void * initCondParams, void * hydroParams)
{
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
	struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;

	double t = hydro->tau_initial;								// initial longitudinal proper time
	double T0 = initCond->initialCentralTemperatureGeV / hbarc;	// central temperature (fm)

	double q = 1.0;												// inverse length size
	double q2 = q * q;
	double q4 = q2 * q2;
	double t2 = t * t;

	// normalize Gubser temperature profile s.t. central temperature = T0
	double T0_hat = T0 * t * pow((1.0 + q2 * t2) / (2.0 * q * t), 2./3.);

	// loop over physical grid points
	for(int i = 2; i < nx + 2; i++)
	{
		double x = (i - 2.0 - (nx - 1.0)/2.0) * dx;

		for(int j = 2; j < ny + 2; j++)
		{
			double y = (j - 2.0 - (ny - 1.0)/2.0) * dy;

			double r = sqrt(x * x  +  y * y);
			double r2 = r * r;

			double kappa = atanh(2.0 * q2 * t * r / (1.0  +  q2 * (t2 + r2)));

			// temperature
			double T = (T0_hat / t) * pow(2.0 * q * t, 2./3.) / pow(1.0  +  2.0 * q2 * (t2 + r2)  +  q4 * (t2 - r2) * (t2 - r2), 1./3.);

			double e_s = equilibriumEnergyDensity(T);

			double ux_s = sinh(kappa) * x / r;
			double uy_s = sinh(kappa) * y / r;
			double ut_s = sqrt(1.0  +  ux_s * ux_s  +  uy_s * uy_s);

			for(int k = 2; k < nz + 2; k++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				e[s] = (precision) e_s;

				u->ut[s] = (precision) ut_s;
				u->ux[s] = (precision) ux_s;
				u->uy[s] = (precision) uy_s;
				u->un[s] = 0.0;

				// also initialize up = u
				up->ut[s] = (precision) ut_s;
				up->ux[s] = (precision) ux_s;
				up->uy[s] = (precision) uy_s;
				up->un[s] = 0.0;
			}
		}
	}
}


// Initial conditions for the T^{\tau\mu}, energy density, equilibrium pressure, fluid velocity and viscous pressures
//////////////////////////////////////
//	1 - Bjorken
//	2 - Gubser
//	3 - Optical Glauber
//	4 - MC Glauber
//////////////////////////////////////
void set_initial_conditions(double t, void * latticeParams, void * initCondParams, void * hydroParams)
{
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	int initialConditionType = initCond->initialConditionType;
	printf("Initial conditions    = ");

	switch(initialConditionType)
	{
		case 1:
		{
			printf("Bjorken ");
		#ifndef CONFORMAL_EOS
			printf("QCD ");
		#else
			printf("Conformal ");
		#endif
			printf("(fluid velocity and viscous pressures initialized to zero)\n");

			set_Bjorken_energy_density_and_flow_profile(nx, ny, nz, initCondParams);
			set_equilibrium_initial_condition(nx, ny, nz);
			set_initial_conserved_variables(t, nx, ny, nz);

			break;
		}
		case 2:
		{
			printf("Gubser Conformal (viscous pressures initialized to zero)\n");
		#ifndef CONFORMAL_EOS
			printf("\nsetInitialConditions error: CONFORMAL_EOS not defined in /rhic/include/EquationOfState.h, exiting...\n");
			exit(-1);
		#endif

			set_Gubser_energy_density_and_flow_profile(nx, ny, nz, dx, dy, dz, initCondParams, hydroParams);
			set_equilibrium_initial_condition(nx, ny, nz);
			set_initial_conserved_variables(t, nx, ny, nz);

			break;
		}
		case 3:
		{
			printf("Optical Glauber (fluid velocity and viscous pressures initialized to zero)\n");

			set_Glauber_energy_density_and_flow_profile(nx, ny, nz, dx, dy, dz, initCondParams);
			set_equilibrium_initial_condition(nx, ny, nz);
			set_initial_conserved_variables(t, nx, ny, nz);

			break;
		}

		case 4:
		{
			printf("MC Glauber");
			set_Glauber_energy_density_and_flow_profile(nx, ny, nz, dx, dy, dz, initCondParams);
			printf("(fluid velocity and viscous pressures initialized to zero)\n");
			set_equilibrium_initial_condition(nx, ny, nz);
			set_initial_conserved_variables(t, nx, ny, nz);

			break;
		}
		case 5:
		{
			printf("Trento + F.S.");
			printf("\t(nothing here yet!)\n");
			exit(-1);

			// for setting trento + fs initial conditions
			//  - reconstruct everything? (no need for constructing Tmunu; what's the f.s. file format?)
			//	- I need to initialize the fluid velocity u
			//	  and previous fluid velocity up (either up = u or somehow get the previous time step)
			break;
		}
		default:
		{
			printf("\nsetInitialConditions error: initial condition type not defined, exiting...\n");
			exit(-1);
		}
	}
}


