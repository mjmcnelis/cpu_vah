/*
 * InitialConditions.c
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#include <stdlib.h> //TEMP
#include <stdio.h> // for printf
#include <math.h> // for math functions
#include <cmath>
#include <iostream>
#include <algorithm>    // for max

#include "../include/DynamicalVariables.h"
#include "../include/InitialConditions.h"
#include "../include/Precision.h"
#include "../include/Parameters.h"
#include "../include/OpticalGlauber.h"
#include "../include/MCGlauber.h"
#include "../include/AnisoGubser.h"
#include "../include/EquationOfState.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>

using namespace std;

#define THETA_FUNCTION(X) ((double)X < (double)0 ? (double)0 : (double)1)

const double hbarc = 0.197326938;


inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}

// initialize (ttt, ttx, tty, ttn)
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

				precision e_s = e[s];

				precision pl = q->pl[s];

			#if (PT_MATCHING == 1)
				precision pt = q->pt[s];
			#else
				precision pt = 0.5 * (e_s - pl);
			#endif

				precision ut = u->ut[s];
				precision ux = u->ux[s];
				precision uy = u->uy[s];
				precision un = u->un[s];

				precision utperp = sqrt(1.0  +  ux * ux  +  uy * uy);

				// z^mu components
				precision zt = t * un / utperp;
				precision zn = ut / t / utperp;

				// only need the time components
			#ifdef PIMUNU
				precision pitt = q->pitt[s];
				precision pitx = q->pitx[s];
				precision pity = q->pity[s];
				precision pitn = q->pitn[s];
			#else
				precision pitt = 0.0;
				precision pitx = 0.0;
				precision pity = 0.0;
				precision pitn = 0.0;
			#endif
			#ifdef WTZMU
				precision WtTz = q->WtTz[s];
				precision WxTz = q->WxTz[s];
				precision WyTz = q->WyTz[s];
				precision WnTz = q->WnTz[s];
			#else
				precision WtTz = 0.0;
				precision WxTz = 0.0;
				precision WyTz = 0.0;
				precision WnTz = 0.0;
			#endif

				// initialize the time components of Tmunu
				q->ttt[s] = (e_s + pt) * ut * ut  -   pt  +  (pl - pt) * zt * zt  +  2.0 * WtTz * zt  +  pitt;
				q->ttx[s] = (e_s + pt) * ut * ux  +  WxTz * zt  +  pitx;
				q->tty[s] = (e_s + pt) * ut * uy  +  WyTz * zt  +  pity;
				q->ttn[s] = (e_s + pt) * ut * un  +  (pl - pt) * zt * zn  +  (WtTz * zn  +  WnTz * zt)  +  pitn;
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

			#if (PT_MATCHING == 1)
				q->pt[s] = p_s;
			#endif

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

				e[s] = e0;

				u->ut[s] = 1.0;
				u->ux[s] = 0.0;
				u->uy[s] = 0.0;
				u->un[s] = 0.0;

				// also initialize up = u
				up->ut[s] = 1.0;
				up->ux[s] = 0.0;
				up->uy[s] = 0.0;
				up->un[s] = 0.0;
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

	// loop over physical grid points
	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				double e_s = max(E_MIN, e0 * eT[i - 2 + (j - 2) * nx] * eL[k - 2]);

				e[s] = (precision) e_s;

				// default initial flow to zero
				u->ut[s] = 1.0;
				u->ux[s] = 0.0;
				u->uy[s] = 0.0;
				u->un[s] = 0.0;

				up->ut[s] = 1.0;
				up->ux[s] = 0.0;
				up->uy[s] = 0.0;
				up->un[s] = 0.0;
			}
		}
	}
}



// initial conditions for the ideal Gubser flow test
void set_ideal_gubser_energy_density_and_flow_profile(int nx, int ny, int nz, double dt, double dx, double dy, double dz, void * initCondParams, void * hydroParams)
{
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
	struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;

	double t = hydro->tau_initial;								// initial longitudinal proper time
	double T0 = initCond->initialCentralTemperatureGeV / hbarc;	// central temperature (fm)

	double q  = 1.0;											// inverse length size (hard coded)
	double q2 = q * q;
	double q4 = q2 * q2;
	double t2 = t * t;

	// normalize Gubser ideal temperature profile so that central temperature = T0
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

			double kappa   = atanh(2.0 * q2 * t * r / (1.0  +  q2 * (t2  +  r2)));
			double kappa_p = atanh(2.0 * q2 * (t - dt) * r / (1.0  +  q2 * ((t - dt) * (t - dt)  +  r2)));

			if(std::isnan(kappa) || std::isnan(kappa_p))
			{
				printf("Gubser initial conditions error: (kappa, kappa_p) = (%lf, %lf)\n", kappa, kappa_p);
				exit(-1);
			}

			// temperature profile
			double T = (T0_hat / t) * pow(2.0 * q * t, 2./3.) / pow(1.0  +  2.0 * q2 * (t2 + r2)  +  q4 * (t2 - r2) * (t2 - r2), 1./3.);

			double e_s = max(E_MIN, equilibriumEnergyDensity(T));

			double ux = sinh(kappa) * x / r;
			double uy = sinh(kappa) * y / r;

			double ux_p = sinh(kappa_p) * x / r;
			double uy_p = sinh(kappa_p) * y / r;

			if(std::isnan(ux)) ux = 0.0;		// remove x/r = 0/0 nan
			if(std::isnan(uy)) uy = 0.0;		// remove y/r = 0/0 nan

			if(std::isnan(ux_p)) ux_p = 0.0;	// remove x/r = 0/0 nan
			if(std::isnan(uy_p)) uy_p = 0.0;	// remove y/r = 0/0 nan

			double ut   = sqrt(1.0  +  ux * ux  +  uy * uy);
			double ut_p = sqrt(1.0  +  ux_p * ux_p  +  uy_p * uy_p);

			for(int k = 2; k < nz + 2; k++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				e[s] = (precision) e_s;

				u->ut[s] = ut;
				u->ux[s] = ux;
				u->uy[s] = uy;
				u->un[s] = 0.0;

				up->ut[s] = ut_p;
				up->ux[s] = ux_p;
				up->uy[s] = uy_p;
				up->un[s] = 0.0;
			}
		}
	}
}


// initial conditions for the anisotropic hydro Gubser flow test
// initial Gubser profile different than the "usual" gubser test
void set_aniso_gubser_energy_density_and_flow_profile(int nx, int ny, int nz, double dt, double dx, double dy, double dz, void * hydroParams)
{
	struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;

	double t    = hydro->tau_initial;		// initial longitudinal proper time
	double etas = hydro->shear_viscosity;	// shear viscosity

	double q0 = 1.0;						// inverse length size (hard coded)
	double plpt_ratio = 0.01;				// initial plpt ratio at transverse corner of grid (hard coded)

	// initial T_hat hard coded so that initial central temperature = 0.6 GeV
	//double T0_hat = 0.0455468;   			// plpt_ratio = 1.0
	//double T0_hat = 0.04296357;				// plpt_ratio = 0.01  (x,y = 5fm  x 5fm)
	//double T0_hat = 0.0328418;				// plpt_ratio = 0.01  (x,y = 6fm  x 6fm)
	double T0_hat = 0.0261391;	  			// low plpt ratio (x,y = 7fm x 7fm)
	//double T0_hat = 0.01537397;	  			// plpt_ratio = 0.01  (x,y = 10fm x 10fm)
	//double T0_hat = 0.01171034;	  			// plpt_ratio = 0.01  (x,y = 12fm x 12fm)

	double x_max = 0.5 * (nx - 1.0) * dx;
	double y_max = 0.5 * (ny - 1.0) * dy;

	// distance to transverse corner
	double r_max = sqrt(x_max * x_max  +  y_max * y_max);

    double rho0 = rho_function(t, r_max, q0);	// min rho at transverse corner
	double rhoP = rho_function(t, 0.0, q0);		// max rho at center

	double drho = 0.0001;

	int rho_pts = ceil(fabs((rhoP - rho0) / drho));

	drho = fabs((rhoP - rho0) / ((double)rho_pts - 1.0));

	double rho_array[rho_pts];

	for(int j = 0; j < rho_pts; j++)
	{
		rho_array[j] = rho0  +  j * drho;
	}

	// pass the e, pl rho_arrays in a module function
	double  e_hat[rho_pts];
	double pl_hat[rho_pts];

	 e_hat[0] = EOS_FACTOR * pow(T0_hat, 4);
	pl_hat[0] = e_hat[0] * plpt_ratio / (2.0 + plpt_ratio);


	// make a separate module (move rho_function too)
	gubser_rho_evolution(e_hat, pl_hat, rho_array, rho_pts, drho, etas);

	// construct the cubic spline interpolations
	gsl_spline * e_hat_spline;
	gsl_spline * pl_hat_spline;

	e_hat_spline = gsl_spline_alloc(gsl_interp_cspline, rho_pts);
	pl_hat_spline = gsl_spline_alloc(gsl_interp_cspline, rho_pts);

	gsl_spline_init(e_hat_spline, rho_array, e_hat, rho_pts);
	gsl_spline_init(pl_hat_spline, rho_array, pl_hat, rho_pts);

	gsl_interp_accel * accel = gsl_interp_accel_alloc();

	double eps = 1.e-5;

	// loop over physical grid points
	for(int i = 2; i < nx + 2; i++)
	{
		double x = (i - 2.0 - (nx - 1.0)/2.0) * dx;

		for(int j = 2; j < ny + 2; j++)
		{
			double y = (j - 2.0 - (ny - 1.0)/2.0) * dy;

			double r = sqrt(x * x  +  y * y);

			double rho = rho_function(t, r, q0);

			// if(!(rho >= rho0 && rho < rhoP))
			// {
			// 	printf("Error: rho = %lf out of interval [%lf, %lf)\n", rho, rho0, rhoP);
			// }
			rho = fmax(rho0, fmin(rho, rhoP - eps));

			// evaluate anisotropic profile
			double e_s  = gsl_spline_eval(e_hat_spline,  rho, accel) / (t * t * t * t);
			double pl_s = gsl_spline_eval(pl_hat_spline, rho, accel) / (t * t * t * t);
			double pt_s = 0.5 * (e_s - pl_s);


			double kappa   = atanh(2.0 * q0 * q0 * t * r / (1.0  +  q0 * q0 * (t * t  +  r * r)));
			double kappa_p = atanh(2.0 * q0 * q0 * (t - dt) * r / (1.0  +  q0 * q0 * ((t - dt) * (t - dt)  +  r * r)));

			if(std::isnan(kappa) || std::isnan(kappa_p))
			{
				printf("Gubser initial conditions error: (kappa, kappa_p) = (%lf, %lf)\n", kappa, kappa_p);
				exit(-1);
			}

			// initial flow profile
			double ux = sinh(kappa) * x / r;
			double uy = sinh(kappa) * y / r;

			double ux_p = sinh(kappa_p) * x / r;
			double uy_p = sinh(kappa_p) * y / r;

			if(std::isnan(ux)) ux = 0.0;		// remove 0/0 nan
			if(std::isnan(uy)) uy = 0.0;

			if(std::isnan(ux_p)) ux_p = 0.0;
			if(std::isnan(uy_p)) uy_p = 0.0;

			double ut   = sqrt(1.0  +  ux * ux  +  uy * uy);
			double ut_p = sqrt(1.0  +  ux_p * ux_p  +  uy_p * uy_p);

			for(int k = 2; k < nz + 2; k++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				e[s] = fmax(E_MIN, e_s);

				q->pl[s] = pl_s;

			#if (PT_MATCHING == 1)
				q->pt[s] = pt_s;
			#endif

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

				u->ut[s] = ut;
				u->ux[s] = ux;
				u->uy[s] = uy;
				u->un[s] = 0.0;

				up->ut[s] = ut_p;
				up->ux[s] = ux_p;
				up->uy[s] = uy_p;
				up->un[s] = 0.0;
			}
		}
	}

	gsl_spline_free(e_hat_spline);
	gsl_spline_free(pl_hat_spline);
	gsl_interp_accel_free(accel);
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

	double dt = lattice->latticeSpacingProperTime;
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
			printf("Ideal Gubser\n");
		#ifndef CONFORMAL_EOS
			printf("\nGubser initial condition error: CONFORMAL_EOS not defined in /rhic/include/EquationOfState.h, exiting...\n");
			exit(-1);
		#endif

			set_ideal_gubser_energy_density_and_flow_profile(nx, ny, nz, dt, dx, dy, dz, initCondParams, hydroParams);
			set_equilibrium_initial_condition(nx, ny, nz);
			set_initial_conserved_variables(t, nx, ny, nz);

			break;
		}
		case 3:
		{
			printf("Anisotropic Gubser (residual shear stress initialized to zero)\n");
		#ifndef CONFORMAL_EOS
			printf("\nGubser initial condition error: CONFORMAL_EOS not defined in /rhic/include/EquationOfState.h, exiting...\n");
			exit(-1);
		#endif

			set_aniso_gubser_energy_density_and_flow_profile(nx, ny, nz, dt, dx, dy, dz, hydroParams);
			set_initial_conserved_variables(t, nx, ny, nz);

			break;
		}
		case 4:
		{
			printf("Optical Glauber (fluid velocity and viscous pressures initialized to zero)\n");

			set_Glauber_energy_density_and_flow_profile(nx, ny, nz, dx, dy, dz, initCondParams);
			set_equilibrium_initial_condition(nx, ny, nz);
			set_initial_conserved_variables(t, nx, ny, nz);

			break;
		}

		case 5:
		{
			printf("MC Glauber");
			set_Glauber_energy_density_and_flow_profile(nx, ny, nz, dx, dy, dz, initCondParams);
			printf("(fluid velocity and viscous pressures initialized to zero)\n");
			set_equilibrium_initial_condition(nx, ny, nz);
			set_initial_conserved_variables(t, nx, ny, nz);

			break;
		}
		case 6:
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


