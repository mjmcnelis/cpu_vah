#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include "../include/Macros.h"
#include "../include/Hydrodynamics.h"
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


inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


void set_initial_T_taumu_variables(double t, int nx, int ny, int nz)
{
	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				precision e_s = e[s];

				precision ux = u[s].ux;
				precision uy = u[s].uy;
				precision un = u[s].un;
				precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t * t * un * un);

			#ifdef ANISO_HYDRO
				precision utperp = sqrt(1.  +  ux * ux  +  uy * uy);
				precision zt = t * un / utperp;
				precision zn = ut / t / utperp;

				precision pl = q[s].pl;

			#if (PT_MATCHING == 1)
				precision pt = q[s].pt;
			#else
				precision pt = (e_s - pl) / 2.;
			#endif

			#else
				precision p = equilibriumPressure(e_s);
			#endif

			#ifdef PIMUNU
				precision pitt = q[s].pitt;
				precision pitx = q[s].pitx;
				precision pity = q[s].pity;
				precision pitn = q[s].pitn;
			#else
				precision pitt = 0;
				precision pitx = 0;
				precision pity = 0;
				precision pitn = 0;
			#endif

			#ifdef WTZMU
				precision WtTz = q[s].WtTz;
				precision WxTz = q[s].WxTz;
				precision WyTz = q[s].WyTz;
				precision WnTz = q[s].WnTz;
			#else
				precision WtTz = 0;
				precision WxTz = 0;
				precision WyTz = 0;
				precision WnTz = 0;
			#endif

			#ifdef PI
				precision Pi = q[s].Pi;
			#else
				precision Pi = 0;
			#endif

			#ifdef ANISO_HYDRO
				q[s].ttt = (e_s + pt) * ut * ut  -   pt  +  (pl - pt) * zt * zt  +  2. * WtTz * zt  +  pitt;
				q[s].ttx = (e_s + pt) * ut * ux  +  WxTz * zt  +  pitx;
				q[s].tty = (e_s + pt) * ut * uy  +  WyTz * zt  +  pity;
				q[s].ttn = (e_s + pt) * ut * un  +  (pl - pt) * zt * zn  +  WtTz * zn  +  WnTz * zt  +  pitn;
			#else
				q[s].ttt = (e_s + p + Pi) * ut * ut  -   p  -  Pi  +  pitt;
				q[s].ttx = (e_s + p + Pi) * ut * ux  +  pitx;
				q[s].tty = (e_s + p + Pi) * ut * uy  +  pity;
				q[s].ttn = (e_s + p + Pi) * ut * un  +  pitn;
			#endif
			}
		}
	}
}


void set_anisotropic_initial_condition(int nx, int ny, int nz, hydro_parameters hydro)
{
	precision plpt_ratio = hydro.plpt_ratio_initial;

	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				precision e_s = e[s];

			#ifdef ANISO_HYDRO
				q[s].pl = e_s * plpt_ratio / (2. + plpt_ratio);		// conformal approximation

			#if (PT_MATCHING == 1)
				q[s].pt = e_s / (2. + plpt_ratio);
			#endif

			#endif

			#ifdef PIMUNU
		  		q[s].pitt = 0;
		  		q[s].pitx = 0;
		  		q[s].pity = 0;
		  		q[s].pitn = 0;
		  		q[s].pixx = 0;
		  		q[s].pixy = 0;
		  		q[s].pixn = 0;
		  		q[s].piyy = 0;
		  		q[s].piyn = 0;
		  		q[s].pinn = 0;
			#endif

			#ifdef WTZMU
		  		q[s].WtTz = 0;
		  		q[s].WxTz = 0;
		  		q[s].WyTz = 0;
		  		q[s].WnTz = 0;
			#endif

		  	#ifdef PI
		  		q[s].Pi = e_s / 3.  -  equilibriumPressure(e_s);	// switching from conformal eos to lattice
		  	#endif
			}
		}
	}
}


// constant initial energy density profile (Bjorken)
void set_Bjorken_energy_density_and_flow_profile(int nx, int ny, int nz, initial_condition_parameters initial, hydro_parameters hydro)
{
	precision conformal_eos_prefactor = hydro.conformal_eos_prefactor;

	precision T0 = initial.initialCentralTemperatureGeV;							// central temperature (GeV)
	precision e0 = equilibriumEnergyDensity(T0 / hbarc, conformal_eos_prefactor);	// energy density

	precision e_min = hydro.energy_min;

	for(int i = 2; i < nx + 2; i++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int k = 2; k < nz + 2; k++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				e[s] = energy_density_cutoff(e_min, e0);

				u[s].ux = 0;	// zero fluid velocity
				u[s].uy = 0;
				u[s].un = 0;

				up[s].ux = 0;	// also initialize up = u
				up[s].uy = 0;
				up[s].un = 0;
			}
		}
	}
}


// Longitudinal energy density profile function eL(eta)
void longitudinal_Energy_Density_Profile(double * const __restrict__ eL, int nz, double dz, initial_condition_parameters initial)
{
	double etaFlat = initial.rapidityMean;
	double etaVariance = initial.rapidityVariance;

	// profile along eta direction is a smooth plateu that exponentially decays when |eta| > etaFlat
	for(int k = 0; k < nz; k++)
	{
		double eta = (k - (nz - 1.)/2.) * dz;				// physical eta points
		double etaScaled = fabs(eta)  -  0.5 * etaFlat;

		eL[k] = exp(- 0.5 * etaScaled * etaScaled / etaVariance * THETA_FUNCTION(etaScaled));
	}
}



// Optical or MC glauber initial energy density profile
void set_Glauber_energy_density_and_flow_profile(int nx, int ny, int nz, double dx, double dy, double dz, initial_condition_parameters initial, hydro_parameters hydro)
{
	int initialConditionType = initial.initialConditionType;

	double T0 = initial.initialCentralTemperatureGeV;									// central temperature (GeV)
	double e0 = equilibriumEnergyDensity(T0 / hbarc, hydro.conformal_eos_prefactor);	// energy density scale factor

	precision e_min = hydro.energy_min;

	double eT[nx * ny];		// normalized transverse profile
	double eL[nz];			// normalized longitudinal profile

	// compute normalized transverse and longitudinal profiles
	if(initialConditionType == 4)
	{
		Optical_Glauber_energy_density_transverse_profile(eT, nx, ny, dx, dy, initial);
	}
	else if(initialConditionType == 5)
	{
		MC_Glauber_energy_density_transverse_profile(eT, nx, ny, dx, dy, initial);
	}
	longitudinal_Energy_Density_Profile(eL, nz, dz, initial);

	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				precision e_s = e0 * eT[i - 2 + (j - 2) * nx] * eL[k - 2];

				e[s] = energy_density_cutoff(e_min, e_s);

				u[s].ux = 0.0;		// zero initial velocity
				u[s].uy = 0.0;
				u[s].un = 0.0;

				up[s].ux = 0.0;		// also set up = u
				up[s].uy = 0.0;
				up[s].un = 0.0;
			}
		}
	}
}



// initial conditions for the ideal Gubser flow test
void set_ideal_gubser_energy_density_and_flow_profile(int nx, int ny, int nz, double dt, double dx, double dy, double dz, initial_condition_parameters initial, hydro_parameters hydro)
{
	double t = hydro.tau_initial;								// initial longitudinal proper time
	double T0 = initial.initialCentralTemperatureGeV / hbarc;	// central temperature (fm)

	precision conformal_eos_prefactor = hydro.conformal_eos_prefactor;
	precision e_min = hydro.energy_min;

	double q  = 1.;												// inverse length size (hard coded)
	double q2 = q * q;
	double q4 = q2 * q2;
	double t2 = t * t;

	// normalize Gubser ideal temperature profile so that central temperature = T0
	double T0_hat = T0 * t * pow((1. + q2 * t2) / (2. * q * t), 2./3.);

	// loop over physical grid points
	for(int i = 2; i < nx + 2; i++)
	{
		double x = (i - 2. - (nx - 1.)/2.) * dx;

		for(int j = 2; j < ny + 2; j++)
		{
			double y = (j - 2. - (ny - 1.)/2.) * dy;

			double r = sqrt(x * x  +  y * y);
			double r2 = r * r;

			double kappa   = atanh(2. * q2 * t * r / (1.  +  q2 * (t2  +  r2)));
			double kappa_p = atanh(2. * q2 * (t - dt) * r / (1.  +  q2 * ((t - dt) * (t - dt)  +  r2)));

			if(std::isnan(kappa) || std::isnan(kappa_p))
			{
				printf("Gubser initial conditions error: (kappa, kappa_p) = (%lf, %lf)\n", kappa, kappa_p);
				exit(-1);
			}

			// temperature profile
			precision T = (T0_hat / t) * pow(2. * q * t, 2./3.) / pow(1.  +  2. * q2 * (t2 + r2)  +  q4 * (t2 - r2) * (t2 - r2), 1./3.);

			precision e_s = equilibriumEnergyDensity(T, conformal_eos_prefactor);

			double ux = sinh(kappa) * x / r;
			double uy = sinh(kappa) * y / r;

			double ux_p = sinh(kappa_p) * x / r;
			double uy_p = sinh(kappa_p) * y / r;

			if(std::isnan(ux)) ux = 0;		// remove x/r = 0/0 nan
			if(std::isnan(uy)) uy = 0;		// remove y/r = 0/0 nan

			if(std::isnan(ux_p)) ux_p = 0;	// remove x/r = 0/0 nan
			if(std::isnan(uy_p)) uy_p = 0;	// remove y/r = 0/0 nan

			for(int k = 2; k < nz + 2; k++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				e[s] = energy_density_cutoff(e_min, e_s);

				u[s].ux = ux;
				u[s].uy = uy;
				u[s].un = 0;

				up[s].ux = ux_p;
				up[s].uy = uy_p;
				up[s].un = 0;
			}
		}
	}
}


// initial conditions for the anisotropic hydro Gubser flow test
// initial Gubser profile different than the "usual" gubser test
void set_aniso_gubser_energy_density_and_flow_profile(int nx, int ny, int nz, double dt, double dx, double dy, double dz, hydro_parameters hydro)
{
#ifdef ANISO_HYDRO
	double t    = hydro.tau_initial;		// initial longitudinal proper time


	double etas = hydro.constant_etas;		// shear viscosity  (NOTE: needs update!)



	precision e_min = hydro.energy_min;



	double q0 = 1.0;						// inverse length size (hard coded)
	double plpt_ratio = 0.01;				// initial plpt ratio at transverse corner of grid (hard coded)

	// initial T_hat hard coded so that initial central temperature = 0.6 GeV
	//double T0_hat = 0.0455468;   			// plpt_ratio = 1.0
	//double T0_hat = 0.04296357;			// plpt_ratio = 0.01  (x,y = 5fm  x 5fm)
	//double T0_hat = 0.0328418;			// plpt_ratio = 0.01  (x,y = 6fm  x 6fm)
	double T0_hat = 0.0261391;	  			// low plpt ratio (x,y = 7fm x 7fm)
	//double T0_hat = 0.01537397;	  		// plpt_ratio = 0.01  (x,y = 10fm x 10fm)
	//double T0_hat = 0.01171034;	  		// plpt_ratio = 0.01  (x,y = 12fm x 12fm)

	double x_max = (nx - 1.) * dx / 2.;
	double y_max = (ny - 1.) * dy / 2.;

	// distance to transverse corner
	double r_max = sqrt(x_max * x_max  +  y_max * y_max);

    double rho0 = rho_function(t, r_max, q0);	// min rho at transverse corner
	double rhoP = rho_function(t, 0, q0);		// max rho at center

	double drho = 0.0001;

	int rho_pts = ceil(fabs((rhoP - rho0) / drho));

	drho = fabs((rhoP - rho0) / ((double)rho_pts - 1.));

	double rho_array[rho_pts];

	for(int j = 0; j < rho_pts; j++)
	{
		rho_array[j] = rho0  +  j * drho;
	}

	// pass the e, pl rho_arrays in a module function
	double  e_hat[rho_pts];
	double pl_hat[rho_pts];

	 e_hat[0] = equilibriumEnergyDensity(T0_hat, hydro.conformal_eos_prefactor);
	pl_hat[0] = e_hat[0] * plpt_ratio / (2. + plpt_ratio);


	// make a separate module (move rho_function too)
	gubser_rho_evolution(e_hat, pl_hat, rho_array, rho_pts, drho, hydro);

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
		double x = (i - 2. - (nx - 1.)/2.) * dx;

		for(int j = 2; j < ny + 2; j++)
		{
			double y = (j - 2. - (ny - 1.)/2.) * dy;

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
			double pt_s = (e_s - pl_s) / 2.;


			double kappa   = atanh(2. * q0 * q0 * t * r / (1.  +  q0 * q0 * (t * t  +  r * r)));
			double kappa_p = atanh(2. * q0 * q0 * (t - dt) * r / (1.  +  q0 * q0 * ((t - dt) * (t - dt)  +  r * r)));

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

			if(std::isnan(ux)) ux = 0;		// remove 0/0 nan
			if(std::isnan(uy)) uy = 0;

			if(std::isnan(ux_p)) ux_p = 0;
			if(std::isnan(uy_p)) uy_p = 0;


			for(int k = 2; k < nz + 2; k++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				e[s] = energy_density_cutoff(e_min, e_s);

				q[s].pl = pl_s;
			#if (PT_MATCHING == 1)
				q[s].pt = pt_s;
			#endif
			#ifdef PIMUNU
		  		q[s].pitt = 0;
		  		q[s].pitx = 0;
		  		q[s].pity = 0;
		  		q[s].pitn = 0;
		  		q[s].pixx = 0;
		  		q[s].pixy = 0;
		  		q[s].pixn = 0;
		  		q[s].piyy = 0;
		  		q[s].piyn = 0;
		  		q[s].pinn = 0;
			#endif
			#ifdef WTZMU
		  		q[s].WtTz = 0;
		  		q[s].WxTz = 0;
		  		q[s].WyTz = 0;
		  		q[s].WnTz = 0;
			#endif

				u[s].ux = ux;
				u[s].uy = uy;
				u[s].un = 0;

				up[s].ux = ux_p;
				up[s].uy = uy_p;
				up[s].un = 0;
			}
		}
	}

	gsl_spline_free(e_hat_spline);
	gsl_spline_free(pl_hat_spline);
	gsl_interp_accel_free(accel);
#endif
}


// Initial conditions for the T^{\tau\mu}, energy density, equilibrium pressure, fluid velocity and viscous pressures
//////////////////////////////////////
//	1 - Bjorken
//	2 - Ideal Gubser
//	3 - Aniso Gubser
//	4 - Optical Glauber
//	5 - MC Glauber
//////////////////////////////////////
void set_initial_conditions(precision t, lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_x;
	int nz = lattice.lattice_points_eta;

	double dx = lattice.lattice_spacing_x;
	double dy = lattice.lattice_spacing_y;
	double dz = lattice.lattice_spacing_eta;

	precision dt = lattice.fixed_time_step;

	if(lattice.adaptive_time_step) dt = lattice.min_time_step;

	printf("\nInitial conditions = ");
	switch(initial.initialConditionType)
	{
		case 1:
		{
			printf("Bjorken ");
		#ifndef CONFORMAL_EOS
			printf("QCD ");
		#else
			printf("Conformal ");
		#endif
			printf("(fluid velocity and viscous pressures initialized to zero)\n\n");

			set_Bjorken_energy_density_and_flow_profile(nx, ny, nz, initial, hydro);
			set_anisotropic_initial_condition(nx, ny, nz, hydro);
			set_initial_T_taumu_variables(t, nx, ny, nz);

			break;
		}
		case 2:
		{
			printf("Ideal Gubser (temporarily out of commission)\n\n");		// should get rid of this (or change it to a viscous Gubser)
			exit(-1);
		// #ifndef CONFORMAL_EOS
		// 	printf("\nGubser initial condition error: CONFORMAL_EOS not defined in /rhic/include/EquationOfState.h, exiting...\n");
		// 	exit(-1);
		// #endif

		// 	set_ideal_gubser_energy_density_and_flow_profile(nx, ny, nz, dt, dx, dy, dz, initial, hydro);
		// 	set_anisotropic_initial_condition(nx, ny, nz);
		// 	set_initial_T_taumu_variables(t, nx, ny, nz);

			break;
		}
		case 3:
		{
		#ifndef ANISO_HYDRO
			printf("Aniso Gubser error: ANISO_HYDRO not defined in /rhic/include/DynamicalVariables.h, exiting...\n");
			exit(-1);
		#endif
			printf("Anisotropic Gubser (residual shear stress initialized to zero)\n\n");
		#ifndef CONFORMAL_EOS
			printf("\nGubser initial condition error: CONFORMAL_EOS not defined in /rhic/include/EquationOfState.h, exiting...\n");
			exit(-1);
		#endif

			set_aniso_gubser_energy_density_and_flow_profile(nx, ny, nz, dt, dx, dy, dz, hydro);
			set_initial_T_taumu_variables(t, nx, ny, nz);

			break;
		}
		case 4:
		{
			printf("Optical Glauber (fluid velocity and viscous pressures initialized to zero)\n\n");

			set_Glauber_energy_density_and_flow_profile(nx, ny, nz, dx, dy, dz, initial, hydro);
			set_anisotropic_initial_condition(nx, ny, nz, hydro);
			set_initial_T_taumu_variables(t, nx, ny, nz);

			break;
		}

		case 5:
		{
			printf("MC Glauber");
			set_Glauber_energy_density_and_flow_profile(nx, ny, nz, dx, dy, dz, initial, hydro);
			printf("(fluid velocity and viscous pressures initialized to zero)\n\n");
			set_anisotropic_initial_condition(nx, ny, nz, hydro);
			set_initial_T_taumu_variables(t, nx, ny, nz);

			break;
		}
		case 6:
		{
			printf("Trento + F.S.");
			printf("\t(nothing here yet! Exiting...)\n");
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


