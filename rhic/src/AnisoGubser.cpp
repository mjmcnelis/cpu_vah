#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include "../include/Hydrodynamics.h"
#include "../include/EquationOfState.h"
#include "../include/TransportAniso.h"
#include "../include/DynamicalVariables.h"
#include "../include/Parameters.h"
#include "../include/Macros.h"
#include "../include/AnisoGubser.h"
using namespace std;

inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


inline int sign(double x)			// for the bisection root search
{
	if(x > 0.) return 1;
	else if(x < 0.) return -1;
	else return 0;
}


// rho(t, r, q)
double rho_function(double t, double r, double q)
{
	return - asinh((1.  -  q * q * (t * t  -  r * r)) / (2. * q * t));
}


// energy conservation law (gubser)
double de_drho(double e, double pl, double rho, double t, hydro_parameters hydro)
{
	return (pl  -  3.*e) * tanh(rho);
}


// pl relaxation equation (gubser)
double dpl_drho(double e, double pl, double rho, double t, hydro_parameters hydro)
{
	double conformal_prefactor = hydro.conformal_eos_prefactor;
	double T_hat = pow(e / conformal_prefactor, 0.25);

	double taupiInv = T_hat / (5. * hydro.constant_etas);

	aniso_transport_coefficients aniso;
	aniso.compute_transport_coefficients(e, pl, (e - pl)/2., conformal_prefactor);

	return - taupiInv * (pl - e/3.)  -  (4. * pl  +  aniso.zeta_LL) * tanh(rho);

}


void gubser_rho_evolution(double * e_hat, double * pl_hat, double * rho_array, int rho_pts, double drho, double t, hydro_parameters hydro)
{
	for(int i = 1; i < rho_pts; i++)		// 4th order Runge-Kutta evolution
	{
		double rho = rho_array[i - 1];

		double ep  =  e_hat[i - 1];
		double plp = pl_hat[i - 1];

		double de1  = drho *  de_drho(ep, plp, rho, t, hydro);
		double dpl1 = drho * dpl_drho(ep, plp, rho, t, hydro);

		double de2  = drho *  de_drho(ep + de1/2., plp + dpl1/2., rho + drho/2., t, hydro);
		double dpl2 = drho * dpl_drho(ep + de1/2., plp + dpl1/2., rho + drho/2., t, hydro);

		double de3  = drho *  de_drho(ep + de2/2., plp + dpl2/2., rho + drho/2., t, hydro);
		double dpl3 = drho * dpl_drho(ep + de2/2., plp + dpl2/2., rho + drho/2., t, hydro);

		double de4  = drho *  de_drho(ep + de3, plp + dpl3, rho + drho, t, hydro);
		double dpl4 = drho * dpl_drho(ep + de3, plp + dpl3, rho + drho, t, hydro);

		e_hat[i]  = ep   +  (de1   +  2. * de2   +  2. * de3   +  de4)  / 6.;
		pl_hat[i] = plp  +  (dpl1  +  2. * dpl2  +  2. * dpl3  +  dpl4) / 6.;
	}
}


double compute_initial_central_temperature(double t0, double rho0, double rhoP, int rho_pts, double * rho_array, double drho, double T0_hat, hydro_parameters hydro)
{
	double conformal_eos_prefactor = hydro.conformal_eos_prefactor;
	double plpt_ratio = hydro.plpt_ratio_initial;

	double  e_hat[rho_pts];
	double pl_hat[rho_pts];

	 e_hat[0] = conformal_eos_prefactor * pow(T0_hat, 4);
	pl_hat[0] = e_hat[0] * plpt_ratio / (2. + plpt_ratio);

	gubser_rho_evolution(e_hat, pl_hat, rho_array, rho_pts, drho, t0, hydro);

	// initial central temperature in GeV
	return hbarc / t0 * pow(e_hat[rho_pts - 1] / conformal_eos_prefactor, 0.25);
}


double search_T0_hat(lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro)
{
	double T0 = initial.initialCentralTemperatureGeV;	// solve the root equation: f(T0_hat) = T0

	double t0 = hydro.tau_initial;
	double q0 = initial.q_gubser;

	int nx = lattice.lattice_points_x;
	double dx = lattice.lattice_spacing_x;

	double r_max = (nx - 1.) * dx / sqrt(2.);			// distance from center to transverse corner

	double rho0 = rho_function(t0, r_max, q0);			// start rho ~ transverse corner at t0
	double rhoP = rho_function(t0, 0, q0);				// final rho ~ center at t0

	int rho_pts = ceil(fabs((rhoP - rho0) / drho_default));
	double drho = fabs((rhoP - rho0) / ((double)rho_pts - 1.));
	double rho_array[rho_pts];

	for(int j = 0; j < rho_pts; j++) rho_array[j] = rho0  +  j * drho;

	double T0_hat_1 = 0.0;								// starting T0_hat values for bisection method
	double T0_hat_2 = 1.0;
	double T0_hat_mid = (T0_hat_1 + T0_hat_2) / 2.;

	int n = 0, max_iterations = 50;

	do 	// bisection root search
	{
		// compute the corresponding initial central temperatures in GeV
		double T0_1   = compute_initial_central_temperature(t0, rho0, rhoP, rho_pts, rho_array, drho, T0_hat_1, hydro);
		double T0_2   = compute_initial_central_temperature(t0, rho0, rhoP, rho_pts, rho_array, drho, T0_hat_2, hydro);
		double T0_mid = compute_initial_central_temperature(t0, rho0, rhoP, rho_pts, rho_array, drho, T0_hat_mid, hydro);

		double sign_1 = sign(T0_1 - T0);
		double sign_2 = sign(T0_2 - T0);
		double sign_mid = sign(T0_mid - T0);

		if(sign_1 == sign_2)	// current boundary points should have opposite signs
		{
			printf("T0_initial_hat_search error: sign_1 = %.1f and sign_2 = %.1f have same sign at interation %d\n", sign_1, sign_2, n);
			exit(-1);
		}

		if(sign_mid == sign_1) T0_hat_1 = T0_hat_mid;	// root bracketing
		else T0_hat_2 = T0_hat_mid;

		T0_hat_mid = (T0_hat_1 + T0_hat_2) / 2.;

		n++;

	} while(fabs(T0_hat_2 - T0_hat_1) > T0_hat_error && n < max_iterations);

	if(n == max_iterations)
	{
		printf("T0_initial_hat_search error: root search failed to convergence to desired tolerance\n");
	}

	return T0_hat_mid;
}


double run_semi_analytic_aniso_gubser(lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro)
{
	int nx = lattice.lattice_points_x;
	double dx = lattice.lattice_spacing_x;

	double t0 = hydro.tau_initial;							// initial time
	double dt_out = lattice.output_interval;				// output time intervals

	double q0 = initial.q_gubser;							// inverse length size
	double T0_hat = search_T0_hat(lattice, initial, hydro);	// initial Gubser temperature
	double plpt_ratio = hydro.plpt_ratio_initial;			// initial plpt ratio at transverse corner of grid

	printf("T0_hat = %lf\n\n", T0_hat);

	double r_max = (nx - 1.) * dx / sqrt(2.);				// distance to transverse corner

	double T_freeze = hydro.freezeout_temperature_GeV;
	double e_freeze = hydro.conformal_eos_prefactor * pow(T_freeze / hbarc, 4);

	double t = t0;

	// increase t until grid below freezeout temperature
	while(true)
	{
		double rho0 = rho_function(t0, r_max, q0);			// min rho at transverse corner
		double rhoP = rho_function(t + dt_out, 0, q0);		// max rho at center (this is the problem)

		int rho_pts = ceil(fabs((rhoP - rho0) / drho_default));

		double drho = fabs((rhoP - rho0) / ((double)rho_pts - 1.));

		double rho_array[rho_pts];

		for(int j = 0; j < rho_pts; j++) rho_array[j] = rho0  +  j * drho;

		double  e_hat[rho_pts];
		double pl_hat[rho_pts];

		 e_hat[0] = hydro.conformal_eos_prefactor * pow(T0_hat, 4);
		pl_hat[0] = e_hat[0] * plpt_ratio / (2. + plpt_ratio);

		gubser_rho_evolution(e_hat, pl_hat, rho_array, rho_pts, drho, t, hydro);

		gsl_spline * e_hat_spline;		// construct the cubic spline interpolations
		gsl_spline * pl_hat_spline;
		gsl_interp_accel * accel = gsl_interp_accel_alloc();

		e_hat_spline = gsl_spline_alloc(gsl_interp_cspline, rho_pts);
		pl_hat_spline = gsl_spline_alloc(gsl_interp_cspline, rho_pts);

		gsl_spline_init(e_hat_spline, rho_array, e_hat, rho_pts);
		gsl_spline_init(pl_hat_spline, rho_array, pl_hat, rho_pts);

		FILE *energy;
		FILE *plptratio;
		char fname1[255];
		char fname2[255];

		sprintf(fname1, "semi_analytic/e_gubser_%.3f.dat", t);
		sprintf(fname2, "semi_analytic/plpt_ratio_gubser_%.3f.dat", t);

		energy 		= fopen(fname1, "w");
		plptratio 	= fopen(fname2, "w");

		bool below_freezeout_surface = true;

		for(int i = 2; i < nx + 2; i++)
		{
			double x = (i - 2. - (nx - 1.)/2.) * dx;

			double rho = rho_function(t, fabs(x), q0);

			rho = fmax(rho0, fmin(rho, rhoP - eps));

			double e  = gsl_spline_eval(e_hat_spline,  rho, accel) / (t * t * t * t);
			double pl = gsl_spline_eval(pl_hat_spline, rho, accel) / (t * t * t * t);

			fprintf(energy, 	"%.3f\t%.8f\n", x, e * hbarc);
			fprintf(plptratio,	"%.3f\t%.8f\n", x, 2. * pl / (e - pl));

			if(e > e_freeze) below_freezeout_surface = false;
		}

		fclose(energy);
		fclose(plptratio);

		gsl_interp_accel_free(accel);
		gsl_spline_free(e_hat_spline);
		gsl_spline_free(pl_hat_spline);

		if(below_freezeout_surface) break;

		t += dt_out;
	}

	return T0_hat;
}


void set_aniso_gubser_energy_density_and_flow_profile(double T0_hat, int nx, int ny, int nz, double dt, double dx, double dy, double dz, hydro_parameters hydro, initial_condition_parameters initial)
{
#ifdef ANISO_HYDRO
	double t = hydro.tau_initial;					// initial longitudinal proper time

	double q0 = initial.q_gubser;					// inverse length size
	double plpt_ratio = hydro.plpt_ratio_initial;	// initial plpt ratio at transverse corner of grid

	double x2_max = pow((nx - 1.) * dx / 2., 2);
	double y2_max = pow((ny - 1.) * dy / 2., 2);

	double r_max = sqrt(x2_max + y2_max);			// distance to transverse cornerx

    double rho0 = rho_function(t, r_max, q0);		// min rho at transverse corner
	double rhoP = rho_function(t, 0, q0);			// max rho at center

	int rho_pts = ceil(fabs((rhoP - rho0) / drho_default));

	double drho = fabs((rhoP - rho0) / ((double)rho_pts - 1.));

	double rho_array[rho_pts];

	for(int j = 0; j < rho_pts; j++) rho_array[j] = rho0  +  j * drho;

	double  e_hat[rho_pts];
	double pl_hat[rho_pts];

	 e_hat[0] = hydro.conformal_eos_prefactor * pow(T0_hat, 4);
	pl_hat[0] = e_hat[0] * plpt_ratio / (2. + plpt_ratio);

	gubser_rho_evolution(e_hat, pl_hat, rho_array, rho_pts, drho, t, hydro);

	gsl_spline * e_hat_spline;						// construct the cubic spline interpolations
	gsl_spline * pl_hat_spline;

	 e_hat_spline = gsl_spline_alloc(gsl_interp_cspline, rho_pts);
	pl_hat_spline = gsl_spline_alloc(gsl_interp_cspline, rho_pts);

	gsl_spline_init(e_hat_spline, rho_array, e_hat, rho_pts);
	gsl_spline_init(pl_hat_spline, rho_array, pl_hat, rho_pts);

	gsl_interp_accel * accel = gsl_interp_accel_alloc();

	double eps = 1.e-5;

	// don't use openmp here
	for(int i = 2; i < nx + 2; i++)
	{
		double x = (i - 2. - (nx - 1.)/2.) * dx;

		for(int j = 2; j < ny + 2; j++)
		{
			double y = (j - 2. - (ny - 1.)/2.) * dy;

			double r = sqrt(x * x  +  y * y);

			double rho = rho_function(t, r, q0);

			rho = fmax(rho0, fmin(rho, rhoP - eps));

			// interpolate anisotropic profile
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

			double ux = sinh(kappa) * x / r;
			double uy = sinh(kappa) * y / r;

			double ux_p = sinh(kappa_p) * x / r;
			double uy_p = sinh(kappa_p) * y / r;

			if(std::isnan(ux)) ux = 0;
			if(std::isnan(uy)) uy = 0;

			if(std::isnan(ux_p)) ux_p = 0;
			if(std::isnan(uy_p)) uy_p = 0;

			for(int k = 2; k < nz + 2; k++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				e[s] = energy_density_cutoff(hydro.energy_min, e_s);

				q[s].pl = pl_s;
				q[s].pt = pt_s;

			#ifdef PIMUNU
		  		q[s].pitt = 0;
		  		q[s].pitx = 0;
		  		q[s].pity = 0;
		  		q[s].pixx = 0;
		  		q[s].pixy = 0;
		  		q[s].piyy = 0;
			#endif

				u[s].ux = ux;
				u[s].uy = uy;

				up[s].ux = ux_p;
				up[s].uy = uy_p;

			}
		}
	}

	gsl_spline_free(e_hat_spline);
	gsl_spline_free(pl_hat_spline);
	gsl_interp_accel_free(accel);
#endif
}




















