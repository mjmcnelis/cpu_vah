#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include "../include/Hydrodynamics.h"
#include "../include/EquationOfState.h"
#include "../include/DynamicalVariables.h"
#include "../include/Parameters.h"
using namespace std;

inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


void set_ideal_gubser_initial_conditions(lattice_parameters lattice, precision dt, initial_condition_parameters initial, hydro_parameters hydro)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	double dx = lattice.lattice_spacing_x;
	double dy = lattice.lattice_spacing_y;

	double t = hydro.tau_initial;								// initial longitudinal proper time
	double t2 = t * t;
	double T0 = initial.initialCentralTemperatureGeV / hbarc;	// central temperature (fm)

	precision conformal_eos_prefactor = hydro.conformal_eos_prefactor;
	precision e_min = hydro.energy_min;

	double q0  = initial.q_gubser;
	double q02 = q0 * q0;
	double q04 = q02 * q02;

	double T0_hat = T0 * t * pow((1. + q02 * t2) / (2. * q0 * t), 2./3.);		// initial central temperature mapped to deSitter space

	printf("T0_hat = %.4f\n\n", T0_hat);

	for(int i = 2; i < nx + 2; i++)
	{
		double x = (i - 2. - (nx - 1.)/2.) * dx;

		for(int j = 2; j < ny + 2; j++)
		{
			double y = (j - 2. - (ny - 1.)/2.) * dy;

			double r = sqrt(x * x  +  y * y);
			double r2 = r * r;

			double kappa   = atanh(2. * q02 * t * r / (1.  +  q02 * (t2  +  r2)));
			double kappa_p = atanh(2. * q02 * (t - dt) * r / (1.  +  q02 * ((t - dt) * (t - dt)  +  r2)));

			if(std::isnan(kappa) || std::isnan(kappa_p))
			{
				printf("set_ideal_gubser_initial_conditions error: (kappa, kappa_p) = (%lf, %lf)\n", kappa, kappa_p);
				exit(-1);
			}

			precision T = (T0_hat / t) * pow(4. * q02 * t2 / (1.  +  2. * q02 * (t2 + r2)  +  q04 * (t2 - r2) * (t2 - r2)), 1./3.);

			precision e_s = equilibrium_energy_density_new(T, conformal_eos_prefactor);

			double ux = sinh(kappa) * x / r;
			double uy = sinh(kappa) * y / r;

			double ux_p = sinh(kappa_p) * x / r;
			double uy_p = sinh(kappa_p) * y / r;

			if(std::isnan(ux)) ux = 0;		// remove 0/0 nans
			if(std::isnan(uy)) uy = 0;

			if(std::isnan(ux_p)) ux_p = 0;
			if(std::isnan(uy_p)) uy_p = 0;

			for(int k = 2; k < nz + 2; k++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				e[s] = energy_density_cutoff(e_min, e_s);

				u[s].ux = ux;
				u[s].uy = uy;

				up[s].ux = ux_p;
				up[s].uy = uy_p;
			}
		}
	}
}


void run_analytic_ideal_gubser(lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro)
{
	int nx = lattice.lattice_points_x;
	double dx = lattice.lattice_spacing_x;

	double t0 = hydro.tau_initial;							// initial time
	double dt_out = lattice.output_interval;				// output time intervals

	double q0 = initial.q_gubser;							// inverse length size
	double q02 = q0 * q0;
	double q04 = q02 * q02;

	double conformal_prefactor = hydro.conformal_eos_prefactor;
	double T_freeze = hydro.freezeout_temperature_GeV / hbarc;

	double t = t0;

	double T_central = initial.initialCentralTemperatureGeV / hbarc;	// central temperature (fm^-1)

	double T_central_hat = T_central * t0 * pow((1. + q02 * t0 * t0) / (2. * q0 * t0), 2./3.);

	while(true)		// increase t until grid below freezeout temperature
	{
		double t2 = t * t;

		FILE *energy;
		char fname1[255];

		sprintf(fname1, "semi_analytic/e_gubser_%.3f.dat", t);

		energy = fopen(fname1, "w");

		bool below_freezeout_surface = true;

		for(int i = 2; i < nx + 2; i++)
		{
			double x = (i - 2. - (nx - 1.)/2.) * dx;
			double r2 = x * x;

			double T = T_central_hat / t * pow(4. * q02 * t2 / (1.  +  2. * q02 * (t2 + r2)  +  q04 * (t2 - r2) * (t2 - r2)), 1./3.);
			double e = conformal_prefactor * T * T * T * T;

			fprintf(energy, "%.3f\t%.8f\n", x, e * hbarc);

			if(T > T_freeze) below_freezeout_surface = false;
		}

		fclose(energy);

		if(below_freezeout_surface) break;

		t += dt_out;
	}
}





















