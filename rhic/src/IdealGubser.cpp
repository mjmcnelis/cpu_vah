#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include "../include/Hydrodynamics.h"
#include "../include/EquationOfState.h"
#include "../include/Parameters.h"
using namespace std;


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



















