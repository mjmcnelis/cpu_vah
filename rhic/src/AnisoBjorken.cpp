#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <cmath>
#include "../include/EquationOfState.h"
#include "../include/TransportCoefficients.h"
#include "../include/Hydrodynamics.h"
#include "../include/DynamicalVariables.h"
#include "../include/Macros.h"
#include "../include/Parameters.h"
using namespace std;

const double dt_eps = 1.e-8;

inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


double de_dt(double e, double pl, double t, hydro_parameters hydro)
{
	return - (e + pl) / t;
}


double dpl_dt(double e, double pl, double t, hydro_parameters hydro)
{
	double conformal_prefactor = hydro.conformal_eos_prefactor;

	equation_of_state eos(e);
	double p = eos.equilibrium_pressure();
	double T = eos.effective_temperature(conformal_prefactor);

	double etas = eta_over_s(T, hydro);
	double taupiInv = T / (5. * etas);

	double pt = (e - pl) / 2.; 	// temporary

	transport_coefficients aniso;
	aniso.compute_transport_coefficients(e, pl, pt, conformal_prefactor);


	return - taupiInv * (pl - p)  +  aniso.zeta_LL / t;

}


void run_semi_analytic_aniso_bjorken(lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro)
{
	// using a highly accurate 4th order numerical solution (fixed time step)
	// should be able to handle both conformal and nonconformal Bjorken (keep it at conformal for now)

	double t  = hydro.tau_initial;												// initial time
	double T0  = initial.initialCentralTemperatureGeV;							// initial temperature

	double dt = lattice.min_time_step;											// use minimum time step ~ 100x smaller than t0
	int decimal = - log10(dt);													// setprecision value for t output

	double conformal_eos_prefactor = hydro.conformal_eos_prefactor;
	double e0 = equilibrium_energy_density(T0 / hbarc, conformal_eos_prefactor);	// initial energy density

	double plpt_ratio = hydro.plpt_ratio_initial;								// initial pl/pt ratio

	double e = e0;
	double pl = e0 * plpt_ratio / (2. + plpt_ratio);							// using the conformal formula I think

	double T_freeze = hydro.freezeout_temperature_GeV;
	double e_freeze = equilibrium_energy_density(T_freeze / hbarc, conformal_eos_prefactor);

	ofstream e_e0_plot;
	ofstream pl_pt_plot;
	e_e0_plot.open("semi_analytic/e_e0_aniso_bjorken.dat", ios::out);
	pl_pt_plot.open("semi_analytic/pl_pt_aniso_bjorken.dat", ios::out);

	int steps = 0;

	while(true)
	{
		double pt = (e - pl) / 2.;	// temporary

		if(steps % 1 == 0)
		{
			e_e0_plot  << fixed << setprecision(decimal + 1) << t << "\t" << scientific << setprecision(12) << e / e0 << endl;
			pl_pt_plot << fixed << setprecision(decimal + 1) << t << "\t" << scientific << setprecision(12) << pl / pt << endl;

			if(e < e_freeze) break;
		}

		double ek1  = dt *  de_dt(e, pl, t, hydro);
		double plk1 = dt * dpl_dt(e, pl, t, hydro);

		double ek2  = dt *  de_dt(e + ek1/2., pl + plk1/2., t + dt/2., hydro);
		double plk2 = dt * dpl_dt(e + ek1/2., pl + plk1/2., t + dt/2., hydro);

		double ek3  = dt *  de_dt(e + ek2/2., pl + plk2/2., t + dt/2., hydro);
		double plk3 = dt * dpl_dt(e + ek2/2., pl + plk2/2., t + dt/2., hydro);

		double ek4  = dt *  de_dt(e + ek3, pl + plk3, t + dt, hydro);
		double plk4 = dt * dpl_dt(e + ek3, pl + plk3, t + dt, hydro);

		e  += (ek1   +  2. * ek2   +  2. * ek3   +  ek4) / 6.;
		pl += (plk1  +  2. * plk2  +  2. * plk3  +  plk4) / 6.;

		t += dt;

		steps++;
	}

	e_e0_plot.close();
	pl_pt_plot.close();
}


void set_aniso_bjorken_initial_condition(int nx, int ny, int nz, initial_condition_parameters initial, hydro_parameters hydro)
{
#ifdef ANISO_HYDRO
	precision plpt_ratio = hydro.plpt_ratio_initial;
	precision e_min = hydro.energy_min;
	precision conformal_eos_prefactor = hydro.conformal_eos_prefactor;

	precision T0 = initial.initialCentralTemperatureGeV;							// central temperature (GeV)
	precision e0 = equilibriumEnergyDensity(T0 / hbarc, conformal_eos_prefactor);	// energy density

	for(int i = 2; i < nx + 2; i++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int k = 2; k < nz + 2; k++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				precision e_s = energy_density_cutoff(e_min, e0);

				e[s] = e_s;

				q[s].pl = e_s * plpt_ratio / (2. + plpt_ratio);		// conformal switch approximation

			#if (PT_MATCHING == 1)
				q[s].pt = e_s / (2. + plpt_ratio);
			#endif

			#ifdef PIMUNU
		  		q[s].pitt = 0;		// zero residual shear stress
		  		q[s].pitx = 0;
		  		q[s].pity = 0;
		  		q[s].pixx = 0;
		  		q[s].pixy = 0;
		  		q[s].piyy = 0;
		  	#ifndef BOOST_INVARIANT
		  		q[s].pitn = 0;
		  		q[s].pixn = 0;
		  		q[s].piyn = 0;
		  		q[s].pinn = 0;
		  	#endif

			#endif

		  		u[s].ux = 0;		// zero fluid velocity
				u[s].uy = 0;

				up[s].ux = 0;		// also initialize up = u
				up[s].uy = 0;
			}
		}
	}
#endif
}




