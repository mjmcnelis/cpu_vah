#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <cmath>
#include "../include/EquationOfState.h"
#include "../include/TransportViscous.h"
#include "../include/DynamicalVariables.h"
#include "../include/Viscosities.h"
#include "../include/Macros.h"
#include "../include/Hydrodynamics.h"
#include "../include/Parameters.h"
using namespace std;

inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


double de_dt_vh(double e, double pi, double bulk, double t, hydro_parameters hydro)
{
	equation_of_state eos(e);
	double p = eos.equilibrium_pressure();

	return - (e + p + bulk - pi) / t;
}


double dpi_dt(double e, double pi, double bulk, double t, hydro_parameters hydro)
{
#ifndef ANISO_HYDRO
#ifdef PIMUNU
	equation_of_state eos(e);
	double p = eos.equilibrium_pressure();
	double T = eos.effective_temperature(hydro.conformal_eos_prefactor);
	double s = (e + p) / T;

	double etas = eta_over_s(T, hydro);
	double eta = s * etas;

	viscous_transport_coefficients viscous(T, e, p, hydro.kinetic_theory_model);

	viscous.compute_shear_transport_coefficients(etas);

	double taupi_inverse = viscous.taupi_inverse;
	double deltapipi = viscous.delta_pipi;
	double taupipi = viscous.tau_pipi;
#ifdef PI
	double lambdapibulk = viscous.lambda_pibulkPi;
#else
	double lambdapibulk = 0;
#endif

	return taupi_inverse * (-pi  +  4./3. * eta / t)  -  (taupipi/3. + deltapipi) * pi / t  +  2./3. * lambdapibulk * bulk / t;
#else
	return 0;
#endif
#else
	return 0;
#endif
}

double dbulk_dt(double e, double pi, double bulk, double t, hydro_parameters hydro)
{
#ifndef ANISO_HYDRO
#ifdef PI
	equation_of_state eos(e);
	double p = eos.equilibrium_pressure();
	double T = eos.effective_temperature(hydro.conformal_eos_prefactor);
	double s = (e + p) / T;
	double cs2 = eos.speed_of_sound_squared();

	double zetas = zeta_over_s(T, hydro);
	double zeta = s * zetas;

	viscous_transport_coefficients viscous(T, e, p, hydro.kinetic_theory_model);

	viscous.compute_bulk_transport_coefficients(zetas, 1./3. - cs2);

	double taubulk_inverse = viscous.taubulk_inverse;
	double deltabulkbulk = viscous.delta_bulkPibulkPi;

#ifdef PIMUNU
	double lambdabulkpi = viscous.lambda_bulkPipi;
#else
	double lambdabulkpi = 0;
#endif

	return - taubulk_inverse * (bulk  +  zeta / t)  +  (- deltabulkbulk * bulk  +  lambdabulkpi * pi) / t;
#else
	return 0;
#endif
#else
	return 0;
#endif
}


void run_semi_analytic_viscous_bjorken(lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro)
{
	// using a highly accurate 4th order numerical solution (fixed time step)
	// should be able to handle both conformal and nonconformal Bjorken (keep it at conformal for now)

	double t  = hydro.tau_initial;												// initial time
	double T0  = initial.initialCentralTemperatureGeV;							// initial temperature

	double dt = lattice.min_time_step;											// use minimum time step ~ 100x smaller than t0
	int decimal = - log10(dt);													// setprecision value for t output

	double conformal_eos_prefactor = hydro.conformal_eos_prefactor;
	double e0 = equilibrium_energy_density(T0/hbarc, conformal_eos_prefactor);	// initial energy density

	double plpt_ratio = hydro.plpt_ratio_initial;								// initial pl/pt ratio

	double e = e0;

	equation_of_state eos(e);
	double p = eos.equilibrium_pressure();
	double pl = e * plpt_ratio / (2. + plpt_ratio);
	double pt = (e - pl) / 2.;

	double pi = 0;
	double bulk = 0;

#ifdef PIMUNU
	pi = - 2. * (pl - pt) / 3.;													// - t2.pinn
#endif
#ifdef PI
	bulk = e/3. - p;
#endif

	double T_freeze = hydro.freezeout_temperature_GeV;
	double e_freeze = equilibrium_energy_density(T_freeze / hbarc, conformal_eos_prefactor);

	ofstream e_e0_plot;
	ofstream pl_pt_plot;
	ofstream Rpi_inv_plot;
	ofstream Rbulk_inv_plot;
	ofstream T_plot; 
	e_e0_plot.open("semi_analytic/e_e0_viscous_bjorken.dat", ios::out);
	pl_pt_plot.open("semi_analytic/pl_pt_viscous_bjorken.dat", ios::out);
	Rpi_inv_plot.open("semi_analytic/Rpi_inv_viscous_bjorken.dat", ios::out);
	Rbulk_inv_plot.open("semi_analytic/Rbulk_inv_viscous_bjorken.dat", ios::out);
	T_plot.open("semi_analytic/T_viscous_bjorken.dat", ios::out);

	//dt = 0.001;

	while(true)
	{
		//p = equilibriumPressure(e);
		equation_of_state EoS(e);
		p = EoS.equilibrium_pressure();
		
		double T = EoS.effective_temperature(hydro.conformal_eos_prefactor);

		pl = p  +  bulk  -  pi;
		pt = p  +  bulk  +  pi/2.;

		e_e0_plot      << fixed << setprecision(decimal + 1) << t << "\t" << scientific << setprecision(12) << e / e0                  << endl;
		pl_pt_plot     << fixed << setprecision(decimal + 1) << t << "\t" << scientific << setprecision(12) << pl / pt                 << endl;
		Rpi_inv_plot   << fixed << setprecision(decimal + 1) << t << "\t" << scientific << setprecision(12) << fabs(pi / sqrt(2.) / p) << endl;
		Rbulk_inv_plot << fixed << setprecision(decimal + 1) << t << "\t" << scientific << setprecision(12) << fabs(bulk / p)          << endl;
		T_plot 		   << fixed << setprecision(decimal + 1) << t << "\t" << scientific << setprecision(12) << T        			   << endl;

		if(e < e_freeze) break;

		double e1    = dt * de_dt_vh(e, pi, bulk, t, hydro);
		double pi1   = dt *   dpi_dt(e, pi, bulk, t, hydro);
		double bulk1 = dt * dbulk_dt(e, pi, bulk, t, hydro);

		double e2    = dt * de_dt_vh(e + e1/2., pi + pi1/2., bulk + bulk1/2., t + dt/2., hydro);
		double pi2   = dt *   dpi_dt(e + e1/2., pi + pi1/2., bulk + bulk1/2., t + dt/2., hydro);
		double bulk2 = dt * dbulk_dt(e + e1/2., pi + pi1/2., bulk + bulk1/2., t + dt/2., hydro);

		double e3    = dt * de_dt_vh(e + e2/2., pi + pi2/2., bulk + bulk2/2., t + dt/2., hydro);
		double pi3   = dt *   dpi_dt(e + e2/2., pi + pi2/2., bulk + bulk2/2., t + dt/2., hydro);
		double bulk3 = dt * dbulk_dt(e + e2/2., pi + pi2/2., bulk + bulk2/2., t + dt/2., hydro);

		double e4    = dt * de_dt_vh(e + e3, pi + pi3, bulk + bulk3/2., t + dt, hydro);
		double pi4   = dt *   dpi_dt(e + e3, pi + pi3, bulk + bulk3/2., t + dt, hydro);
		double bulk4 = dt * dbulk_dt(e + e3, pi + pi3, bulk + bulk3/2., t + dt, hydro);

		e    += (e1     +  2. * e2     +  2. * e3     +  e4)    / 6.;
		pi   += (pi1    +  2. * pi2    +  2. * pi3    +  pi4)   / 6.;
		bulk += (bulk1  +  2. * bulk2  +  2. * bulk3  +  bulk4) / 6.;

		t += dt;
	}

	e_e0_plot.close();
	pl_pt_plot.close();
	Rpi_inv_plot.close();
	Rbulk_inv_plot.close();
	T_plot.close();
}


void set_viscous_bjorken_initial_condition(int nx, int ny, int nz, initial_condition_parameters initial, hydro_parameters hydro)
{
#ifndef ANISO_HYDRO
	precision t = hydro.tau_initial;
	precision t2 = t * t;
	precision plpt_ratio = hydro.plpt_ratio_initial;
	precision e_min = hydro.energy_min;
	precision conformal_eos_prefactor = hydro.conformal_eos_prefactor;

	precision T0 = initial.initialCentralTemperatureGeV;							// central temperature (GeV)
	precision e0 = equilibrium_energy_density(T0 / hbarc, conformal_eos_prefactor);	// energy density

	for(int i = 2; i < nx + 2; i++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int k = 2; k < nz + 2; k++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				precision e_s = energy_density_cutoff(e_min, e0);

				e[s] = e_s;

				equation_of_state eos(e_s);
				precision p = eos.equilibrium_pressure();
				precision pl = e_s * plpt_ratio / (2. + plpt_ratio);		// conformal switch approximation
				precision pt = e_s / (2. + plpt_ratio);

			#ifdef PIMUNU
		  		q[s].pitt = 0;
		  		q[s].pitx = 0;
		  		q[s].pity = 0;
		  		q[s].pixx = - (pl - pt) / 3.;
		  		q[s].pixy = 0;
		  		q[s].piyy = - (pl - pt) / 3.;
		  		q[s].pinn = 2. * (pl - pt) / (3. * t2);
			#endif

		  	#ifdef PI
		  		q[s].Pi = (2. * pt  +  pl) / 3.  -  p;
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


