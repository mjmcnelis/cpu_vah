#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <cmath>
#include "../include/EquationOfState.h"
#include "../include/TransportAniso.h"
#include "../include/TransportViscous.h"
#include "../include/Viscosities.h"
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


double de_dt_conformal(double e, double pl, double t, hydro_parameters hydro)
{
	return - (e + pl) / t;
}

double de_dt(double e, double pl, double pt, double t, hydro_parameters hydro)
{
	return - (e + pl) / t;
}


double dpl_dt_conformal(double e, double pl, double t, hydro_parameters hydro)
{
	double conformal_prefactor = hydro.conformal_eos_prefactor;

	equation_of_state eos(e);
	double p = eos.equilibrium_pressure();
	double T = eos.effective_temperature(conformal_prefactor);

	double etas = eta_over_s(T, hydro);
	double taupiInv = T / (5. * etas);

	double pt = (e - pl) / 2.; 	// temporary

	aniso_transport_coefficients aniso;
	aniso.compute_transport_coefficients(e, pl, pt, conformal_prefactor);

	return - taupiInv * (pl - p)  +  aniso.zeta_LL / t;

}


void run_semi_analytic_aniso_bjorken(lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro)
{
	// using a highly accurate 4th order numerical solution (fixed time step)
	// kinetic theory model should be quasiparticle but it can run either lattice + quasiparticle or conformal (or switch_eos mode)

	if(hydro.kinetic_theory_model != 1)
	{
		printf("run_semi_analytic_aniso_bjorken error: set kinetic_theory_model = 1 (default setting is quasiparticle)\n");
		exit(-1);
	}

	double t  = hydro.tau_initial;												// initial time
	double dt = lattice.min_time_step;											// use minimum time step ~ 100x smaller than t0
	int decimal = - log10(dt);													// setprecision value for t output


	// initial state information
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::
	double T0  = initial.initialCentralTemperatureGeV;							// initial temperature
	double T = T0 / hbarc;
	double plpt_ratio = hydro.plpt_ratio_initial;								// initial pl/pt ratio
	double conformal_prefactor = hydro.conformal_eos_prefactor;
	double e0 = equilibrium_energy_density(T, conformal_prefactor);				// initial energy density

	equation_of_state eos(e0);
	double p0 = eos.equilibrium_pressure();										// initial thermal pressure
	double s0 = (e0 + p0) / T;													// initial entropy density
	double cs2 = eos.speed_of_sound_squared();
	double B0 = eos.equilibrium_mean_field(T);									// initial thermal mean field
	double m = T * eos.z_quasi(T);												// quasiparticle mass
	double m2 = m * m;
	double mdmde = eos.mdmde_quasi();

	double zetas = zeta_over_s(T, hydro);										// compute bulk relaxation time for asymptotic dB
	double betabulk = eos.beta_bulk(T);
	double taubulk = zetas * s0 / betabulk;
	
	printf("initial energy density = %lf fm^-4\n", e0);
	printf("initial thermal pressure = %lf fm^-4\n", p0);
	printf("initial entropy density = %lf fm^-4\n", s0);
	printf("initial thermal mean field = %lf fm^-4\n", B0);
	printf("quasiparticle mass = %lf\n", m);
	printf("quasiparticle mdmde = %lf\n", mdmde);
	printf("zeta/S = %lf\n", zetas);
	printf("betabulk = %lf fm^-4\n", betabulk);
	printf("taubulk = %lf fm\n\n", taubulk);
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::

	
	// initialize hydrodynamic variables (e, pl, pt, B, lambda, aL, aT)
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::
	
	double e = e0;
	double pl, pt;

#ifdef CONFORMAL_EOS
	pl = e0 * plpt_ratio / (2. + plpt_ratio);									// using the conformal formula I think
	pt = e0 / (2. + plpt_ratio);
#endif

#ifdef LATTICE_QCD
	pl = p0;																	// start with equilibrium initial conditions for now
	pt = p0;
#endif

	// asymptotic formula for non-equilibrium mean field component dB
	double dB0 = -3. * taubulk * mdmde * (e + pl) * (2.*pt/3. + pl/3. - p0) / (t * m2) / (1.  +  4. * taubulk * mdmde * (e + pl) / (t * m2));														// keep dB = 0 for now
																				// simulation needs different initialization formula than Bjorken

	double B = B0 + dB0;														// initial mean field
	double lambda = T;															// equilibrium anisotropic variables for now 
	double aL = 1.;
	double aT = 1.;

	printf("initial pl = %lf fm^-4\n", pl);
	printf("initial pt = %lf fm^-4\n", pt);
	printf("initial dB = %lf fm^-4\n", dB0);
	printf("initial B  = %lf fm^-4\n", B);
	printf("initial lambda = %lf fm^-1\n", lambda);
	exit(-1);

	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::


	

	// freezeout condition
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::
	double T_freeze = hydro.freezeout_temperature_GeV;
	double e_freeze = equilibrium_energy_density(T_freeze / hbarc, conformal_prefactor);
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::



	// plot observables
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::
	ofstream e_e0_plot;
	ofstream pl_pt_plot;
	ofstream B_plot;
	e_e0_plot.open("semi_analytic/e_e0_aniso_bjorken.dat", ios::out);
	pl_pt_plot.open("semi_analytic/pl_pt_aniso_bjorken.dat", ios::out);
	B_plot.open("semi_analytic/B_aniso_bjorken.dat", ios::out);
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::


	int steps = 0;

	// start evolution
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::
	while(true)
	{
		double pt = (e - pl) / 2.;	// temporary

		if(steps % 1 == 0)			// write to file 
		{
			e_e0_plot  << fixed << setprecision(decimal + 1) << t << "\t" << scientific << setprecision(12) << e / e0 << endl;
			pl_pt_plot << fixed << setprecision(decimal + 1) << t << "\t" << scientific << setprecision(12) << pl / pt << endl;

			if(e < e_freeze) break;
		}

	// I need to think about how to use these if statements correctly
	#ifdef CONFORMAL_EOS
		double ek1  = dt *  de_dt_conformal(e, pl, t, hydro);
		double plk1 = dt * dpl_dt_conformal(e, pl, t, hydro);

		double ek2  = dt *  de_dt_conformal(e + ek1/2., pl + plk1/2., t + dt/2., hydro);
		double plk2 = dt * dpl_dt_conformal(e + ek1/2., pl + plk1/2., t + dt/2., hydro);

		double ek3  = dt *  de_dt_conformal(e + ek2/2., pl + plk2/2., t + dt/2., hydro);
		double plk3 = dt * dpl_dt_conformal(e + ek2/2., pl + plk2/2., t + dt/2., hydro);

		double ek4  = dt *  de_dt_conformal(e + ek3, pl + plk3, t + dt, hydro);
		double plk4 = dt * dpl_dt_conformal(e + ek3, pl + plk3, t + dt, hydro);

		e  += (ek1   +  2. * ek2   +  2. * ek3   +  ek4) / 6.;
		pl += (plk1  +  2. * plk2  +  2. * plk3  +  plk4) / 6.;
	#endif

	#ifdef LATTICE_QCD

		// let's use the RK2 method for now since the old code uses that
		double de1  = dt *  de_dt(e, pl, pt, t, hydro);
		double dpl1 = dt * dpl_dt(e, pl, pt, t, hydro);
		double dpt1 = dt * dpt_dt(e, pl, pt, t, hydro);

		double e_mid = e + de1;
		double pl_mid = pl + dpl1;
		double pt_mid = pt + dpt1;

		// compute anisotropic variables from intermediate values (e_mid, pl_mid, pt_mid)

		double de2 = dt *  de_dt(e_mid, pl_mid, pt_mid, t + dt, hydro);
		double dpl2 = dt *  dpl_dt(e_mid, pl_mid, pt_mid, t + dt, hydro);
		double dpt2 = dt *  dpt_dt(e_mid, pl_mid, pt_mid, t + dt, hydro);

		e += (de1 + de2) / 2.;
		pl += (dpl1 + dpl2) / 2.;
		pt += (dpt1 + dpt2) / 2.;

		// compute anisotropic variables from updated values (e,pl,pt)

	#endif

		t += dt;

		steps++;
	}
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::

	e_e0_plot.close();								// close files
	pl_pt_plot.close();
	B_plot.close();
}


void set_aniso_bjorken_initial_condition(int nx, int ny, int nz, initial_condition_parameters initial, hydro_parameters hydro)
{
#ifdef ANISO_HYDRO
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

				q[s].pl = e_s * plpt_ratio / (2. + plpt_ratio);		// conformal eos initialization
				q[s].pt = e_s / (2. + plpt_ratio);

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




