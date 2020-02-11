#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <cmath>
#include "../include/EquationOfState.h"
#include "../include/TransportAniso.h"
#include "../include/TransportAnisoNonconformal.h"
#include "../include/AnisoVariables.h"
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

double de_dt(double e, double pl, double pt, double B, double lambda, double aT, double aL, double t, hydro_parameters hydro)
{
	return - (e + pl) / t;
}


double dpl_dt_conformal(double e, double pl, double t, hydro_parameters hydro)
{
	double conformal = hydro.conformal_eos_prefactor;

	equation_of_state eos(e);
	double p = eos.equilibrium_pressure();
	double T = eos.effective_temperature(conformal);

	double etas = eta_over_s(T, hydro);
	double taupiInv = T / (5. * etas);

	double pt = (e - pl) / 2.; 	// temporary

	aniso_transport_coefficients aniso;
	aniso.compute_transport_coefficients(e, pl, pt, conformal);

	return - taupiInv * (pl - p)  +  aniso.zeta_LL / t;

}


double dpl_dt(double e, double pl, double pt, double B, double lambda, double aT, double aL, double t, hydro_parameters hydro)
{
	double conformal = hydro.conformal_eos_prefactor;

	equation_of_state eos(e);
	double p = eos.equilibrium_pressure();
	double T = eos.effective_temperature(conformal);
	double s = (e + p) / T;
	double mass = T * eos.z_quasi(T);				     // m(T)
	double mbar = mass / lambda;       					 // m(T) / lambda
	double mdmde = eos.mdmde_quasi();

	double taupiInv = eos.beta_shear(T, conformal) / (s * eta_over_s(T, hydro));
	double taubulkInv = eos.beta_bulk(T) / (s * zeta_over_s(T, hydro));

	aniso_transport_coefficients_nonconformal aniso;
	aniso.compute_transport_coefficients(e, pl, pt, B, lambda, aT, aL, mbar, mass, mdmde);

	return taubulkInv * (p - (2.*pt + pl) / 3.)  -  2./3. * taupiInv * (pl - pt)  +  aniso.zeta_LL / t;
}


double dpt_dt(double e, double pl, double pt, double B, double lambda, double aT, double aL, double t, hydro_parameters hydro)
{
	double conformal = hydro.conformal_eos_prefactor;

	equation_of_state eos(e);
	double p = eos.equilibrium_pressure();
	double T = eos.effective_temperature(conformal);
	double s = (e + p) / T;
	double mass = T * eos.z_quasi(T);				     // m(T)
	double mbar = mass / lambda;       					 // m(T) / lambda
	double mdmde = eos.mdmde_quasi();

	double taupiInv = eos.beta_shear(T, conformal) / (s * eta_over_s(T, hydro));
	double taubulkInv = eos.beta_bulk(T) / (s * zeta_over_s(T, hydro));

	aniso_transport_coefficients_nonconformal aniso;
	aniso.compute_transport_coefficients(e, pl, pt, B, lambda, aT, aL, mbar, mass, mdmde);

	return taubulkInv * (p - (2.*pt + pl) / 3.)  +  taupiInv * (pl - pt) / 3.  +  aniso.zeta_LT / t;
}


double dB_dt(double e, double pl, double pt, double B, double lambda, double aT, double aL, double t, hydro_parameters hydro)
{
	double conformal = hydro.conformal_eos_prefactor;

	equation_of_state eos(e);
	double p = eos.equilibrium_pressure();
	double T = eos.effective_temperature(conformal);
	double s = (e + p) / T;
	double mass = T * eos.z_quasi(T);
	double Beq = eos.equilibrium_mean_field(T);
	double mdmde = eos.mdmde_quasi();
	double taubulkInv = eos.beta_bulk(T) / (s * zeta_over_s(T, hydro));

	return - taubulkInv * (B - Beq)  -  (2.*pt + pl - e + 4.*B) * mdmde * (e + pl) / (t * mass * mass);
}


void run_semi_analytic_aniso_bjorken(lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro)
{
	// using a highly accurate 4th order numerical solution (fixed time step)
	// kinetic theory model should be quasiparticle but it can run either lattice + quasiparticle or conformal (or switch_eos mode)

	if(hydro.kinetic_theory_model != 1)
	{
		printf("run_semi_analytic_aniso_bjorken error: set kinetic_theory_model = 1 to default value\n");
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
	double mass = T * eos.z_quasi(T);											// quasiparticle mass
	double mass2 = mass * mass;
	double mdmde = eos.mdmde_quasi();

	double zetas = zeta_over_s(T, hydro);										// compute bulk relaxation time for asymptotic dB
	double betabulk = eos.beta_bulk(T);
	double taubulk = (zetas * s0) / betabulk;

	// printf("\n");
	// printf("e = %lf fm^-4\n", e0);
	// printf("p = %lf fm^-4\n", p0);
	// printf("s = %lf fm^-4\n", s0);
	//printf("Beq = %lf fm^-4\n", B0);
	// printf("mass = %lf\n", mass);
	// printf("mdmde = %lf\n", mdmde);
	// printf("zetaS = %lf\n", zetas);
	// printf("betabulk = %lf fm^-4\n", betabulk);
	// printf("taubulk = %lf fm\n\n", taubulk);
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::


	// initialize hydrodynamic variables (e, pl, pt, B, lambda, aL, aT)
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::

	double e = e0;

#ifdef CONFORMAL_EOS
	double pl = e0 * plpt_ratio / (2. + plpt_ratio);								// using the conformal formula I think
	double pt = e0 / (2. + plpt_ratio);
#endif

#ifdef LATTICE_QCD
	// double pl = p0;																	// start with equilibrium initial conditions for now
	// double pt = p0;
	double pl = e0 * plpt_ratio / (2. + plpt_ratio);								// using the conformal formula I think
	double pt = e0 / (2. + plpt_ratio);

	// asymptotic formula for non-equilibrium mean field component dB
	double dB0asy = -3. * taubulk * mdmde * (e + pl) * (2.*pt/3. + pl/3. - p0) / (t * mass2) / (1.  +  4. * taubulk * mdmde * (e + pl) / (t * mass2));														// keep dB = 0 for now
																				// simulation needs different initialization formula than Bjorken

	double B = B0 + dB0asy;														// initial mean field

	double lambda = T;															// initial guess for anisotropic variables
	double aT = 1.0;
	double aL = 1.0;

	aniso_variables X0 = find_anisotropic_variables(e, pl, pt, B, mass, lambda, aT, aL);

	lambda = X0.lambda;
	aT = X0.aT;
	aL = X0.aL;
#endif

	// printf("pl = %lf fm^-4\n", pl);
	// printf("pt = %lf fm^-4\n", pt);
	// printf("dB = %lf fm^-4\n", dB0asy);
	// printf("B  = %lf fm^-4\n", B);
	// printf("lambda = %lf fm^-1\n\n", lambda);
	//exit(-1);

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

	clock_t begin;
    double duration;
    begin = clock();

	// start evolution
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::
	while(true)
	{
		if(steps % 1 == 0)			// write to file
		{
			e_e0_plot  << fixed << setprecision(decimal + 1) << t << "\t" << scientific << setprecision(12) << e / e0 << endl;
			pl_pt_plot << fixed << setprecision(decimal + 1) << t << "\t" << scientific << setprecision(12) << pl / pt << endl;
			B_plot     << fixed << setprecision(decimal + 1) << t << "\t" << scientific << setprecision(12) << B << endl;
		
			if(e < e_freeze) break;
		}

	// I need to think about how to use these macro if statements correctly
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
		pt = (e - pl) / 2.;
	#endif

	#ifdef LATTICE_QCD

		// let's just use the RK2 method for now since the old code uses that
		double de1  = dt *  de_dt(e, pl, pt, B, lambda, aT, aL, t, hydro);
		double dpl1 = dt * dpl_dt(e, pl, pt, B, lambda, aT, aL, t, hydro);
		double dpt1 = dt * dpt_dt(e, pl, pt, B, lambda, aT, aL, t, hydro);
		double dB1  = dt *  dB_dt(e, pl, pt, B, lambda, aT, aL, t, hydro);

		double e_mid = e + de1;
		double pl_mid = pl + dpl1;
		double pt_mid = pt + dpt1;
		double B_mid = B + dB1;

		// printf("e_mid = %lf\n", e_mid);
		// printf("pl_mid = %lf\n", pl_mid);
		// printf("pt_mid = %lf\n", pt_mid);
		// printf("B_mid = %lf\n\n", B_mid);

		equation_of_state eos_mid(e_mid);
		double T_mid = eos_mid.effective_temperature(hydro.conformal_eos_prefactor);
		double mass_mid = T_mid * eos_mid.z_quasi(T_mid);

		aniso_variables X_mid = find_anisotropic_variables(e_mid, pl_mid, pt_mid, B_mid, mass_mid, lambda, aT, aL);

		lambda = X_mid.lambda;
		aT = X_mid.aT;
		aL = X_mid.aL;
		// exit(-1);

		// compute anisotropic variables from intermediate values (e_mid, pl_mid, pt_mid)

		double de2  = dt *  de_dt(e_mid, pl_mid, pt_mid, B_mid, lambda, aT, aL, t + dt, hydro);
		double dpl2 = dt * dpl_dt(e_mid, pl_mid, pt_mid, B_mid, lambda, aT, aL, t + dt, hydro);
		double dpt2 = dt * dpt_dt(e_mid, pl_mid, pt_mid, B_mid, lambda, aT, aL, t + dt, hydro);
		double dB2  = dt *  dB_dt(e_mid, pl_mid, pt_mid, B_mid, lambda, aT, aL, t + dt, hydro);

		e  += (de1  + de2)  / 2.;
		pl += (dpl1 + dpl2) / 2.;
		pt += (dpt1 + dpt2) / 2.;
		B  += (dB1  + dB2)  / 2.;

		// printf("e_end = %lf\n", e);
		// printf("pl_end = %lf\n", pl);
		// printf("pt_end = %lf\n", pt);
		// printf("B_end = %lf\n\n", B);
		// exit(-1);

		equation_of_state eos_end(e);
		double T_end = eos_end.effective_temperature(hydro.conformal_eos_prefactor);
		double mass_end = T_end * eos_mid.z_quasi(T_end);

		aniso_variables X_end = find_anisotropic_variables(e, pl, pt, B, mass_end, lambda, aT, aL);

		lambda = X_end.lambda;
		aT = X_end.aT;
		aL = X_end.aL;
	#endif

		t += dt;

		steps++;

		if(steps > 100000)
		{
			printf("Aniso bjorken failed to reach freezeout temperature\n");
			exit(-1);
		}
	}

	duration = (clock() - begin) / (double)CLOCKS_PER_SEC;
    printf("Average time/step = %.4g ms\n", duration / (double)steps * 1000.);
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::

	e_e0_plot.close();								// close files
	pl_pt_plot.close();
	B_plot.close();
}


void set_aniso_bjorken_initial_condition(int nx, int ny, int nz, initial_condition_parameters initial, hydro_parameters hydro)
{
#ifdef ANISO_HYDRO
	precision t = hydro.tau_initial;	
	precision plpt_ratio = hydro.plpt_ratio_initial;
	precision e_min = hydro.energy_min;
	precision conformal_eos_prefactor = hydro.conformal_eos_prefactor;

	precision T0 = initial.initialCentralTemperatureGeV / hbarc;				// central temperature (fm^-1)
	precision e0 = equilibrium_energy_density(T0, conformal_eos_prefactor);		// energy density

	if(e0 < e_min)
	{
		printf("set_aniso_bjorken_initial_condition error: initial temperature is too low\n");
		exit(-1);
	}
	else
	{
		e0 = energy_density_cutoff(e_min, e0);
	}
	
	precision pl0 = e0 * plpt_ratio / (2. + plpt_ratio);						// conformal initialization of pl, pt 
	precision pt0 =  e0 / (2. + plpt_ratio);

	equation_of_state eos(e0);													// compute mean field components
	T0 = eos.effective_temperature(hydro.conformal_eos_prefactor);				// readjust T0 (doesn't really do anything)
	precision b0 = eos.equilibrium_mean_field(T0);								// thermal mean field
	precision p0 = eos.equilibrium_pressure();																		
	precision mass = T0 * eos.z_quasi(T0);										// quasiparticle mass and derivative
	precision mdmde = eos.mdmde_quasi();
	precision taubulk = zeta_over_s(T0, hydro) * (e0 + p0) / (T0 * eos.beta_bulk(T0));		// bulk relaxation time 

	// only for testing equilibrium initial conditions with lattice eos only (this sets db0_asy = 0)
	// pl0 = p0;
	// pt0 = p0;

	// bjorken asymptotic non-equilibrium mean field component
	precision db0_asy = -3. * taubulk * mdmde * (e0 + pl0) * (2.*pt0/3. + pl0/3. - p0) / (t * mass * mass) / (1.  +  4. * taubulk * mdmde * (e0 + pl0) / (t * mass * mass));		

	// printf("b0  = %lf fm^-4\n", b0);
	// printf("db = %lf fm^-4\n", db0_asy);
	// printf("b  = %lf fm^-4\n", b0 + db0_asy);

	// initialize anisotropic variables
	precision lambda0 = T0;														// initial guess
	precision aT0 = 1.0;
	precision aL0 = 1.0;

	aniso_variables X = find_anisotropic_variables(e0, pl0, pt0, b0 + db0_asy, mass, lambda0, aT0, aL0);

	lambda0 = X.lambda;
	aT0 = X.aT;
	aL0 = X.aL;

	printf("lambda  = %lf fm^-1\n", lambda0);
	printf("aT      = %lf\n", aT0);
	printf("aL      = %lf\n", aL0);

	for(int i = 2; i < nx + 2; i++)
	{
		//printf("i = %d\n", i);

		for(int j = 2; j < ny + 2; j++)
		{
			for(int k = 2; k < nz + 2; k++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				e[s] = e0;
				q[s].pl = pl0; 	
				q[s].pt = pt0;

			#ifdef B_FIElD
				q[s].b = b0 + db0_asy;
			#endif

				// need to initialize anisotropic variables
			#ifdef LATTICE_QCD
				lambda[s] = lambda0;
				aT[s] = aT0;
				aL[s] = aL0;
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




