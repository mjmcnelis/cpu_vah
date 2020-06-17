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

double de_dt(double e, double pl, double pt, double b, double lambda, double aT, double aL, double t, hydro_parameters hydro)
{
	return - (e + pl) / t;
}


double dpl_dt_conformal(double e, double pl, double t, hydro_parameters hydro)
{
	double conformal = hydro.conformal_eos_prefactor;
	double p = e / 3.;
	double T = pow(e / conformal, 0.25);

	double etas = eta_over_s(T, hydro);
	double taupiInv = T / (5. * etas);

	double pt = (e - pl) / 2.;

	aniso_transport_coefficients aniso;
	aniso.compute_transport_coefficients(e, pl, pt, conformal);

	return - taupiInv * (pl - p)  +  aniso.zeta_LL / t;

}


double dpl_dt(double e, double pl, double pt, double db, double lambda, double aT, double aL, double t, hydro_parameters hydro)
{
	double conformal = hydro.conformal_eos_prefactor;

	equation_of_state_new eos(e, conformal);
	double p = eos.equilibrium_pressure();
	double T = eos.T;
	double s = (e + p) / T;
	double beq = eos.equilibrium_mean_field();
	double mass = T * eos.z_quasi();				     // m(T)
	double mbar = mass / lambda;       					 // m(T) / lambda
	double mdmde = eos.mdmde_quasi();

	double taupiInv = eos.beta_shear() / (s * eta_over_s(T, hydro));
	double taubulkInv = eos.beta_bulk() / (s * zeta_over_s(T, hydro));

	double b = beq + db;

	aniso_transport_coefficients_nonconformal aniso;
	aniso.compute_transport_coefficients(e, pl, pt, b, lambda, aT, aL, mbar, mass, mdmde);

	return taubulkInv * (p - (2.*pt + pl) / 3.)  -  2./3. * taupiInv * (pl - pt)  +  aniso.zeta_LL / t;
}


double dpt_dt(double e, double pl, double pt, double db, double lambda, double aT, double aL, double t, hydro_parameters hydro)
{
	double conformal = hydro.conformal_eos_prefactor;

	equation_of_state_new eos(e, conformal);
	double p = eos.equilibrium_pressure();
	double T = eos.T;
	double s = (e + p) / T;
	double beq = eos.equilibrium_mean_field();
	double mass = T * eos.z_quasi();				     // m(T)
	double mbar = mass / lambda;       					 // m(T) / lambda
	double mdmde = eos.mdmde_quasi();

	double taupiInv = eos.beta_shear() / (s * eta_over_s(T, hydro));
	double taubulkInv = eos.beta_bulk() / (s * zeta_over_s(T, hydro));

	double b = beq + db;

	aniso_transport_coefficients_nonconformal aniso;
	aniso.compute_transport_coefficients(e, pl, pt, b, lambda, aT, aL, mbar, mass, mdmde);

	return taubulkInv * (p - (2.*pt + pl) / 3.)  +  taupiInv * (pl - pt) / 3.  +  aniso.zeta_LT / t;
}


double ddb_dt(double e, double pl, double pt, double db, double lambda, double aT, double aL, double t, hydro_parameters hydro)
{
	equation_of_state_new eos(e, hydro.conformal_eos_prefactor);
	double p = eos.equilibrium_pressure();
	double T = eos.T;
	double s = (e + p) / T;
	double mass = T * eos.z_quasi();
	double beq = eos.equilibrium_mean_field();
	double mdmde = eos.mdmde_quasi();
	double taubulkInv = eos.beta_bulk() / (s * zeta_over_s(T, hydro));

	double Pi = (pl + 2.*pt)/3. - p;

	return - taubulkInv * db  -  (3.*Pi + 4.*db) * (e + pl) * mdmde / (t * mass * mass);
}


void run_semi_analytic_aniso_bjorken(lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro)
{
	// using a highly accurate 4th order numerical solution (fixed time step)
	// kinetic theory model should be quasiparticle but it can run either lattice + quasiparticle or conformal (or switch_eos mode)

	double t  = hydro.tau_initial;												// initial time
	double dt = t / 20.;														// use minimum time step 20x smaller than t0
	int decimal = - log10(dt);													// setprecision value for t output


	// initial state information
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::
	double T0  = initial.initialCentralTemperatureGeV;							// initial temperature
	double T = T0 / hbarc;
	double plpt_ratio = hydro.plpt_ratio_initial;								// initial pl/pt ratio
	double conformal_prefactor = hydro.conformal_eos_prefactor;

	double e0 = equilibrium_energy_density_new(T, conformal_prefactor);
	equation_of_state_new eos(e0, conformal_prefactor);
	double p0 = eos.equilibrium_pressure();										// initial thermal pressure
	double s0 = (e0 + p0) / T;													// initial entropy density
	double cs2 = eos.speed_of_sound_squared();
	double b0 = eos.equilibrium_mean_field();									// initial thermal mean field
	double mass = T * eos.z_quasi();											// quasiparticle mass
	double mass2 = mass * mass;
	double mdmde = eos.mdmde_quasi();

	double zetas = zeta_over_s(T, hydro);										// compute bulk relaxation time for asymptotic dB
	double betabulk = eos.beta_bulk();
	double taubulk = (zetas * s0) / betabulk;

	// initialize hydrodynamic variables (e, pl, pt, B, lambda, aL, aT)
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::

	double e = e0;
	double pl = 3. * p0 * plpt_ratio / (2. + plpt_ratio);
	double pt = 3. * p0 / (2. + plpt_ratio);
	double db = 0;

#ifdef CONFORMAL_EOS
	double b = 0;
#else
	double b = b0;																// initial mean field
	double lambda = T;															// initial guess for anisotropic variables
	double aT = 1.0;
	double aL = 1.0;

	aniso_variables X0 = find_anisotropic_variables(e, pl, pt, b, mass, lambda, aT, aL);

	if(X0.did_not_find_solution)
	{
		printf("run_semi_analytic_aniso_bjorken error: couldn't initialize anisotropic variables\n");
		exit(-1);
	}

	lambda = X0.lambda;
	aT = X0.aT;
	aL = X0.aL;
#endif
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::




	// freezeout condition
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::
	double T_freeze = hydro.freezeout_temperature_GeV;
	double e_freeze = equilibrium_energy_density_new(T_freeze / hbarc, conformal_prefactor);
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::



	// plot observables
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::
	ofstream e_e0_plot;
	ofstream pl_pt_plot;
	ofstream db_peq_plot;
	ofstream db2_peq_plot;
	ofstream dbasy_peq_plot;
	ofstream shear_peq_plot;
	ofstream bulk_peq_plot;
	ofstream pl_pt_NS_plot;
	ofstream shear_NS_plot;
	ofstream bulk_NS_plot;
	e_e0_plot.open("semi_analytic/e_e0.dat", ios::out);
	pl_pt_plot.open("semi_analytic/pl_pt.dat", ios::out);
	db_peq_plot.open("semi_analytic/db_peq.dat", ios::out);
	db2_peq_plot.open("semi_analytic/db2_peq.dat", ios::out);
	dbasy_peq_plot.open("semi_analytic/dbasy_peq.dat", ios::out);
	bulk_peq_plot.open("semi_analytic/bulk_peq.dat", ios::out);
	shear_peq_plot.open("semi_analytic/shear_peq.dat", ios::out);
	pl_pt_NS_plot.open("semi_analytic/pl_pt_NS.dat", ios::out);
	bulk_NS_plot.open("semi_analytic/bulk_NS.dat", ios::out);
	shear_NS_plot.open("semi_analytic/shear_NS.dat", ios::out);
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::


	// start evolution
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::
	while(true)
	{
		equation_of_state_new eos_start(e, conformal_prefactor);
		double T_s = eos_start.T;
		double p = eos_start.equilibrium_pressure();
		double beq = eos_start.equilibrium_mean_field();
		double s = (e + p) / T_s;
		double eta = s * eta_over_s(T_s, hydro);
		double zeta = s * zeta_over_s(T_s, hydro);
		double mdmde_s = eos_start.mdmde_quasi();
		double m = T_s * eos_start.z_quasi();
		double taubulk_s = zeta / eos_start.beta_bulk();

		double pl_NS = p - zeta/t - 4.*eta/(3.*t);
		double pt_NS = p - zeta/t + 2.*eta/(3.*t);

		double db2nd = - taubulk_s * mdmde_s * (e + pl) * (2.*pt + pl - 3.*p) / (t * m * m);

		double dbasy = - taubulk_s * mdmde_s * (e + pl) * (2.*pt + pl - 3.*p) / (t * m * m) / (1. + 4. * taubulk_s * mdmde_s * (e + pl) / (t * m * m));

		e_e0_plot  		<< fixed << setprecision(decimal + 2) << t << "\t" << scientific << setprecision(12) << e / e0 << endl;
		pl_pt_plot 		<< fixed << setprecision(decimal + 2) << t << "\t" << scientific << setprecision(12) << pl / pt << endl;
		db_peq_plot		<< fixed << setprecision(decimal + 2) << t << "\t" << scientific << setprecision(12) << db / p << endl;
		db2_peq_plot	<< fixed << setprecision(decimal + 2) << t << "\t" << scientific << setprecision(12) << db2nd / p << endl;
		dbasy_peq_plot	<< fixed << setprecision(decimal + 2) << t << "\t" << scientific << setprecision(12) << dbasy / p << endl;
		// db_peq_plot		<< fixed << setprecision(decimal + 2) << t << "\t" << scientific << setprecision(12) << (B - beq) * hbarc << endl;
		// dbasy_peq_plot	<< fixed << setprecision(decimal + 2) << t << "\t" << scientific << setprecision(12) << dbasy * hbarc << endl;
		bulk_peq_plot	<< fixed << setprecision(decimal + 2) << t << "\t" << scientific << setprecision(12) << (pl + 2.*pt) / (3.*p) - 1. << endl;
		shear_peq_plot	<< fixed << setprecision(decimal + 2) << t << "\t" << scientific << setprecision(12) << 2.*(pl - pt) / (3.*p) << endl;
		pl_pt_NS_plot	<< fixed << setprecision(decimal + 2) << t << "\t" << scientific << setprecision(12) << pl_NS / pt_NS << endl;
		bulk_NS_plot	<< fixed << setprecision(decimal + 2) << t << "\t" << scientific << setprecision(12) << -zeta / (p*t) << endl;
		shear_NS_plot	<< fixed << setprecision(decimal + 2) << t << "\t" << scientific << setprecision(12) << -4.*eta / (3.*p*t) << endl;

		if(e < e_freeze) break;


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
	#else

		double de1  = dt *  de_dt(e, pl, pt, db, lambda, aT, aL, t, hydro);
		double dpl1 = dt * dpl_dt(e, pl, pt, db, lambda, aT, aL, t, hydro);
		double dpt1 = dt * dpt_dt(e, pl, pt, db, lambda, aT, aL, t, hydro);
		double ddb1  = dt *  ddb_dt(e, pl, pt, db, lambda, aT, aL, t, hydro);

		double e_mid = e + de1;
		double pl_mid = pl + dpl1;
		double pt_mid = pt + dpt1;
		double db_mid = db + ddb1;





		equation_of_state_new eos_mid(e_mid, conformal_prefactor);
		double T_mid = eos_mid.T;
		double mass_mid = T_mid * eos_mid.z_quasi();
		double beq_mid = eos_mid.equilibrium_mean_field();



		// regulating the mean field
		//db_mid *= fmin(1., fabs(beq_mid / db_mid));

		if((beq_mid >= 0 && db_mid > beq_mid) || (beq_mid < 0 && db_mid < beq_mid))
		{
			db_mid = beq_mid;
		}


		double b_mid = beq_mid + db_mid;

		aniso_variables X_mid = find_anisotropic_variables(e_mid, pl_mid, pt_mid, b_mid, mass_mid, lambda, aT, aL);

		lambda = X_mid.lambda;
		aT = X_mid.aT;
		aL = X_mid.aL;



		double de2  = dt *  de_dt(e_mid, pl_mid, pt_mid, db_mid, lambda, aT, aL, t + dt, hydro);
		double dpl2 = dt * dpl_dt(e_mid, pl_mid, pt_mid, db_mid, lambda, aT, aL, t + dt, hydro);
		double dpt2 = dt * dpt_dt(e_mid, pl_mid, pt_mid, db_mid, lambda, aT, aL, t + dt, hydro);
		double ddb2  = dt *  ddb_dt(e_mid, pl_mid, pt_mid, db_mid, lambda, aT, aL, t + dt, hydro);

		e  += (de1  + de2)  / 2.;
		pl += (dpl1 + dpl2) / 2.;
		pt += (dpt1 + dpt2) / 2.;
		db += (ddb1  + ddb2)  / 2.;



		equation_of_state_new eos_end(e, conformal_prefactor);
		double T_end = eos_end.T;
		double mass_end = T_end * eos_end.z_quasi();
		double beq_end = eos_end.equilibrium_mean_field();

		// regulating the mean field
		//db *= fmin(1., fabs(beq_end / db));
		if((beq_end >= 0 && db > beq_end) || (beq_end < 0 && db < beq_end))
		{
			db = beq_end;
		}

		b = beq_end + db;

		aniso_variables X_end = find_anisotropic_variables(e, pl, pt, b, mass_end, lambda, aT, aL);

		lambda = X_end.lambda;
		aT = X_end.aT;
		aL = X_end.aL;

	#endif

		t += dt;
	}

	e_e0_plot.close();								// close files
	pl_pt_plot.close();
	db_peq_plot.close();
	db2_peq_plot.close();
	dbasy_peq_plot.close();
	shear_peq_plot.close();
	bulk_peq_plot.close();
	pl_pt_NS_plot.close();
	shear_NS_plot.close();
	bulk_NS_plot.close();
}


void set_aniso_bjorken_initial_condition(int nx, int ny, int nz, initial_condition_parameters initial, hydro_parameters hydro)
{
#ifdef ANISO_HYDRO
	precision t = hydro.tau_initial;
	precision plpt_ratio = hydro.plpt_ratio_initial;
	precision e_min = hydro.energy_min;
	precision conformal_eos_prefactor = hydro.conformal_eos_prefactor;

	precision T0 = initial.initialCentralTemperatureGeV / hbarc;				// central temperature (fm^-1)
	//precision e0 = equilibrium_energy_density(T0, conformal_eos_prefactor);		// energy density
	precision e0 = equilibrium_energy_density_new(T0, conformal_eos_prefactor);		// energy density


	equation_of_state_new eos(e0, conformal_eos_prefactor);													// compute mean field components
	precision b0 = eos.equilibrium_mean_field();								// thermal mean field
	precision p0 = eos.equilibrium_pressure();
	precision mass = T0 * eos.z_quasi();

	precision pl0 = 3. * p0 * plpt_ratio / (2. + plpt_ratio);						// conformal initialization of pl, pt
	precision pt0 = 3. * p0 / (2. + plpt_ratio);

	// equation_of_state eos(e0);													// compute mean field components
	// precision b0 = eos.equilibrium_mean_field(T0);								// thermal mean field
	// precision p0 = eos.equilibrium_pressure();
	// precision mass = T0 * eos.z_quasi(T0);										// quasiparticle mass and derivative
	// precision mdmde = eos.mdmde_quasi();
	// precision taubulk = zeta_over_s(T0, hydro) * (e0 + p0) / (T0 * eos.beta_bulk(T0));		// bulk relaxation time

	// only for testing equilibrium initial conditions with lattice eos only (this sets db0_asy = 0)
	//precision pl0 = p0;
	//precision pt0 = p0;

	// bjorken asymptotic non-equilibrium mean field component
	//precision db0 = -3. * taubulk * mdmde * (e0 + pl0) * (2.*pt0/3. + pl0/3. - p0) / (t * mass * mass) / (1.  +  4. * taubulk * mdmde * (e0 + pl0) / (t * mass * mass));
	//precision db0 = - 3. * taubulk * mdmde * (e0 + p0) * (2.*pt0/3. + pl0/3. - p0) / (t * mass * mass);

	// printf("b0  = %lf fm^-4\n", b0);
	// printf("db = %lf fm^-4\n", db0_asy);
	// printf("b  = %lf fm^-4\n", b0 + db0);

	// printf("\ne = %lf\n", e0);
	// printf("T = %lf\n", T0);
	// printf("pl = %lf\n", pl0);
	// printf("pt = %lf\n", pt0);
	// printf("beq = %lf\n", b0);
	// printf("db = %lf\n", db0_asy);
	// printf("b  = %lf\n", b0 + db0);


	// initialize anisotropic variables
	precision lambda0 = T0;														// initial guess
	precision aT0 = 1.0;
	precision aL0 = 1.0;

	aniso_variables X = find_anisotropic_variables(e0, pl0, pt0, b0, mass, lambda0, aT0, aL0);

	lambda0 = X.lambda;
	aT0 = X.aT;
	aL0 = X.aL;

	// printf("lambda  = %lf\n", lambda0);
	// printf("aT      = %lf\n", aT0);
	// printf("aL      = %lf\n", aL0);

	// leave out openmp for now
	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				e[s] = e0;
				q[s].pl = pl0;
				q[s].pt = pt0;

			#ifdef B_FIELD
				//q[s].b = b0 + db0;
				q[s].b = b0;
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




