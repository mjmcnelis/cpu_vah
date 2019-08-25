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
#include "../include/Parameters.h"
using namespace std;


double de_dt(double e, double pl, double t, hydro_parameters hydro)
{
	return - (e + pl) / t;
}


double dpl_dt(double e, double pl, double t, hydro_parameters hydro)
{
	double p = equilibriumPressure(e);
	double T = effectiveTemperature(e);

	double etas = eta_over_s(T, hydro);
	double taupiInv = T / (5. * etas);

	double pt = (e - pl) / 2.; 	// temporary

	transport_coefficients aniso;
	aniso.compute_transport_coefficients(e, pl, pt);


	return - taupiInv * (pl - p)  +  aniso.zeta_LL / t;

}


void run_semi_analytic_aniso_bjorken(lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro)
{
	// using a highly accurate 4th order numerical solution (fixed time step)
	// should be able to handle both conformal and nonconformal Bjorken (keep it at conformal for now)

	double t  = hydro.tau_initial;						// initial time
	double T0  = initial.initialCentralTemperatureGeV;	// initial temperature
	double e0 = equilibriumEnergyDensity(T0 / hbarc);	// initial energy density
	double dt = pow(10., round(log10(t)) - 2.);			// time step ~ 100x smaller than t0 (rounded to nearest decimal)

	double plpt_ratio = hydro.plpt_ratio_initial;		// initial pl/pt ratio

	double e = e0;
	double pl = e0 * plpt_ratio / (2. + plpt_ratio);
	
	double T_freeze = hydro.freezeout_temperature_GeV;
	double e_freeze = equilibriumEnergyDensity(T_freeze / hbarc);

	ofstream e_e0, pl_pt;
	e_e0.open("semi_analytic/e_e0_aniso_bjorken.dat", ios::out);
	pl_pt.open("semi_analytic/pl_pt_aniso_bjorken.dat", ios::out);

	while(true)
	{
		double pt = (e - pl) / 2.;	// temporary

		e_e0  << setprecision(8) << t << "\t" << e / e0 << endl;
		pl_pt << setprecision(8) << t << "\t" << pl / pt << endl;

		if(e < e_freeze) break;

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
	}

	e_e0.close();
	pl_pt.close();
}



