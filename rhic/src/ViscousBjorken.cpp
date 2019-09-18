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


double de_dt_vh(double e, double pi, double t, hydro_parameters hydro)
{
	double p = equilibriumPressure(e);

	return - (e + p - pi) / t;
}


double dpi_dt(double e, double pi, double t, hydro_parameters hydro)
{	
	double p = equilibriumPressure(e);
	double T = effectiveTemperature(e, hydro.conformal_eos_prefactor);
	double s = (e + p) / T;

	double etas = eta_over_s(T, hydro);
	double eta = s * etas;
	double taupiInv = T / (5. * etas);

	double taupipi = 10./7.;
	double deltapipi = 4./3.;
	double lambdapiPi = 1.2;

	// 2.0*lambdapiPi*Pi/3.0/tau;

	return taupiInv * (-pi  +  4./3. * eta / t)  -  (taupipi/3. + deltapipi) * pi / t;
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
	double e0 = equilibriumEnergyDensity(T0 / hbarc, conformal_eos_prefactor);	// initial energy density

	double plpt_ratio = hydro.plpt_ratio_initial;								// initial pl/pt ratio

	double e = e0;
	double p = equilibriumPressure(e);
	double pl = e0 * plpt_ratio / (2. + plpt_ratio);			
	double pt = (e - pl) / 2.;				

	double pi = 2./3. * (pl - pt);												// - t2.pinn
	double bulkPi = 0;															// set bulk pressure = 0 for now
	
	double T_freeze = hydro.freezeout_temperature_GeV;
	double e_freeze = equilibriumEnergyDensity(T_freeze / hbarc, conformal_eos_prefactor);

	ofstream e_e0_plot;
	ofstream pl_pt_plot;
	e_e0_plot.open("semi_analytic/e_e0_viscous_bjorken.dat", ios::out);
	pl_pt_plot.open("semi_analytic/pl_pt_viscous_bjorken.dat", ios::out);

	while(true)
	{
		p = equilibriumPressure(e);

		pl = p + bulkPi - pi;
		pt = p + bulkPi + pi/2.; 
		
		e_e0_plot  << fixed << setprecision(decimal + 1) << t << "\t" << scientific << setprecision(12) << e / e0 << endl;
		pl_pt_plot << fixed << setprecision(decimal + 1) << t << "\t" << scientific << setprecision(12) << pl / pt << endl;

		if(e < e_freeze) break;
		
		double e1  = dt *  de_dt_vh(e, pi, t, hydro);
		double pi1 = dt * dpi_dt(e, pi, t, hydro);

		double e2  = dt *  de_dt_vh(e + e1/2., pi + pi1/2., t + dt/2., hydro);
		double pi2 = dt * dpi_dt(e + e1/2., pi + pi1/2., t + dt/2., hydro);

		double e3  = dt *  de_dt_vh(e + e2/2., pi + pi2/2., t + dt/2., hydro);
		double pi3 = dt * dpi_dt(e + e2/2., pi + pi2/2., t + dt/2., hydro);

		double e4  = dt *  de_dt_vh(e + e3, pi + pi3, t + dt, hydro);
		double pi4 = dt * dpi_dt(e + e3, pi + pi3, t + dt, hydro);

		e  += (e1   +  2. * e2   +  2. * e3   +  e4)  / 6.;
		pl += (pi1  +  2. * pi2  +  2. * pi3  +  pi4) / 6.;		

		t += dt;
	}

	e_e0_plot.close();
	pl_pt_plot.close();
}



