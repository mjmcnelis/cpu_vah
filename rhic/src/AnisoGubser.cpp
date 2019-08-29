#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include "../include/Hydrodynamics.h"
#include "../include/EquationOfState.h"
#include "../include/TransportCoefficients.h"
#include "../include/Parameters.h"
#include "../include/AnisoGubser.h"
using namespace std;



double rho_function(double t, double r, double q)
{
	return - asinh((1.  -  q * q * (t * t  -  r * r)) / (2. * q * t));
}


double de_drho(double e, double pl, double rho, hydro_parameters hydro)
{
	return (pl  -  3.*e) * tanh(rho);
}


double dpl_drho(double e, double pl, double rho, hydro_parameters hydro)
{
	double conformal_eos_prefactor = hydro.conformal_eos_prefactor;
	double T = effectiveTemperature(e, conformal_eos_prefactor);
	double etas = eta_over_s(T, hydro);
	double taupiInv = T / (5. * etas);

	transport_coefficients aniso;
	aniso.compute_transport_coefficients(e, pl, (e - pl)/2., conformal_eos_prefactor);

	return - taupiInv * (pl - e/3.)  -  (4. * pl  +  aniso.zeta_LL) * tanh(rho);

}


void gubser_rho_evolution(double * e_hat, double * pl_hat, double * rho_array, int rho_pts, double drho, hydro_parameters hydro)
{
	for(int i = 1; i < rho_pts; i++)
	{
		double rho = rho_array[i - 1];

		double ep  =  e_hat[i - 1];
		double plp = pl_hat[i - 1];

		double de1  = drho *  de_drho(ep, plp, rho, hydro);
		double dpl1 = drho * dpl_drho(ep, plp, rho, hydro);

		double de2  = drho *  de_drho(ep + de1/2., plp + dpl1/2., rho + drho/2., hydro);
		double dpl2 = drho * dpl_drho(ep + de1/2., plp + dpl1/2., rho + drho/2., hydro);

		double de3  = drho *  de_drho(ep + de2/2., plp + dpl2/2., rho + drho/2., hydro);
		double dpl3 = drho * dpl_drho(ep + de2/2., plp + dpl2/2., rho + drho/2., hydro);

		double de4  = drho *  de_drho(ep + de3, plp + dpl3, rho + drho, hydro);
		double dpl4 = drho * dpl_drho(ep + de3, plp + dpl3, rho + drho, hydro);

		e_hat[i]  = ep   +  (de1   +  2. * de2   +  2. * de3   +  de4)  / 6.;
		pl_hat[i] = plp  +  (dpl1  +  2. * dpl2  +  2. * dpl3  +  dpl4) / 6.;
	}
}


void run_semi_analytic_aniso_gubser(lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;

	double dx = lattice.lattice_spacing_x;
	double dy = lattice.lattice_spacing_y;

	double t = hydro.tau_initial;					// initial longitudinal proper time

	double dt_out = lattice.output_interval;		// output times

	double q0 = 1.0;								// inverse length size (hard coded)
	double T0_hat = initial.T_hat_initial;			// value hard coded in parameters/initial.properties
	double plpt_ratio = hydro.plpt_ratio_initial;	// initial plpt ratio at transverse corner of grid

	// initial T_hat hard coded so that initial central temperature = 0.6 GeV (make a parameter)
	//double T0_hat = 0.0455468;   			// plpt_ratio = 1.0
	//double T0_hat = 0.04296357;			// plpt_ratio = 0.01  (x,y = 5fm  x 5fm)
	//double T0_hat = 0.0328418;			// plpt_ratio = 0.01  (x,y = 6fm  x 6fm)
	//double T0_hat = 0.0261391;	  		// plpt ratio = 0.01  (x,y = 7fm x 7fm)
	//double T0_hat = 0.01537397;	  		// plpt_ratio = 0.01  (x,y = 10fm x 10fm)
	//double T0_hat = 0.01171034;	  		// plpt_ratio = 0.01  (x,y = 12fm x 12fm)

	double x_max = (nx - 1.) * dx / 2.;
	double y_max = (ny - 1.) * dy / 2.;

	double r_max = sqrt(x_max * x_max  +  y_max * y_max);	// distance to transverse corner

    double rho0 = rho_function(t, r_max, q0);				// min rho at transverse corner
	double rhoP = rho_function(t, 0, q0);					// max rho at center

	double drho = 0.0001;

	int rho_pts = ceil(fabs((rhoP - rho0) / drho));

	drho = fabs((rhoP - rho0) / ((double)rho_pts - 1.));

	double rho_array[rho_pts];

	for(int j = 0; j < rho_pts; j++)
	{
		rho_array[j] = rho0  +  j * drho;
	}

	double  e_hat[rho_pts];
	double pl_hat[rho_pts];

	 e_hat[0] = equilibriumEnergyDensity(T0_hat, hydro.conformal_eos_prefactor);
	pl_hat[0] = e_hat[0] * plpt_ratio / (2. + plpt_ratio);

	gubser_rho_evolution(e_hat, pl_hat, rho_array, rho_pts, drho, hydro);

	gsl_spline * e_hat_spline;		// construct the cubic spline interpolations
	gsl_spline * pl_hat_spline;

	e_hat_spline = gsl_spline_alloc(gsl_interp_cspline, rho_pts);
	pl_hat_spline = gsl_spline_alloc(gsl_interp_cspline, rho_pts);

	gsl_spline_init(e_hat_spline, rho_array, e_hat, rho_pts);
	gsl_spline_init(pl_hat_spline, rho_array, pl_hat, rho_pts);



	// freezeout conditions
	double T_freeze = hydro.freezeout_temperature_GeV;
	double e_freeze = equilibriumEnergyDensity(T_freeze / hbarc, hydro.conformal_eos_prefactor);

	double eps = 1.e-5;

	//t += dt_out;

	while(true)
	{
		gsl_interp_accel * accel = gsl_interp_accel_alloc();

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

			double e_r  = gsl_spline_eval(e_hat_spline,  rho, accel) / (t * t * t * t);
			double pl = gsl_spline_eval(pl_hat_spline, rho, accel) / (t * t * t * t);

			double T = effectiveTemperature(e_r, hydro.conformal_eos_prefactor);

			if((i - 2 == (nx - 1) / 2)) printf("%lf\t%lf\t%lf\n", t, e_r * hbarc, T * hbarc);

			fprintf(energy, 	"%.3f\t%.8f\n", x, e_r);
			fprintf(plptratio,	"%.3f\t%.8f\n", x, 2. * pl / (e_r - pl));

			if(e_r > e_freeze) below_freezeout_surface = false;
		}

		fclose(energy);
		fclose(plptratio);

		if(below_freezeout_surface) break;

		t += dt_out;

		gsl_interp_accel_free(accel);
	}

	gsl_spline_free(e_hat_spline);
	gsl_spline_free(pl_hat_spline);

}



















