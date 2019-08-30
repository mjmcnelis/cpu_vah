#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <libconfig.h>
#include "../include/Parameters.h"
using namespace std;


void getIntegerProperty(config_t * cfg, const char * name, int * value)
{
	config_lookup_int(cfg, name, value);
}


void getDoubleProperty(config_t * cfg, const char * name, double * value)
{
	config_lookup_float(cfg, name, value);
}


bool starting_time_step_within_CFL_bound(double dt, lattice_parameters lattice)
{
	double dx = lattice.lattice_spacing_x;
	double dy = lattice.lattice_spacing_y;
	double dn = lattice.lattice_spacing_eta;

	double dt_CFL = 0.125 * fmin(dx, fmin(dy, dn));

	if(dt <= dt_CFL) return true;
	else return false;
}


double compute_conformal_prefactor(double flavors)
{
	double colors = 3.;

	return M_PI * M_PI * (2. * (colors * colors  - 1.)  +  3.5 * colors * flavors) / 30.;
}


hydro_parameters load_hydro_parameters()
{
	config_t cfg;
	config_init(&cfg);

	char fname[255] = "parameters/hydro.properties";
	if(!config_read_file(&cfg, fname))
	{
		fprintf(stderr, "No configuration file  %s found for hydrodynamic parameters - %s.\n", fname, config_error_text(&cfg));
		exit(-1);
	}

	hydro_parameters hydro;

	double quark_flavors;

	getIntegerProperty(&cfg,"run_hydro",				&hydro.run_hydro);

	getDoubleProperty(&cfg, "tau_initial", 				&hydro.tau_initial);
	getDoubleProperty(&cfg, "plpt_ratio_initial",		&hydro.plpt_ratio_initial);
	getDoubleProperty(&cfg, "quark_flavors",			&quark_flavors);

	getIntegerProperty(&cfg,"temperature_etas", 		&hydro.temperature_etas);

	getDoubleProperty(&cfg, "etas_min", 				&hydro.etas_min);
	getDoubleProperty(&cfg, "etas_slope", 				&hydro.etas_slope);
	getDoubleProperty(&cfg, "constant_etas", 			&hydro.constant_etas);
	getDoubleProperty(&cfg, "freezeout_temperature_GeV",&hydro.freezeout_temperature_GeV);
	getDoubleProperty(&cfg, "flux_limiter",				&hydro.flux_limiter);

	getIntegerProperty(&cfg,"include_vorticity", 		&hydro.include_vorticity);

	getDoubleProperty(&cfg, "energy_min",				&hydro.energy_min);
	getDoubleProperty(&cfg, "pressure_min",				&hydro.pressure_min);

	getIntegerProperty(&cfg,"regulation_scheme", 		&hydro.regulation_scheme);
	getIntegerProperty(&cfg,"reprojection", 			&hydro.reprojection);

	getDoubleProperty(&cfg, "rho_max",					&hydro.rho_max);
	getDoubleProperty(&cfg, "xi0",						&hydro.xi0);

	config_destroy(&cfg);

	hydro.conformal_eos_prefactor = compute_conformal_prefactor(quark_flavors);

	printf("\nHydro parameters:");
	printf("\n-----------------\n");
	printf("run_hydro                 = %d\n",		hydro.run_hydro);
	printf("tau_initial               = %.3g\n",	hydro.tau_initial);
	printf("plpt_ratio_initial        = %.2g\n", 	hydro.plpt_ratio_initial);
	printf("quark_flavors             = %.1f\n", 	quark_flavors);
	printf("conformal_eos_prefactor   = %.6g\n", 	hydro.conformal_eos_prefactor);
	printf("temperature_etas          = %d\n", 		hydro.temperature_etas);
	printf("etas_min                  = %.3g\n", 	hydro.etas_min);
	printf("etas_slope                = %.3g\n", 	hydro.etas_slope);
	printf("constant_etas             = %.3g\n", 	hydro.constant_etas);
	printf("freezeout_temperature_GeV = %.3f\n", 	hydro.freezeout_temperature_GeV);
	printf("flux_limiter              = %.3g\n", 	hydro.flux_limiter);
	printf("include_vorticity         = %d\n", 		hydro.include_vorticity);
	printf("energy_min                = %.1e\n", 	hydro.energy_min);
	printf("pressure_min              = %.1e\n", 	hydro.pressure_min);
	printf("regulation_scheme         = %d\n", 		hydro.regulation_scheme);
	printf("reprojection              = %d\n", 		hydro.reprojection);
	printf("rho_max                   = %.2f\n", 	hydro.rho_max);
	printf("xi0                       = %.2f\n", 	hydro.xi0);
	printf("\n");


	if(hydro.tau_initial == 0)
	{
		printf("load_hydro_parameters error: tau_initial = %.3f is not allowed\n", hydro.tau_initial);
		exit(-1);
	}

	if(fabs(hydro.flux_limiter - 1.5) > 0.5)
	{
		printf("load_hydro_parameters error: flux_limiter = %.3g is out of bounds\n", hydro.flux_limiter);
		exit(-1);
	}

	return hydro;
}


lattice_parameters load_lattice_parameters(hydro_parameters hydro)
{
	config_t cfg;
	config_init(&cfg);

	char fname[255] = "parameters/lattice.properties";
	if(!config_read_file(&cfg, fname))
	{
		fprintf(stderr, "No configuration file %s found for lattice parameters - %s.\n", fname, config_error_text(&cfg));
		exit(-1);
	}

	lattice_parameters lattice;

	getIntegerProperty(&cfg, "lattice_points_x", 		&lattice.lattice_points_x);
	getIntegerProperty(&cfg, "lattice_points_y", 		&lattice.lattice_points_y);
	getIntegerProperty(&cfg, "lattice_points_eta", 		&lattice.lattice_points_eta);

	getDoubleProperty(&cfg,  "lattice_spacing_x", 		&lattice.lattice_spacing_x);
	getDoubleProperty(&cfg,  "lattice_spacing_y", 		&lattice.lattice_spacing_y);
	getDoubleProperty(&cfg,  "lattice_spacing_eta", 	&lattice.lattice_spacing_eta);

	getIntegerProperty(&cfg, "max_time_steps",			&lattice.max_time_steps);

	getDoubleProperty(&cfg,  "output_interval", 		&lattice.output_interval);

	getDoubleProperty(&cfg,  "fixed_time_step", 		&lattice.fixed_time_step);

	getIntegerProperty(&cfg, "adaptive_time_step", 		&lattice.adaptive_time_step);

	config_destroy(&cfg);

	lattice.min_time_step = pow(10., round(log10(hydro.tau_initial)) - 2.);		// min time step ~ 100x smaller than t0 

	printf("Lattice parameters:");
	printf("\n-------------------\n");
	printf("lattice_points_x         = %d\n", 	lattice.lattice_points_x);
	printf("lattice_points_y         = %d\n", 	lattice.lattice_points_y);
	printf("lattice_points_eta       = %d\n", 	lattice.lattice_points_eta);
	printf("lattice_spacing_x        = %.3g\n", lattice.lattice_spacing_x);
	printf("lattice_spacing_y        = %.3g\n", lattice.lattice_spacing_y);
	printf("lattice_spacing_eta      = %.3g\n", lattice.lattice_spacing_eta);

	printf("max_time_steps           = %d\n", 	lattice.max_time_steps);
	printf("output_interval          = %.2f\n", lattice.output_interval);

	printf("fixed_time_step          = %.3g\n", lattice.fixed_time_step);

	printf("adaptive_time_step       = %d\n", 	lattice.adaptive_time_step);
	printf("min_time_step            = %.2e\n", lattice.min_time_step);
	printf("\n");


	double dt = lattice.fixed_time_step;						// dt = starting time step
	if(lattice.adaptive_time_step) dt = lattice.min_time_step;		

	if(!starting_time_step_within_CFL_bound(dt, lattice))
	{
		printf("load_lattice_parameters error: starting time step dt = %.2g greater than CFL bound\n", dt);
		exit(-1);
	}

	return lattice;
}


initial_condition_parameters load_initial_condition_parameters()
{
	config_t cfg;
	config_init(&cfg);

	char fname[255] = "parameters/initial.properties";
	if(!config_read_file(&cfg, fname))
	{
		fprintf(stderr, "No configuration file  %s found for initial condition parameters - %s.\n", fname, config_error_text(&cfg));
		exit(-1);
	}

	initial_condition_parameters initial;

	getIntegerProperty(&cfg, "initialConditionType", 		&initial.initialConditionType);
	getIntegerProperty(&cfg, "numberOfNucleonsPerNuclei", 	&initial.numberOfNucleonsPerNuclei);

	getDoubleProperty(&cfg, "initialCentralTemperatureGeV",	&initial.initialCentralTemperatureGeV);
	getDoubleProperty(&cfg, "scatteringCrossSectionNN", 	&initial.scatteringCrossSectionNN);
	getDoubleProperty(&cfg, "impactParameter", 				&initial.impactParameter);
	getDoubleProperty(&cfg,	"fractionOfBinaryCollisions", 	&initial.fractionOfBinaryCollisions);
	getDoubleProperty(&cfg, "rapidityVariance", 			&initial.rapidityVariance);
	getDoubleProperty(&cfg, "rapidityMean", 				&initial.rapidityMean);
	getDoubleProperty(&cfg, "q_gubser", 					&initial.q_gubser);

	config_destroy(&cfg);

	printf("Initial condition parameters:");
	printf("\n-----------------------------\n");
	printf("initialConditionType         = %d\n", 	initial.initialConditionType);
	printf("numberOfNucleonsPerNuclei    = %d\n", 	initial.numberOfNucleonsPerNuclei);
	printf("initialCentralTemperatureGeV = %.3g\n", initial.initialCentralTemperatureGeV);
	printf("scatteringCrossSectionNN     = %.3g\n", initial.scatteringCrossSectionNN);
	printf("impactParameter              = %.3g\n", initial.impactParameter);
	printf("fractionOfBinaryCollisions   = %.3g\n", initial.fractionOfBinaryCollisions);
	printf("rapidityVariance             = %.3g\n", initial.rapidityVariance);
	printf("rapidityMean                 = %.3g\n", initial.rapidityMean);
	printf("q_gubser                     = %.3g\n", initial.q_gubser);
	printf("\n");

	return initial;
}





