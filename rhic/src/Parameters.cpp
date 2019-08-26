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


lattice_parameters load_lattice_parameters()
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
	getIntegerProperty(&cfg, "max_number_of_time_steps",&lattice.max_number_of_time_steps);

	getDoubleProperty(&cfg,  "lattice_spacing_x", 		&lattice.lattice_spacing_x);
	getDoubleProperty(&cfg,  "lattice_spacing_y", 		&lattice.lattice_spacing_y);
	getDoubleProperty(&cfg,  "lattice_spacing_eta", 	&lattice.lattice_spacing_eta);
	getDoubleProperty(&cfg,  "fixed_time_step", 		&lattice.fixed_time_step);

	getIntegerProperty(&cfg, "output_period", 			&lattice.output_period);
	getIntegerProperty(&cfg, "adaptive_time_step", 		&lattice.adaptive_time_step);

	getDoubleProperty(&cfg,  "min_time_step", 			&lattice.min_time_step);

	config_destroy(&cfg);

	printf("\nLattice parameters:");
	printf("\n-------------------\n");
	printf("lattice_points_x         = %d\n", 	lattice.lattice_points_x);
	printf("lattice_points_y         = %d\n", 	lattice.lattice_points_y);
	printf("lattice_points_eta       = %d\n", 	lattice.lattice_points_eta);

	printf("max_number_of_time_steps = %d\n", 	lattice.max_number_of_time_steps);

	printf("lattice_spacing_x        = %.3g\n", lattice.lattice_spacing_x);
	printf("lattice_spacing_y        = %.3g\n", lattice.lattice_spacing_y);
	printf("lattice_spacing_eta      = %.3g\n", lattice.lattice_spacing_eta);

	printf("fixed_time_step          = %.3g\n", lattice.fixed_time_step);

	printf("output_period            = %d\n", 	lattice.output_period);
	printf("adaptive_time_step       = %d\n", 	lattice.adaptive_time_step);

	printf("min_time_step            = %.2e\n", lattice.min_time_step);
	printf("\n");

	return lattice;
}


initial_condition_parameters load_initial_condition_parameters()
{
	config_t cfg;
	config_init(&cfg);

	char fname[255] = "parameters/ic.properties";
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
	printf("\n");

	return initial;
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

	getIntegerProperty(&cfg,"run_hydro",				&hydro.run_hydro);
	getDoubleProperty(&cfg, "tau_initial", 				&hydro.tau_initial);
	getDoubleProperty(&cfg, "plpt_ratio_initial",		&hydro.plpt_ratio_initial);
	getIntegerProperty(&cfg,"temperature_etas", 		&hydro.temperature_etas);
	getDoubleProperty(&cfg, "constant_etas", 			&hydro.constant_etas);
	getDoubleProperty(&cfg, "freezeout_temperature_GeV",&hydro.freezeout_temperature_GeV);
	getDoubleProperty(&cfg, "flux_limiter",				&hydro.flux_limiter);
	getDoubleProperty(&cfg, "energy_min",				&hydro.energy_min);
	getDoubleProperty(&cfg, "pressure_min",				&hydro.pressure_min);
	getIntegerProperty(&cfg,"regulation_scheme", 		&hydro.regulation_scheme);
	getIntegerProperty(&cfg,"reprojection", 			&hydro.reprojection);
	getDoubleProperty(&cfg, "rho_max",					&hydro.rho_max);
	getDoubleProperty(&cfg, "xi0",						&hydro.xi0);

	config_destroy(&cfg);

	printf("Hydro parameters:");
	printf("\n-----------------\n");
	printf("run_hydro                 = %d\n",		hydro.run_hydro);
	printf("tau_initial               = %.3g\n",	hydro.tau_initial);
	printf("plpt_ratio_initial        = %.2g\n", 	hydro.plpt_ratio_initial);
	printf("temperature_etas          = %d\n", 		hydro.temperature_etas);
	printf("constant_etas             = %.3g\n", 	hydro.constant_etas);
	printf("freezeout_temperature_GeV = %.3f\n", 	hydro.freezeout_temperature_GeV);
	printf("flux_limiter              = %.3g\n", 	hydro.flux_limiter);
	printf("energy_min                = %.1e\n", 	hydro.energy_min);
	printf("pressure_min              = %.1e\n", 	hydro.pressure_min);
	printf("regulation_scheme         = %d\n", 		hydro.regulation_scheme);
	printf("reprojection              = %d\n", 		hydro.reprojection);
	printf("rho_max                   = %.2f\n", 	hydro.rho_max);
	printf("xi0                       = %.2f\n", 	hydro.xi0);
	printf("\n");

	if(fabs(hydro.flux_limiter - 1.5) > 0.5)
	{
		printf("load_hydro_parameters error: flux_limiter = %.3g is out of bounds\n", hydro.flux_limiter);
		exit(-1);
	}

	return hydro;
}



