#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "../include/Parameters.h"
using namespace std;

void getIntegerProperty(config_t * cfg, const char * propName, int * propValue)
{
	config_lookup_int(cfg, propName, propValue);
}


void getDoubleProperty(config_t * cfg, const char * propName, double * propValue)
{
	config_lookup_float(cfg, propName, propValue);
}


void load_lattice_parameters(config_t *cfg, void * params)
{
	char fname[255];
	sprintf(fname, "parameters/lattice.properties");

	if(!config_read_file(cfg, fname))
	{
		fprintf(stderr, "No configuration file %s found for lattice parameters - %s.\n", fname, config_error_text(cfg));
		exit(-1);
	}

	int lattice_points_x;
	int lattice_points_y;
	int lattice_points_eta;

	int max_number_of_time_steps;

	double lattice_spacing_x;
	double lattice_spacing_y;
	double lattice_spacing_eta;
	double fixed_time_step;

	int adaptive_time_step;
	double min_time_step;

	// get the lattice parameters
	getIntegerProperty(cfg, "lattice_points_x", 		&lattice_points_x);
	getIntegerProperty(cfg, "lattice_points_y", 		&lattice_points_y);
	getIntegerProperty(cfg, "lattice_points_eta", 		&lattice_points_eta);
	getIntegerProperty(cfg, "max_number_of_time_steps", &max_number_of_time_steps);
	getDoubleProperty(cfg,  "lattice_spacing_x", 		&lattice_spacing_x);
	getDoubleProperty(cfg,  "lattice_spacing_y", 		&lattice_spacing_y);
	getDoubleProperty(cfg,  "lattice_spacing_eta", 		&lattice_spacing_eta);
	getDoubleProperty(cfg,  "fixed_time_step", 			&fixed_time_step);
	getIntegerProperty(cfg, "adaptive_time_step", 		&adaptive_time_step);
	getDoubleProperty(cfg,  "min_time_step", 			&min_time_step);

#ifdef PRINT_PARAMETERS
	printf("\nLattice parameters:");
	printf("\n-------------------\n");
	printf("lattice_points_x         = %d\n", 	lattice_points_x);
	printf("lattice_points_y         = %d\n", 	lattice_points_y);
	printf("lattice_points_eta       = %d\n", 	lattice_points_eta);
	printf("max_number_of_time_steps = %d\n", 	max_number_of_time_steps);
	printf("lattice_spacing_x        = %.3g\n", lattice_spacing_x);
	printf("lattice_spacing_y        = %.3g\n", lattice_spacing_y);
	printf("lattice_spacing_eta      = %.3g\n", lattice_spacing_eta);
	printf("fixed_time_step          = %.3g\n", fixed_time_step);
	printf("adaptive_time_step       = %d\n", 	adaptive_time_step);
	printf("min_time_step            = %.3g\n", min_time_step);
	printf("\n");
#endif

	struct LatticeParameters * lattice = (struct LatticeParameters *) params;
	lattice->lattice_points_x 			= lattice_points_x;
	lattice->lattice_points_y 			= lattice_points_y;
	lattice->lattice_points_eta			= lattice_points_eta;
	lattice->max_number_of_time_steps 	= max_number_of_time_steps;
	lattice->lattice_spacing_x 			= lattice_spacing_x;
	lattice->lattice_spacing_y 			= lattice_spacing_y;
	lattice->lattice_spacing_eta 		= lattice_spacing_eta;
	lattice->fixed_time_step 			= fixed_time_step;
	lattice->adaptive_time_step 		= adaptive_time_step;
	lattice->min_time_step 				= min_time_step;
}


void loadInitialConditionParameters(config_t *cfg, void * params)
{
	char fname[255];
	sprintf(fname, "parameters/ic.properties");

	if(!config_read_file(cfg, fname))
	{
		fprintf(stderr, "No configuration file  %s found for initial condition parameters - %s.\n", fname, config_error_text(cfg));
		exit(-1);
	}

	int initialConditionType, numberOfNucleonsPerNuclei;
	double initialCentralTemperatureGeV, scatteringCrossSectionNN, impactParameter, fractionOfBinaryCollisions, rapidityVariance, rapidityMean;

	// get initial condition parameters
	getIntegerProperty(cfg, "initialConditionType", 		&initialConditionType);
	getIntegerProperty(cfg, "numberOfNucleonsPerNuclei", 	&numberOfNucleonsPerNuclei);
	getDoubleProperty(cfg, 	"initialCentralTemperatureGeV",	&initialCentralTemperatureGeV);
	getDoubleProperty(cfg, 	"scatteringCrossSectionNN", 	&scatteringCrossSectionNN);
	getDoubleProperty(cfg, 	"impactParameter", 				&impactParameter);
	getDoubleProperty(cfg,	"fractionOfBinaryCollisions", 	&fractionOfBinaryCollisions);
	getDoubleProperty(cfg, 	"rapidityVariance", 			&rapidityVariance);
	getDoubleProperty(cfg, 	"rapidityMean", 				&rapidityMean);

#ifdef PRINT_PARAMETERS
	printf("Initial condition parameters:");
	printf("\n-----------------------------\n");
	printf("initialConditionType         = %d\n", initialConditionType);
	printf("numberOfNucleonsPerNuclei    = %d\n", numberOfNucleonsPerNuclei);
	printf("initialCentralTemperatureGeV = %.3g\n", initialCentralTemperatureGeV);
	printf("scatteringCrossSectionNN     = %.3g\n", scatteringCrossSectionNN);
	printf("impactParameter              = %.3g\n", impactParameter);
	printf("fractionOfBinaryCollisions   = %.3g\n", fractionOfBinaryCollisions);
	printf("rapidityVariance             = %.3g\n", rapidityVariance);
	printf("rapidityMean                 = %.3g\n", rapidityMean);
	printf("\n");
#endif

	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) params;
	initCond->initialConditionType 			= initialConditionType;
	initCond->numberOfNucleonsPerNuclei 	= numberOfNucleonsPerNuclei;
	initCond->initialCentralTemperatureGeV 	= initialCentralTemperatureGeV;
	initCond->scatteringCrossSectionNN 		= scatteringCrossSectionNN;
	initCond->impactParameter 				= impactParameter;
	initCond->fractionOfBinaryCollisions 	= fractionOfBinaryCollisions;
	initCond->rapidityVariance 				= rapidityVariance;
	initCond->rapidityMean 					= rapidityMean;
}


void loadHydroParameters(config_t *cfg, void * params)
{
	char fname[255];
	sprintf(fname, "parameters/hydro.properties");

	if(!config_read_file(cfg, fname))
	{
		fprintf(stderr, "No configuration file  %s found for hydrodynamic parameters - %s.\n", fname, config_error_text(cfg));
		exit(-1);
	}

	double tau_initial;
	double shear_viscosity;
	double freezeout_temperature_GeV;
	double flux_limiter;

	// get hydro parameters
	getDoubleProperty(cfg, "tau_initial", 				&tau_initial);
	getDoubleProperty(cfg, "shear_viscosity", 			&shear_viscosity);
	getDoubleProperty(cfg, "freezeout_temperature_GeV",	&freezeout_temperature_GeV);
	getDoubleProperty(cfg, "flux_limiter",				&flux_limiter);

#ifdef PRINT_PARAMETERS
	printf("Hydro parameters:");
	printf("\n-----------------\n");
	printf("tau_initial             = %.3f\n", tau_initial);
	printf("shear_viscosity         = %.3f\n", shear_viscosity);
	printf("freezeout_temperature_GeV = %.3f\n", freezeout_temperature_GeV);
	printf("flux_limiter            = %.3f\n", flux_limiter);
	printf("\n");
#endif

	// set hydro struct
	struct HydroParameters * hydro = (struct HydroParameters *) params;
	hydro->tau_initial				= tau_initial;
	hydro->shear_viscosity			= shear_viscosity;
	hydro->freezeout_temperature_GeV	= freezeout_temperature_GeV;
	hydro->flux_limiter				= flux_limiter;
}


hydro_parameters load_hydro_parameters(config_t *cfg)
{
	char filename[255] = "parameters/hydro.properties";
	if(!config_read_file(cfg, filename))
	{
		fprintf(stderr, "No configuration file  %s found for hydrodynamic parameters - %s.\n", filename, config_error_text(cfg));
		exit(-1);
	}

	hydro_parameters hydro;

	getDoubleProperty(cfg, "tau_initial", 				&hydro.tau_initial);
	getDoubleProperty(cfg, "shear_viscosity", 			&hydro.shear_viscosity);
	getDoubleProperty(cfg, "freezeout_temperature_GeV",	&hydro.freezeout_temperature_GeV);
	getDoubleProperty(cfg, "flux_limiter",				&hydro.flux_limiter);

#ifdef PRINT_PARAMETERS
	printf("Hydro parameters:");
	printf("\n-----------------\n");
	printf("tau_initial               = %.3f\n", hydro.tau_initial);
	printf("shear_viscosity           = %.3f\n", hydro.shear_viscosity);
	printf("freezeout_temperature_GeV = %.3f\n", hydro.freezeout_temperature_GeV);
	printf("flux_limiter              = %.3f\n", hydro.flux_limiter);
	printf("\n");
#endif

	return hydro;
}



