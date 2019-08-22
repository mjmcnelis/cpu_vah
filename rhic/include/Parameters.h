
#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <libconfig.h>

#define PRINT_PARAMETERS

void getIntegerProperty(config_t * cfg, const char * propName, int * propValue);

void getDoubleProperty(config_t * cfg, const char * propName, double * propValue);


struct LatticeParameters
{
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
};


struct InitialConditionParameters
{
	int initialConditionType;

	int numberOfNucleonsPerNuclei;			// this is an A = B collision

	double initialCentralTemperatureGeV;	// initial central temperature for a central A = B collision (b = 0)
	double scatteringCrossSectionNN;
	double impactParameter;
	double fractionOfBinaryCollisions;

	// longitudinal energy density profile parameters
	double rapidityVariance; 	// \sigma^{2}_{\eta}
	double rapidityMean; 		// flat region around \ets_s = 0
};


struct HydroParameters
{
	double tau_initial;
	double shear_viscosity;
	double freezeout_temperature_GeV;
	double flux_limiter;
};


typedef struct 
{
	double tau_initial;
	double shear_viscosity;
	double freezeout_temperature_GeV;
	double flux_limiter;

} hydro_parameters;


void load_lattice_parameters(config_t *cfg, void * params);

void loadInitialConditionParameters(config_t *cfg, void * params);

void loadHydroParameters(config_t *cfg, void * params);

hydro_parameters load_hydro_parameters(config_t *cfg);

#endif




