
#ifndef PARAMETERS_H_
#define PARAMETERS_H_

typedef struct
{
	int lattice_points_x;
	int lattice_points_y;
	int lattice_points_eta;
	int max_number_of_time_steps;

	double lattice_spacing_x;
	double lattice_spacing_y;
	double lattice_spacing_eta;
	double fixed_time_step;

	int output_period;
	int adaptive_time_step;

	double min_time_step;

} lattice_parameters;


typedef struct
{
	int initialConditionType;
	int numberOfNucleonsPerNuclei;			// can only do A = B collision

	double initialCentralTemperatureGeV;	// initial central temperature for a central collision
	double scatteringCrossSectionNN;
	double impactParameter;
	double fractionOfBinaryCollisions;

	// longitudinal energy density profile parameters
	double rapidityVariance; 	// \sigma^{2}_{\eta}
	double rapidityMean; 		// flat region around \ets_s = 0

} initial_condition_parameters;


typedef struct
{
	int run_hydro;
	double tau_initial;
	int temperature_etas;
	double constant_etas;
	double freezeout_temperature_GeV;
	double flux_limiter;
	double energy_min;
	double plpt_ratio_initial;

} hydro_parameters;


lattice_parameters load_lattice_parameters();
initial_condition_parameters load_initial_condition_parameters();
hydro_parameters load_hydro_parameters();

#endif




