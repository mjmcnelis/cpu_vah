
#ifndef PARAMETERS_H_
#define PARAMETERS_H_


typedef struct
{
	int run_hydro;

	double tau_initial;
	double plpt_ratio_initial;

	double conformal_eos_prefactor;	 // e = conformal_eos_prefactor * T^4

	int temperature_etas;
	double etas_min;
	double etas_slope;
	double constant_etas;

	double zetas_normalization_factor;
	double zetas_peak_temperature_GeV;

	double freezeout_temperature_GeV;
	double flux_limiter;
	int include_vorticity;

	double energy_min;
	double pressure_min;

	int regulation_scheme;

	double rho_max;

} hydro_parameters;


typedef struct
{
	int lattice_points_x;
	int lattice_points_y;
	int lattice_points_eta;
	double lattice_spacing_x;
	double lattice_spacing_y;
	double lattice_spacing_eta;

	int max_time_steps;
	double output_interval;

	double fixed_time_step;

	int adaptive_time_step;
	double min_time_step;		// minimum time step set for the program

	double delta_0;
	double alpha;

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

	double q_gubser;

} initial_condition_parameters;


hydro_parameters load_hydro_parameters();
lattice_parameters load_lattice_parameters(hydro_parameters hydro);
initial_condition_parameters load_initial_condition_parameters();

#endif




