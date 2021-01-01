
#ifndef PARAMETERS_H_
#define PARAMETERS_H_

typedef struct
{
	double impact_parameter;				// initial condition parameters
	double trento_normalization_GeV;
	double trento_geometric_parameter;
	double trento_nucleon_width;
	double trento_min_nucleon_distance;
	double trento_gamma_standard_deviation;

	double freezeout_temperature_GeV;		// hydro parameters
	double etas_Tk_GeV;
	double etas_etask;
	double etas_aL;
	double etas_aH;
	double zetas_normalization_factor;
	double zetas_peak_temperature_GeV;
	double zetas_width_GeV;
	double zetas_skew;

} random_model_parameters;


typedef struct
{
	int run_hydro;

	double tau_initial;
	double plpt_ratio_initial;

	int kinetic_theory_model;

	double conformal_eos_prefactor;	 // for e = conformal_eos_prefactor.T^4

	int temperature_etas;

	double constant_etas;
	double etas_min;

	double etas_aL;
	double etas_aH;
	double etas_etask;
	double etas_Tk_GeV;

	double zetas_normalization_factor;
	double zetas_peak_temperature_GeV;
	double zetas_width_GeV;
	double zetas_skew;

	double freezeout_temperature_GeV;
	int freezeout_finder_period;

	double flux_limiter;

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

	int resolve_nucleons;
	int fit_rapidity_plateau;

	int training_grid;
	double train_coarse_factor;
	int auto_grid;

	double sigma_factor;
	double buffer;

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
	int initial_condition_type;
	int nucleus_A;
	int nucleus_B;

	double initialCentralTemperatureGeV;
	double impact_parameter;
	double rapidity_variance;
	double rapidity_flat;

	double q_gubser;

	double trento_normalization_GeV;
	double trento_nucleon_width;
	double trento_min_nucleon_distance;
	double trento_geometric_parameter;
	double trento_gamma_standard_deviation;

	int trento_average_over_events;
	int trento_number_of_average_events;
	int trento_fixed_seed;

} initial_condition_parameters;


random_model_parameters load_random_model_parameters(int sample);
hydro_parameters load_hydro_parameters(bool sample_parameters, random_model_parameters random);
lattice_parameters load_lattice_parameters(hydro_parameters hydro, initial_condition_parameters initial, bool sample_parameters, int sample);
initial_condition_parameters load_initial_condition_parameters(bool sample_parameters, random_model_parameters random);

#endif




