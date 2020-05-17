#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <string>
#include "../include/Macros.h"
#include "../include/Parameters.h"
using namespace std;


bool time_step_within_CFL_bound(double dt, lattice_parameters lattice)
{
	double dx = lattice.lattice_spacing_x;
	double dy = lattice.lattice_spacing_y;
	double dn = lattice.lattice_spacing_eta;

	double dt_CFL = 0.125 * fmin(dx, fmin(dy, dn));

	if(dt <= dt_CFL)
	{
		return true;
	}
	else
	{
		return false;
	}
}



double compute_conformal_prefactor(double flavors)
{
	double colors = 3.;

	return M_PI * M_PI * (2. * (colors * colors  - 1.)  +  3.5 * colors * flavors) / 30.;
}


hydro_parameters load_hydro_parameters()
{
	char fname[255] = "parameters/hydro.properties";

	hydro_parameters hydro;

	double quark_flavors;

	std::ifstream cFile(fname);
	if(cFile.is_open())
	{
		std::string line;

		// can't I just make a function?...

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		auto delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		hydro.run_hydro = atoi(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		hydro.output = atoi(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		hydro.tau_initial = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		hydro.plpt_ratio_initial = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		hydro.kinetic_theory_model = atoi(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		quark_flavors = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		hydro.temperature_etas = atoi(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		hydro.constant_etas = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		hydro.etas_aL = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		hydro.etas_aH = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		hydro.etas_Tk_GeV = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		hydro.etas_etask = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		hydro.zetas_normalization_factor = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		hydro.zetas_peak_temperature_GeV = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		hydro.zetas_width_GeV = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		hydro.zetas_skew = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		hydro.freezeout_temperature_GeV = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		hydro.flux_limiter = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		hydro.include_vorticity = atoi(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		hydro.energy_min = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		hydro.pressure_min = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		hydro.regulation_scheme = atoi(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		hydro.rho_max = atof(line.c_str());
	}
	else
	{
		std::cerr << "No configuration file  %s found for hydro parameters\n";
	}

	hydro.conformal_eos_prefactor = compute_conformal_prefactor(quark_flavors);

	// printf("\nHydro parameters:");
	// printf("\n-----------------\n");
	// printf("run_hydro                  = %d\n",		hydro.run_hydro);
	// printf("output                     = %d\n",		hydro.output);
	// printf("tau_initial                = %.3g\n",	hydro.tau_initial);
	// printf("plpt_ratio_initial         = %.2g\n", 	hydro.plpt_ratio_initial);
	// printf("kinetic_theory_model       = %d\n", 	hydro.kinetic_theory_model);
	// printf("quark_flavors              = %.1f\n", 	quark_flavors);
	// printf("conformal_eos_prefactor    = %.6g\n", 	hydro.conformal_eos_prefactor);
	// printf("temperature_etas           = %d\n", 	hydro.temperature_etas);
	// printf("etas_min                   = %.3g\n", 	hydro.etas_min);
	// printf("etas_slope                 = %.3g\n", 	hydro.etas_slope);
	// printf("constant_etas              = %.3g\n", 	hydro.constant_etas);
	// printf("zetas_normalization_factor = %.3g\n", 	hydro.zetas_normalization_factor);
	// printf("zetas_peak_temperature_GeV = %.3g\n", 	hydro.zetas_peak_temperature_GeV);
	// printf("freezeout_temperature_GeV  = %.3f\n", 	hydro.freezeout_temperature_GeV);
	// printf("flux_limiter               = %.3g\n", 	hydro.flux_limiter);
	// printf("include_vorticity          = %d\n", 	hydro.include_vorticity);
	// printf("energy_min                 = %.1e\n", 	hydro.energy_min);
	// printf("pressure_min               = %.1e\n", 	hydro.pressure_min);
	// printf("regulation_scheme          = %d\n", 	hydro.regulation_scheme);
	// printf("rho_max                    = %.2f\n", 	hydro.rho_max);
	// printf("\n");

	#ifdef LATTICE_QCD
	double hotqcd_e_min = 0.00175;			// min energy density in fm^-4 rounded up to 3 SFs (see notebook in eos/hotqcd_smash)

	if(hydro.energy_min < hotqcd_e_min)
	{
		printf("load_hydro_parameters error: energy_min = %.3f is smaller than minimum energy density = %lf in hotqcd eos table\n", hydro.energy_min, hotqcd_e_min);
		exit(-1);
	}
	#endif

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
	char fname[255] = "parameters/lattice.properties";

	lattice_parameters lattice;

	std::ifstream cFile(fname);
	if(cFile.is_open())
	{
		std::string line;

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		auto delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		lattice.lattice_points_x = atoi(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		lattice.lattice_points_y = atoi(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		lattice.lattice_points_eta = atoi(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		lattice.lattice_spacing_x = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		lattice.lattice_spacing_y = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		lattice.lattice_spacing_eta = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		lattice.max_time_steps = atoi(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		lattice.output_interval = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		lattice.fixed_time_step = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		lattice.adaptive_time_step = atoi(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		lattice.delta_0 = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		lattice.alpha = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		lattice.tau_coarse_factor = atoi(line.c_str());
	}
	else
	{
		std::cerr << "No configuration file  %s found for hydro parameters\n";
	}

	//lattice.min_time_step = 5. * pow(10., round(log10(hydro.tau_initial)) - 2.);		// min time step ~ 20x smaller than t0

	lattice.min_time_step = hydro.tau_initial / 20.;

#ifdef BOOST_INVARIANT
	if(lattice.lattice_points_eta > 1)
	{
		printf("load_lattice_parameters: BOOST_INVARIANT is defined but you set lattice_spacing_eta = %d. Setting lattice_spacing_eta = 1\n", lattice.lattice_points_eta);
		lattice.lattice_points_eta = 1;		// automatic default
	}
#endif

	// printf("Lattice parameters:");
	// printf("\n-------------------\n");
	// printf("lattice_points_x    = %d\n", 	lattice.lattice_points_x);
	// printf("lattice_points_y    = %d\n", 	lattice.lattice_points_y);
	// printf("lattice_points_eta  = %d\n", 	lattice.lattice_points_eta);
	// printf("lattice_spacing_x   = %.3g\n",	lattice.lattice_spacing_x);
	// printf("lattice_spacing_y   = %.3g\n", 	lattice.lattice_spacing_y);
	// printf("lattice_spacing_eta = %.3g\n", 	lattice.lattice_spacing_eta);
	// printf("max_time_steps      = %d\n", 	lattice.max_time_steps);
	// printf("output_interval     = %.2f\n", 	lattice.output_interval);
	// printf("fixed_time_step     = %.3g\n", 	lattice.fixed_time_step);
	// printf("adaptive_time_step  = %d\n", 	lattice.adaptive_time_step);
	// printf("min_time_step       = %.2e\n", 	lattice.min_time_step);
	// printf("delta_0             = %.3g\n", 	lattice.delta_0);
	// printf("alpha               = %.3g\n", 	lattice.alpha);
	// printf("tau_coarse_factor   = %d\n", 	lattice.tau_coarse_factor);
	// printf("\n");

	//exit(-1);

	double dt = lattice.fixed_time_step;						// dt = default starting time step

	if(!time_step_within_CFL_bound(dt, lattice))
	{
		printf("load_lattice_parameters error: fixed time step dt = %.2g greater than strict CFL bound\n", dt);
		exit(-1);
	}

	if(lattice.adaptive_time_step)
	{
		dt = lattice.min_time_step;								// if adaptive, starting time step is min
	}

	if(!time_step_within_CFL_bound(dt, lattice))
	{
		printf("load_lattice_parameters error: starting time step dt = %.2g greater than strict CFL bound\n", dt);
		exit(-1);
	}

	return lattice;
}


initial_condition_parameters load_initial_condition_parameters()
{
	char fname[255] = "parameters/initial.properties";

	initial_condition_parameters initial;

	std::ifstream cFile(fname);
	if(cFile.is_open())
	{
		std::string line;

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		auto delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initial.initial_condition_type = atoi(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initial.nucleus_A = atoi(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initial.nucleus_B = atoi(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initial.initialCentralTemperatureGeV = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initial.impactParameter = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initial.rapidity_variance = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initial.rapidity_mean = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initial.q_gubser = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initial.trento_normalization_GeV = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initial.trento_nucleon_width = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initial.trento_min_nucleon_distance = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initial.trento_geometric_parameter = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initial.trento_gamma_standard_deviation = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initial.trento_average_over_events = atoi(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initial.trento_number_of_average_events = atoi(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initial.trento_fixed_seed = atoi(line.c_str());
	}
	else
	{
		std::cerr << "No configuration file  %s found for hydro parameters\n";
	}

	// printf("Initial condition parameters:");
	// printf("\n-----------------------------\n");
	// printf("initial_condition_type          = %d\n", 	initial.initial_condition_type);
	// printf("nucleus_A                       = %d\n", 	initial.nucleus_A);
	// printf("nucleus_B                       = %d\n", 	initial.nucleus_B);
	// printf("initial_central_temperature_GeV = %.3g\n",	initial.initialCentralTemperatureGeV);
	// printf("impactParameter                 = %.2f\n", 	initial.impactParameter);
	// printf("rapidityVariance                = %.3g\n", 	initial.rapidity_variance);
	// printf("rapidityMean                    = %.3g\n", 	initial.rapidity_mean);
	// printf("q_gubser                        = %.2f\n", 	initial.q_gubser);
	// printf("trento_normalization_GeV        = %.2f\n", 	initial.trento_normalization_GeV);
	// printf("trento_nucleon_width            = %.2f\n", 	initial.trento_nucleon_width);
	// printf("trento_min_nucleon_distance     = %.2f\n", 	initial.trento_min_nucleon_distance);
	// printf("trento_geometric_parameter      = %.2f\n", 	initial.trento_geometric_parameter);
	// printf("trento_gamma_standard_deviation = %.2f\n", 	initial.trento_gamma_standard_deviation);
	// printf("trento_average_over_events      = %d\n", 	initial.trento_average_over_events);
	// printf("trento_number_of_average_events = %d\n", 	initial.trento_number_of_average_events);
	// printf("trento_fixed_seed               = %d\n", 	initial.trento_fixed_seed);
	// printf("\n");

	return initial;
}












