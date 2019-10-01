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
	char fname[255] = "parameters/hydro.properties";

	hydro_parameters hydro;

	double quark_flavors;

	std::ifstream cFile(fname);
	if(cFile.is_open())
	{
		std::string line;

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		auto delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		hydro.run_hydro = atoi(line.c_str());

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
		hydro.etas_min = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		hydro.etas_slope = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		hydro.constant_etas = atof(line.c_str());

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

	printf("\nHydro parameters:");
	printf("\n-----------------\n");
	printf("run_hydro                  = %d\n",		hydro.run_hydro);
	printf("tau_initial                = %.3g\n",	hydro.tau_initial);
	printf("plpt_ratio_initial         = %.2g\n", 	hydro.plpt_ratio_initial);
	printf("quark_flavors              = %.1f\n", 	quark_flavors);
	printf("conformal_eos_prefactor    = %.6g\n", 	hydro.conformal_eos_prefactor);
	printf("temperature_etas           = %d\n", 	hydro.temperature_etas);
	printf("etas_min                   = %.3g\n", 	hydro.etas_min);
	printf("etas_slope                 = %.3g\n", 	hydro.etas_slope);
	printf("constant_etas              = %.3g\n", 	hydro.constant_etas);
	printf("zetas_normalization_factor = %.3g\n", 	hydro.zetas_normalization_factor);
	printf("zetas_peak_temperature_GeV = %.3g\n", 	hydro.zetas_peak_temperature_GeV);
	printf("freezeout_temperature_GeV  = %.3f\n", 	hydro.freezeout_temperature_GeV);
	printf("flux_limiter               = %.3g\n", 	hydro.flux_limiter);
	printf("include_vorticity          = %d\n", 	hydro.include_vorticity);
	printf("energy_min                 = %.1e\n", 	hydro.energy_min);
	printf("pressure_min               = %.1e\n", 	hydro.pressure_min);
	printf("regulation_scheme          = %d\n", 	hydro.regulation_scheme);
	printf("rho_max                    = %.2f\n", 	hydro.rho_max);
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
	}
	else
	{
		std::cerr << "No configuration file  %s found for hydro parameters\n";
	}

	lattice.min_time_step = 2. * pow(10., round(log10(hydro.tau_initial)) - 2.);		// min time step ~ 20x smaller than t0

#ifdef BOOST_INVARIANT
	lattice.lattice_points_eta = 1;		// automatic default
#endif

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
	printf("delta_0                  = %.3g\n", lattice.delta_0);
	printf("alpha                    = %.3g\n", lattice.alpha);
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
		initial.initialConditionType = atoi(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initial.numberOfNucleonsPerNuclei = atoi(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initial.initialCentralTemperatureGeV = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initial.scatteringCrossSectionNN = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initial.impactParameter = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initial.fractionOfBinaryCollisions = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initial.rapidityVariance = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initial.rapidityMean = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initial.q_gubser = atof(line.c_str());
	}
	else
	{
		std::cerr << "No configuration file  %s found for hydro parameters\n";
	}

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





