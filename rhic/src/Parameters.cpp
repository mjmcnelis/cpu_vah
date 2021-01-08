#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include "../include/Macros.h"
#include "../include/Parameters.h"


double round_two_decimals(double number)
{
	return ceil((number * 100.) - 0.4999999) / 100.;
}

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



double compute_conformal_prefactor()
{
	double colors = 3.;
	double flavors = 3.;

	return M_PI * M_PI * (2. * (colors * colors  - 1.)  +  3.5 * colors * flavors) / 30.;

	// return 14.004;	// temporary (for freezeout surface test)
}


random_model_parameters load_random_model_parameters(int sample)
{
	printf("\nLoading model parameters from python/model_parameters/model_parameters_%d.dat\n\n", sample);

	FILE * parameter_file;
	char fname[255];
	sprintf(fname, "python/model_parameters/model_parameters_%d.dat", sample);
  	parameter_file = fopen(fname, "r");

  	if(parameter_file == NULL)
  	{
  		printf("load_random_model_parameters error: could not open model_parameters_%d.dat\n", sample);
  		exit(-1);
  	}

  	double b, N, p, w, dmin, sigmak;
  	double Tsw, Tk, etask, aL, aH, zetasN, Tp, wz, lambdaz;

  	fscanf(parameter_file, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf", &b, &N, &p, &w, &dmin, &sigmak, &Tsw, &Tk, &etask, &aL, &aH, &zetasN, &Tp, &wz, &lambdaz);
	fclose(parameter_file);

	printf("b       = %.3e\n", b);
	printf("N       = %.3e\n", N);
	printf("p       = %.3e\n", p);
	printf("w       = %.3e\n", w);
	printf("dmin    = %.3e\n", dmin);
	printf("sigmak  = %.3e\n", sigmak);
	printf("");
	printf("Tsw     = %.3e\n", Tsw);
	printf("Tk      = %.3e\n", Tk);
	printf("etask   = %.3e\n", etask);
	printf("aL      = %.3e\n", aL);
	printf("aH      = %.3e\n", aH);
	printf("zetasN  = %.3e\n", zetasN);
	printf("Tp      = %.3e\n", Tp);
	printf("wz      = %.3e\n", wz);
	printf("lambdaz = %.3e\n", lambdaz);

	random_model_parameters random;

	random.impact_parameter 				= b;
	random.trento_normalization_GeV 		= N;
	random.trento_geometric_parameter 		= p;
	random.trento_nucleon_width 			= w;
	random.trento_min_nucleon_distance 		= dmin;
	random.trento_gamma_standard_deviation	= sigmak;

	random.freezeout_temperature_GeV		= Tsw;
	random.etas_Tk_GeV						= Tk;
	random.etas_etask 						= etask;
	random.etas_aL							= aL;
	random.etas_aH							= aH;
	random.zetas_normalization_factor		= zetasN;
	random.zetas_peak_temperature_GeV		= Tp;
	random.zetas_width_GeV					= wz;
	random.zetas_skew						= lambdaz;

	return random;
}


hydro_parameters load_hydro_parameters(bool sample_parameters, random_model_parameters random)
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
		hydro.etas_min = atof(line.c_str());

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
		hydro.freezeout_finder_period = atoi(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		hydro.flux_limiter = atof(line.c_str());

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
		std::cerr << "No configuration file %s found for hydro parameters\n";
	}

	hydro.conformal_eos_prefactor = compute_conformal_prefactor();

#ifdef RANDOM_MODEL_PARAMETERS    				// overwrite hydro model parameters
	if(sample_parameters)
	{
		hydro.etas_aL 						= random.etas_aL;
		hydro.etas_aH 						= random.etas_aH;
		hydro.etas_Tk_GeV 					= random.etas_Tk_GeV;
		hydro.etas_etask 					= random.etas_etask;
		hydro.zetas_normalization_factor 	= random.zetas_normalization_factor;
		hydro.zetas_peak_temperature_GeV 	= random.zetas_peak_temperature_GeV;
		hydro.zetas_width_GeV 				= random.zetas_width_GeV;
		hydro.zetas_skew 					= random.zetas_skew;
		hydro.freezeout_temperature_GeV 	= random.freezeout_temperature_GeV;
	}
#endif

#ifdef PRINT_PARAMETERS
	printf("\nHydro parameters:");
	printf("\n-----------------\n");
	printf("run_hydro                  = %d\n",		hydro.run_hydro);
	printf("tau_initial                = %.3g\n",	hydro.tau_initial);
	printf("plpt_ratio_initial         = %.2g\n", 	hydro.plpt_ratio_initial);
	printf("kinetic_theory_model       = %d\n", 	hydro.kinetic_theory_model);
	printf("conformal_eos_prefactor    = %.6g\n", 	hydro.conformal_eos_prefactor);
	printf("temperature_etas           = %d\n", 	hydro.temperature_etas);
	printf("constant_etas              = %.3g\n", 	hydro.constant_etas);
	printf("etas_min                   = %.3g\n", 	hydro.etas_min);
	printf("etas_aL                    = %.3g\n", 	hydro.etas_aL);
	printf("etas_aH                    = %.3g\n", 	hydro.etas_aH);
	printf("etas_Tk_GeV                = %.3g\n", 	hydro.etas_Tk_GeV);
	printf("etas_etask                 = %.3g\n", 	hydro.etas_etask);
	printf("zetas_normalization_factor = %.3g\n", 	hydro.zetas_normalization_factor);
	printf("zetas_peak_temperature_GeV = %.3g\n", 	hydro.zetas_peak_temperature_GeV);
	printf("zetas_width_GeV            = %.3g\n", 	hydro.zetas_width_GeV);
	printf("zetas_skew                 = %.3g\n", 	hydro.zetas_skew);
	printf("freezeout_temperature_GeV  = %.3f\n", 	hydro.freezeout_temperature_GeV);
	printf("freezeout_finder_period    = %d\n", 	hydro.freezeout_finder_period);
	printf("flux_limiter               = %.3g\n", 	hydro.flux_limiter);
	printf("energy_min                 = %.1e\n", 	hydro.energy_min);
	printf("pressure_min               = %.1e\n", 	hydro.pressure_min);
	printf("regulation_scheme          = %d\n", 	hydro.regulation_scheme);
	printf("rho_max                    = %.2f\n", 	hydro.rho_max);
	printf("\n");
#endif

	#ifdef LATTICE_QCD
	double hotqcd_e_min = 0.00175;			// min energy density in fm^-4 rounded up to 3 SFs (see notebook in eos/hotqcd_smash)

	if(hydro.energy_min < hotqcd_e_min)
	{
		printf("load_hydro_parameters error: energy_min = %.3e is smaller than minimum energy density = %3e in hotqcd eos table\n", hydro.energy_min, hotqcd_e_min);
		exit(-1);
	}
	#endif

	if(hydro.tau_initial == 0)
	{
		printf("load_hydro_parameters error: set tau_initial > 0\n");
		exit(-1);
	}

	if(fabs(hydro.flux_limiter - 1.5) > 0.5)
	{
		printf("load_hydro_parameters error: flux_limiter = %.3g is not between [1,2]\n", hydro.flux_limiter);
		exit(-1);
	}

	return hydro;
}


lattice_parameters load_lattice_parameters(hydro_parameters hydro, initial_condition_parameters initial, bool sample_parameters, int sample)
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
		lattice.resolve_nucleons = atoi(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		lattice.fit_rapidity_plateau = atoi(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		lattice.training_grid = atoi(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		lattice.train_coarse_factor = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		lattice.auto_grid = atoi(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		lattice.sigma_factor = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		lattice.buffer = atof(line.c_str());

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
		std::cerr << "No configuration file %s found for hydro parameters\n";
	}


	if(lattice.resolve_nucleons)
	{
		double Lx = lattice.lattice_spacing_x * (lattice.lattice_points_x - 1);
		double Ly = lattice.lattice_spacing_y * (lattice.lattice_points_y - 1);

		double w = initial.trento_nucleon_width;

		double dx = round_two_decimals(w / 5.);
		double dy = round_two_decimals(w / 5.);

		lattice.lattice_spacing_x = dx;
		lattice.lattice_spacing_y = dy;

		int nx = (int)ceil(Lx / dx)  +  (1 + (int)ceil(Lx / dx)) % 2;
		int ny = (int)ceil(Ly / dy)  +  (1 + (int)ceil(Ly / dy)) % 2;

		if((nx - 1) * dx < Lx)
		{
			nx += 2;
		}

		if((ny - 1) * dy < Ly)
		{
			ny += 2;
		}

		lattice.lattice_points_x = nx;
		lattice.lattice_points_y = ny;
	}

#ifndef BOOST_INVARIANT
	if(lattice.fit_rapidity_plateau)
	{
		double eta_flat = initial.rapidity_flat;
		double eta_sigma = sqrt(initial.rapidity_variance);

		double Lz = eta_flat  +  10. * eta_sigma;			 // hard-coded to cover 5 eta_sigmas
		double dz = round_two_decimals(eta_sigma / 5.);

		lattice.lattice_spacing_eta = dz;

		int nz = (int)ceil(Lz / dz)  +  (1 + (int)ceil(Lz / dz)) % 2;

		if((nz - 1) * dz < Lz)
		{
			nz += 2;
		}
		lattice.lattice_points_eta = nz;
	}
#endif

	if(sample_parameters)
	{
		if(lattice.training_grid)
		{
			printf("\nRun hydro simulation on the training grid\n\n");

			double Lx = 30;
			double Ly = 30;
			double Lz = lattice.lattice_spacing_eta * (lattice.lattice_points_eta - 1);

			if(lattice.train_coarse_factor < 1)
			{
				printf("load_lattice_parameters error: for traning grid, set train_coarse_factor >= 1\n");
				exit(-1);
			}

		#ifndef BOOST_INVARIANT
			if(!lattice.fit_rapidity_plateau)
			{
				printf("load_lattice_parameters error: for training grid in 3+1d, set fit_rapidity_plateau = 1\n");
				exit(-1);
			}
		#endif

			double dx = lattice.lattice_spacing_x   * lattice.train_coarse_factor;
			double dy = lattice.lattice_spacing_y   * lattice.train_coarse_factor;
			double dz = lattice.lattice_spacing_eta * lattice.train_coarse_factor;

			lattice.lattice_spacing_x   = dx;
			lattice.lattice_spacing_y   = dy;
			lattice.lattice_spacing_eta = dz;

			int nx = (int)ceil(Lx / dx)  +  (1 + (int)ceil(Lx / dx)) % 2;
			int ny = (int)ceil(Ly / dy)  +  (1 + (int)ceil(Ly / dy)) % 2;
			int nz = (int)ceil(Lz / dz)  +  (1 + (int)ceil(Lz / dz)) % 2;

			if((nx - 1) * dx < Lx)
			{
				nx += 2;
			}

			if((ny - 1) * dy < Ly)
			{
				ny += 2;
			}

			if((nz - 1) * dz < Lz)
			{
				nz += 2;
			}

			lattice.lattice_points_x   = nx;
			lattice.lattice_points_y   = ny;
		#ifndef BOOST_INVARIANT
			lattice.lattice_points_eta = nz;
		#endif

		}
		else if(lattice.auto_grid)
		{
			printf("\nRun hydro simulation on the auto grid\n");

			printf("\nLoading predicted fireball size from python/fireball_size_predictions/fireball_size_%d.dat\n\n", sample);

			FILE * fireball_size;
			char fname[255];
			sprintf(fname, "python/fireball_size_predictions/fireball_size_%d.dat", sample);
		  	fireball_size = fopen(fname, "r");

		  	if(fireball_size == NULL)
		  	{
		  		printf("load_lattice_parameters error: could not open fireball_size_%d.dat. Need to run regression model\n", sample);
		  		exit(-1);
		  	}

		  	double fireball_radius_mean;
			double fireball_radius_std;

	  		fscanf(fireball_size, "%lf\t%lf", &fireball_radius_mean, &fireball_radius_std);
			fclose(fireball_size);

			printf("Predicted mean radius = %.2f fm\n", fireball_radius_mean);
			printf("Predicted std radius  = %.2f fm\n", fireball_radius_std);

			double sigma_factor = lattice.sigma_factor;
			double buffer 		= lattice.buffer;

			double L = 2. * (fireball_radius_mean  +  sigma_factor * fireball_radius_std  +  buffer);

			printf("\nAuto transverse grid size length L = %.3f fm\n\n", L);

			double dx = lattice.lattice_spacing_x;
			double dy = lattice.lattice_spacing_y;

			// ensure odd number of points
			int nx = (int)ceil(L / dx)  +  (1 + (int)ceil(L / dx)) % 2;
			int ny = (int)ceil(L / dy)  +  (1 + (int)ceil(L / dy)) % 2;

			if(((nx - 1) * dx < L) || ((ny - 1) * dy < L))
			{
				nx += 2;
				ny += 2;
			}

			lattice.lattice_points_x = nx;
			lattice.lattice_points_y = ny;

			printf("Reconfiguring transverse grid size length to L ~ %.3f fm\n\n", (nx - 1) * dx);
		}
	}

#ifdef BOOST_INVARIANT
	if(lattice.lattice_points_eta > 1)
	{
		printf("load_lattice_parameters: BOOST_INVARIANT is defined but you set lattice_points_eta = %d. Setting lattice_points_eta = 1\n\n", lattice.lattice_points_eta);
		lattice.lattice_points_eta = 1;		// automatic default
	}
#endif

	lattice.min_time_step = hydro.tau_initial / 20.;			// min time step ~ 20x smaller than t0

	double dt = lattice.fixed_time_step;						// dt = default starting time step is fixed_time_step

	if(!time_step_within_CFL_bound(dt, lattice))
	{
		printf("load_lattice_parameters flag: fixed time step dt = %.2g greater than fixed CFL bound: flooring fixed_time_step...\n", dt);
		double dx = lattice.lattice_spacing_x;
		double dy = lattice.lattice_spacing_y;
		double dz = lattice.lattice_spacing_eta;

		lattice.fixed_time_step = 0.125 * fmin(dx, fmin(dy, dz));
	}

	if(lattice.adaptive_time_step)
	{
		dt = lattice.min_time_step;								// if adaptive, starting time step is min_time_step
	}

	if(!time_step_within_CFL_bound(dt, lattice))
	{
		printf("load_lattice_parameters flag: starting time step dt = %.2g greater than fixed CFL bound: flooring min_time_step...\n", dt);
		double dx = lattice.lattice_spacing_x;
		double dy = lattice.lattice_spacing_y;
		double dz = lattice.lattice_spacing_eta;

		lattice.min_time_step = 0.125 * fmin(dx, fmin(dy, dz));
	}


#ifdef PRINT_PARAMETERS
	printf("Lattice parameters:");
	printf("\n-------------------\n");
	printf("lattice_points_x    = %d\n", 	lattice.lattice_points_x);
	printf("lattice_points_y    = %d\n", 	lattice.lattice_points_y);
	printf("lattice_points_eta  = %d\n", 	lattice.lattice_points_eta);
	printf("lattice_spacing_x   = %.3g\n",	lattice.lattice_spacing_x);
	printf("lattice_spacing_y   = %.3g\n", 	lattice.lattice_spacing_y);
	printf("lattice_spacing_eta = %.3g\n", 	lattice.lattice_spacing_eta);
	printf("resolve_nucleons    = %d\n", 	lattice.resolve_nucleons);
	printf("fit_rapidity_plateau= %d\n", 	lattice.fit_rapidity_plateau);
	printf("training_grid       = %d\n", 	lattice.training_grid);
	printf("train_coarse_factor = %.2g\n", 	lattice.train_coarse_factor);
	printf("auto_grid           = %d\n", 	lattice.auto_grid);
	printf("sigma_factor        = %.3g\n", 	lattice.sigma_factor);
	printf("buffer              = %.3g\n", 	lattice.buffer);
	printf("max_time_steps      = %d\n", 	lattice.max_time_steps);
	printf("output_interval     = %.2f\n", 	lattice.output_interval);
	printf("fixed_time_step     = %.3g\n", 	lattice.fixed_time_step);
	printf("adaptive_time_step  = %d\n", 	lattice.adaptive_time_step);
	printf("min_time_step       = %.2e\n", 	lattice.min_time_step);
	printf("delta_0             = %.3g\n", 	lattice.delta_0);
	printf("alpha               = %.3g\n", 	lattice.alpha);
	printf("\n");
#endif

	return lattice;
}


initial_condition_parameters load_initial_condition_parameters(bool sample_parameters, random_model_parameters random)
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
		initial.impact_parameter = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initial.rapidity_variance = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initial.rapidity_flat = atof(line.c_str());

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
		std::cerr << "No configuration file %s found for hydro parameters\n";
	}

#ifdef RANDOM_MODEL_PARAMETERS    				// overwrite hydro model parameters
	if(sample_parameters)
	{
		initial.impact_parameter				= random.impact_parameter;
		initial.trento_normalization_GeV 		= random.trento_normalization_GeV;
		initial.trento_nucleon_width 			= random.trento_nucleon_width;
		initial.trento_min_nucleon_distance 	= random.trento_min_nucleon_distance;
		initial.trento_geometric_parameter 		= random.trento_geometric_parameter;
		initial.trento_gamma_standard_deviation = random.trento_gamma_standard_deviation;
	}
#endif

#ifdef PRINT_PARAMETERS
	printf("Initial condition parameters:");
	printf("\n-----------------------------\n");
	printf("initial_condition_type          = %d\n", 	initial.initial_condition_type);
	printf("nucleus_A                       = %d\n", 	initial.nucleus_A);
	printf("nucleus_B                       = %d\n", 	initial.nucleus_B);
	printf("initial_central_temperature_GeV = %.3g\n",	initial.initialCentralTemperatureGeV);
	printf("impact_parameter                = %.2f\n", 	initial.impact_parameter);
	printf("rapidityVariance                = %.3g\n", 	initial.rapidity_variance);
	printf("rapidityMean                    = %.3g\n", 	initial.rapidity_flat);
	printf("q_gubser                        = %.2f\n", 	initial.q_gubser);
	printf("trento_normalization_GeV        = %.2f\n", 	initial.trento_normalization_GeV);
	printf("trento_nucleon_width            = %.2f\n", 	initial.trento_nucleon_width);
	printf("trento_min_nucleon_distance     = %.2f\n", 	initial.trento_min_nucleon_distance);
	printf("trento_geometric_parameter      = %.2f\n", 	initial.trento_geometric_parameter);
	printf("trento_gamma_standard_deviation = %.2f\n", 	initial.trento_gamma_standard_deviation);
	printf("trento_average_over_events      = %d\n", 	initial.trento_average_over_events);
	printf("trento_number_of_average_events = %d\n", 	initial.trento_number_of_average_events);
	printf("trento_fixed_seed               = %d\n", 	initial.trento_fixed_seed);
	printf("\n");
#endif

	return initial;
}












