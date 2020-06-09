#include <iostream>
#include <sstream>
#include <string>
#include "../include/HydroWrapper.h"
#include "../include/EnergyVector.h"
#include "../include/Parameters.h"
#include "../include/Print.h"
#include "../include/Macros.h"
#include "../include/Hydrodynamics.h"
#include "../include/FileIO.h"
//#include "../include/RuntimeParameters.h"		// Derek: what are the runtime parameters for?
using namespace std;

// this is a C++ wrapper for OSU hydro codes (cpu-vh, gpu-vh, cpu_vah, etc)
// originally written by Derek Everett 2018


HYDRO::HYDRO()
{

}


HYDRO::~HYDRO()
{

}


void HYDRO::load_initial_energy_density_vector(std::vector<double> energy_vector)
{
	// set initial_energy_density from a vector (useful for JETSCAPE)
    // note: argument should be [GeV/fm^3], later we convert it to [fm^-4]

	initial_energy_density_vector = energy_vector;
}


void HYDRO::store_freezeout_surface(freezeout_surface surface)
{
	printf("\nStoring freezeout surface in memory for particlization module...\n\n");

	long number_of_cells = surface.tau.size();

	if(number_of_cells == 0)
	{
		printf("HYDRO::save_freezeout_surface flag: freezeout surface is empty\n");
		return;
	}
	else
	{
		printf("Number of freezeout cells = %ld\n", number_of_cells);

		for(long i = 0; i < number_of_cells; i++)
		{
			tau.push_back((double)surface.tau[i]);
			x.push_back(  (double)surface.x[i]);
			y.push_back(  (double)surface.y[i]);
			eta.push_back((double)surface.eta[i]);

			dsigma_tau.push_back((double)surface.dsigma_tau[i]);
			dsigma_x.push_back(  (double)surface.dsigma_x[i]);
			dsigma_y.push_back(  (double)surface.dsigma_y[i]);
			dsigma_eta.push_back((double)surface.dsigma_eta[i]);

			ux.push_back((double)surface.ux[i]);
			uy.push_back((double)surface.uy[i]);
			un.push_back((double)surface.un[i]);

			E.push_back((double)surface.E[i]);
			T.push_back((double)surface.T[i]);
			P.push_back((double)surface.P[i]);

			pixx.push_back((double)surface.pixx[i]);
			pixy.push_back((double)surface.pixy[i]);
			pixn.push_back((double)surface.pixn[i]);
			piyy.push_back((double)surface.piyy[i]);
			piyn.push_back((double)surface.piyn[i]);

			Pi.push_back((double)surface.Pi[i]);
		}
	}
}


void HYDRO::start_hydro(int argc, char **argv)
{
	//HydroInitialTmunu init_tmunu;						// Derek: renamed this class initial_energy_vector
  														// does not appear to be initialized in cpu-vh

	initial_energy_vector initial_energy_profile;



	// the rest of my code:

	bool sample_parameters = false;						// default: use model parameters in parameters/
	int sample = 0;
	random_model_parameters random;

#ifdef RANDOM_MODEL_PARAMETERS
	if(argc == 2)										// load sample model parameters
	{
		sample_parameters = true;

		istringstream iss(argv[1]);

        if(!(iss >> sample))
        {
  			printf("main error: parameter sample index %s invalid\n", argv[1]);
  			exit(-1);
  		}
		else if(!iss.eof())
		{
			printf("main error: trailing characeters after sample index %s.\n", argv[1]);
			exit(-1);
		}
		if(sample <= 0)
		{
			printf("main error: parameter sample index %d must be greater than zero\n", sample);
			exit(-1);
		}

		random = load_random_model_parameters(sample);
	}
#endif

	hydro_parameters hydro = load_hydro_parameters(sample_parameters, random);
	initial_condition_parameters initial = load_initial_condition_parameters(sample_parameters, random);
	lattice_parameters lattice = load_lattice_parameters(hydro, initial, sample_parameters, sample);


	if(hydro.run_hydro)									// main hydro simulation (freezeout surface empty by default)
	{
		print_hydro_mode(hydro);
		store_freezeout_surface(run_hydro(lattice, initial, hydro, sample));
	}
	else
	{
		output_semi_analytic_solution_if_any(lattice, initial, hydro);
	}

	printf("\nEnd of hydro simulation\n");
}




