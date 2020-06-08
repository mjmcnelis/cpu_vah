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


void HYDRO::set_initial_energy_density_vector(std::vector<double> energy_vector)
{
	initial_energy_density_vector = energy_vector;
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

	if(hydro.run_hydro)
	{
		print_hydro_mode(hydro);
		run_hydro(lattice, initial, hydro, sample);		// main hydro simulation
	}
	else
	{
		output_semi_analytic_solution_if_any(lattice, initial, hydro);
	}
}




