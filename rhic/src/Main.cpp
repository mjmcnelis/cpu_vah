/*
-------------------------------------------------------
| Code        	| CPU VAH
-------------------------------------------------------
| Authors      	| Dennis Bazow, Mike McNelis
-------------------------------------------------------
| Date created	| 10/12/2015 by Dennis Bazow
-------------------------------------------------------
| Last edited	| 5/19/2020 by Mike McNelis
-------------------------------------------------------
| Description 	| A viscous anisotropic hydrodynamic
|				| simulation of a heavy-ion collision
-------------------------------------------------------
 */

#include <iostream>
#include <sstream>
#include <string>
#include "../include/Parameters.h"
#include "../include/Print.h"
#include "../include/Macros.h"
#include "../include/Hydrodynamics.h"
#include "../include/FileIO.h"
using namespace std;

int main(int argc, char *argv[])
{
	bool sample_parameters = false;						// default
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
		run_hydro(lattice, initial, hydro, sample);		// main function
	}
	else
	{
		output_semi_analytic_solution_if_any(lattice, initial, hydro);
	}

	return 0;
}

