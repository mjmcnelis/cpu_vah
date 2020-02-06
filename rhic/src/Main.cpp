/*
-------------------------------------------------------
| Code        	| CPU VAH
-------------------------------------------------------
| Authors      	| Dennis Bazow, Mike McNelis
-------------------------------------------------------
| Date created	| 10/12/15
-------------------------------------------------------
| Last edited	| 2/6/20
-------------------------------------------------------
| Description 	| A viscous anisotropic hydrodynamic
|				| simulation of a heavy-ion collision
-------------------------------------------------------
 */

#include <iostream>
#include "../include/Parameters.h"
#include "../include/Print.h"
#include "../include/Hydrodynamics.h"
#include "../include/FileIO.h"

int main()
{
	hydro_parameters hydro = load_hydro_parameters();
	lattice_parameters lattice = load_lattice_parameters(hydro);
	initial_condition_parameters initial = load_initial_condition_parameters();

	if(hydro.run_hydro)
	{
		print_hydro_mode(hydro);
		run_hydro(lattice, initial, hydro);		// main function
	}
	else
	{
		output_semi_analytic_solution_if_any(lattice, initial, hydro);
	}

	return 0;
}

