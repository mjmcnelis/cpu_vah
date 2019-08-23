/*
-------------------------------------------------------
| Code        	| CPU VAH
-------------------------------------------------------
| Authors      	| Dennis Bazow, Mike McNelis
-------------------------------------------------------
| Date created	| 10/12/15
-------------------------------------------------------
| Last edited	| 8/2/19
-------------------------------------------------------
| Description 	| A viscous anisotropic hydrodynamic
|				| simulation of a heavy-ion collision
-------------------------------------------------------
 */
#include <stdlib.h>
#include <stdio.h>
#include "../include/Parameters.h"
#include "../include/Print.h"
#include "../include/Hydrodynamics.h"

int main()
{
	lattice_parameters lattice = load_lattice_parameters();
	initial_condition_parameters initial = load_initial_condition_parameters();
	hydro_parameters hydro = load_hydro_parameters();

	print_hydro_mode(hydro);

	run_hydro(lattice, initial, hydro);		// main function

	printf("\nFinished hydro\n");
	return 0;
}

