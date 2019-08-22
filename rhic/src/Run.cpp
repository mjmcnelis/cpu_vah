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
#include <libconfig.h>
#include "../include/DynamicalVariables.h"
#include "../include/Parameters.h"
#include "../include/Hydrodynamics.h"


int main(int argc, char **argv)
{
	struct LatticeParameters latticeParams;
	struct InitialConditionParameters initCondParams;
	struct HydroParameters hydroParams;

	config_t latticeConfig, initCondConfig, hydroConfig, hydroConfig_2;

	config_init(&latticeConfig);
	config_init(&initCondConfig);
	config_init(&hydroConfig);
	config_init(&hydroConfig_2);

	// load parameters
	load_lattice_parameters(&latticeConfig, &latticeParams);
	loadInitialConditionParameters(&initCondConfig, &initCondParams);
	loadHydroParameters(&hydroConfig, &hydroParams);

	// test new loading parameters
	hydro_parameters hydro = load_hydro_parameters(&hydroConfig_2);

#ifdef ANISO_HYDRO
	printf("\n:::::::::::::::::::::::::::::::::::::::::::\n");
	printf(":::  Running viscous anisotropic hydro  :::\n");
	printf(":::::::::::::::::::::::::::::::::::::::::::\n\n");
#else
	printf("\n:::::::::::::::::::::::::::::::::::::::::::\n");
	printf(":::   Running 2nd order viscous hydro   :::\n");
	printf(":::::::::::::::::::::::::::::::::::::::::::\n\n");
#endif

	// main function
	run_hydro(&latticeParams, &initCondParams, &hydroParams, hydro);

	printf("\nFinished hydro\n");

	config_destroy(&latticeConfig);
	config_destroy(&initCondConfig);
	config_destroy(&hydroConfig);
	return 0;
}


