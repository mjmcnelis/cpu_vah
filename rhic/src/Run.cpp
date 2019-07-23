

/*
-------------------------------------------------------
| Code        	| CPU VAH
-------------------------------------------------------
| Authors      	| Dennis Bazow, Mike McNelis
-------------------------------------------------------
| Date created	| 10/12/15
-------------------------------------------------------
| Last edited	| 7/23/19
-------------------------------------------------------
| Version     	| 2.0
-------------------------------------------------------
| Description 	| A viscous anisotropic hydrodynamic
|				| simulation of a heavy-ion collision
-------------------------------------------------------
 */


#include <stdlib.h>
#include <stdio.h> 			// for printf
#include <sys/time.h> 		// for timing
// #include <unistd.h>		// for current working directory
#include <libconfig.h>
#include <iostream>

using namespace std;

#include "../include/Parameters.h"
#include "../include/Hydrodynamics.h"


int main(int argc, char **argv)
{
	struct LatticeParameters latticeParams;
	struct InitialConditionParameters initCondParams;
	struct HydroParameters hydroParams;

	config_t latticeConfig;
	config_t initCondConfig;
	config_t hydroConfig;

	config_init(&latticeConfig);
	config_init(&initCondConfig);
	config_init(&hydroConfig);

	// load parameters from rhic-conf files
	loadLatticeParameters(&latticeConfig, &latticeParams);
	loadInitialConditionParameters(&initCondConfig, &initCondParams);
	loadHydroParameters(&hydroConfig, &hydroParams);


	printf("\n:::::::::::::::::::::::::::::::::::::::::::\n");
	printf(":::  Running viscous anisotropic hydro  :::\n");
	printf(":::::::::::::::::::::::::::::::::::::::::::\n\n");

	// main function
	run_hydro(&latticeParams, &initCondParams, &hydroParams);


	printf("\nFinished hydro\n");

	config_destroy(&latticeConfig);
	config_destroy(&initCondConfig);
	config_destroy(&hydroConfig);

	return 0;
}
