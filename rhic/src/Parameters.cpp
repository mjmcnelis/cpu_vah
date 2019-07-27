#include <stdio.h>
#include <stdlib.h>
#include "../include/Parameters.h"

void getIntegerProperty(config_t * cfg, const char * propName, int * propValue)
{
	config_lookup_int(cfg, propName, propValue);
}


void getDoubleProperty(config_t * cfg, const char * propName, double * propValue)
{
	config_lookup_float(cfg, propName, propValue);
}


void loadLatticeParameters(config_t *cfg, void * params)
{
	char fname[255];
	sprintf(fname, "parameters/lattice.properties");

	if(!config_read_file(cfg, fname))
	{
		fprintf(stderr, "No configuration file %s found for lattice parameters - %s.\n", fname, config_error_text(cfg));
		exit(-1);
	}

	int numLatticePointsX, numLatticePointsY, numLatticePointsRapidity, numProperTimePoints;
	double latticeSpacingX, latticeSpacingY, latticeSpacingRapidity, latticeSpacingProperTime;

	// get the lattice parameters
	getIntegerProperty(cfg, "numLatticePointsX", 		&numLatticePointsX);
	getIntegerProperty(cfg, "numLatticePointsY", 		&numLatticePointsY);
	getIntegerProperty(cfg, "numLatticePointsRapidity", &numLatticePointsRapidity);
	getIntegerProperty(cfg, "numProperTimePoints", 		&numProperTimePoints);
	getDoubleProperty(cfg,  "latticeSpacingX", 			&latticeSpacingX);
	getDoubleProperty(cfg,  "latticeSpacingY", 			&latticeSpacingY);
	getDoubleProperty(cfg,  "latticeSpacingRapidity", 	&latticeSpacingRapidity);
	getDoubleProperty(cfg,  "latticeSpacingProperTime", &latticeSpacingProperTime);

	// printf("\nLattice parameters:");
	// printf("\n-------------------\n");
	// printf("numLatticePointsX        = %d\n", numLatticePointsX);
	// printf("numLatticePointsY        = %d\n", numLatticePointsY);
	// printf("numLatticePointsRapidity = %d\n", numLatticePointsRapidity);
	// printf("numProperTimePoints      = %d\n", numProperTimePoints);
	// printf("latticeSpacingX          = %.3f\n", latticeSpacingX);
	// printf("latticeSpacingY          = %.3f\n", latticeSpacingY);
	// printf("latticeSpacingRapidity   = %.3f\n", latticeSpacingRapidity);
	// printf("latticeSpacingProperTime = %.3f\n", latticeSpacingProperTime);
	// printf("\n");

	struct LatticeParameters * lattice = (struct LatticeParameters *) params;
	lattice->numLatticePointsX 			= numLatticePointsX;
	lattice->numLatticePointsY 			= numLatticePointsY;
	lattice->numLatticePointsRapidity	= numLatticePointsRapidity;
	lattice->numProperTimePoints 		= numProperTimePoints;
	lattice->latticeSpacingX 			= latticeSpacingX;
	lattice->latticeSpacingY 			= latticeSpacingY;
	lattice->latticeSpacingRapidity 	= latticeSpacingRapidity;
	lattice->latticeSpacingProperTime 	= latticeSpacingProperTime;
}


void loadInitialConditionParameters(config_t *cfg, void * params)
{
	char fname[255];
	sprintf(fname, "parameters/ic.properties");

	if(!config_read_file(cfg, fname))
	{
		fprintf(stderr, "No configuration file  %s found for initial condition parameters - %s.\n", fname, config_error_text(cfg));
		exit(-1);
	}

	int initialConditionType, numberOfNucleonsPerNuclei;
	double initialCentralTemperatureGeV, scatteringCrossSectionNN, impactParameter, fractionOfBinaryCollisions, rapidityVariance, rapidityMean;

	// get initial condition parameters
	getIntegerProperty(cfg, "initialConditionType", 		&initialConditionType);
	getIntegerProperty(cfg, "numberOfNucleonsPerNuclei", 	&numberOfNucleonsPerNuclei);
	getDoubleProperty(cfg, 	"initialCentralTemperatureGeV",	&initialCentralTemperatureGeV);
	getDoubleProperty(cfg, 	"scatteringCrossSectionNN", 	&scatteringCrossSectionNN);
	getDoubleProperty(cfg, 	"impactParameter", 				&impactParameter);
	getDoubleProperty(cfg,	"fractionOfBinaryCollisions", 	&fractionOfBinaryCollisions);
	getDoubleProperty(cfg, 	"rapidityVariance", 			&rapidityVariance);
	getDoubleProperty(cfg, 	"rapidityMean", 				&rapidityMean);

	// printf("Initial condition parameters:");
	// printf("\n-----------------------------\n");
	// printf("initialConditionType         = %d\n", initialConditionType);
	// printf("numberOfNucleonsPerNuclei    = %d\n", numberOfNucleonsPerNuclei);
	// printf("initialCentralTemperatureGeV = %.3f\n", initialCentralTemperatureGeV);
	// printf("scatteringCrossSectionNN     = %.2f\n", scatteringCrossSectionNN);
	// printf("impactParameter              = %.2f\n", impactParameter);
	// printf("fractionOfBinaryCollisions   = %.2f\n", fractionOfBinaryCollisions);
	// printf("rapidityVariance             = %.3f\n", rapidityVariance);
	// printf("rapidityMean                 = %.2f\n", rapidityMean);
	// printf("\n");

	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) params;
	initCond->initialConditionType 			= initialConditionType;
	initCond->numberOfNucleonsPerNuclei 	= numberOfNucleonsPerNuclei;
	initCond->initialCentralTemperatureGeV 	= initialCentralTemperatureGeV;
	initCond->scatteringCrossSectionNN 		= scatteringCrossSectionNN;
	initCond->impactParameter 				= impactParameter;
	initCond->fractionOfBinaryCollisions 	= fractionOfBinaryCollisions;
	initCond->rapidityVariance 				= rapidityVariance;
	initCond->rapidityMean 					= rapidityMean;
}


void loadHydroParameters(config_t *cfg, void * params)
{
	char fname[255];
	sprintf(fname, "parameters/hydro.properties");

	if(!config_read_file(cfg, fname))
	{
		fprintf(stderr, "No configuration file  %s found for hydrodynamic parameters - %s.\n", fname, config_error_text(cfg));
		exit(-1);
	}

	double tau_initial, shear_viscosity, freezeoutTemperatureGeV;

	// get hydro parameters
	getDoubleProperty(cfg, "tau_initial", 				&tau_initial);
	getDoubleProperty(cfg, "shear_viscosity", 			&shear_viscosity);
	getDoubleProperty(cfg, "freezeoutTemperatureGeV",	&freezeoutTemperatureGeV);

	// printf("Hydro parameters:");
	// printf("\n-----------------\n");
	// printf("tau_initial             = %.3f\n", tau_initial);
	// printf("shear_viscosity         = %.3f\n", shear_viscosity);
	// printf("freezeoutTemperatureGeV = %.3f\n", freezeoutTemperatureGeV);
	// printf("\n");

	// set hydro struct
	struct HydroParameters * hydro = (struct HydroParameters *) params;
	hydro->tau_initial             = tau_initial;
	hydro->shear_viscosity         = shear_viscosity;
	hydro->freezeoutTemperatureGeV = freezeoutTemperatureGeV;
}



