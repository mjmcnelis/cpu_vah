/*
 * InitialConditionParameters.c
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

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
	// read the file
	char fname[255];
	sprintf(fname, "parameters/lattice.properties");

	if(!config_read_file(cfg, fname))
	{
		fprintf(stderr, "No configuration file %s found for lattice parameters - %s.\n", fname, config_error_text(cfg));
		exit(-1);
	}

	int numLatticePointsX;
	int numLatticePointsY;
	int numLatticePointsRapidity;
	int numProperTimePoints;

	double latticeSpacingX;
	double latticeSpacingY;
	double latticeSpacingRapidity;
	double latticeSpacingProperTime;


	// get the lattice parameters
	getIntegerProperty(cfg, "numLatticePointsX", &numLatticePointsX);
	getIntegerProperty(cfg, "numLatticePointsY", &numLatticePointsY);
	getIntegerProperty(cfg, "numLatticePointsRapidity", &numLatticePointsRapidity);
	getIntegerProperty(cfg, "numProperTimePoints", &numProperTimePoints);

	getDoubleProperty(cfg, "latticeSpacingX", &latticeSpacingX);
	getDoubleProperty(cfg, "latticeSpacingY", &latticeSpacingY);
	getDoubleProperty(cfg, "latticeSpacingRapidity", &latticeSpacingRapidity);
	getDoubleProperty(cfg, "latticeSpacingProperTime", &latticeSpacingProperTime);


	printf("\nLattice parameters:");
	printf("\n-------------------\n");
	printf("numLatticePointsX        = %d\n", numLatticePointsX);
	printf("numLatticePointsY        = %d\n", numLatticePointsY);
	printf("numLatticePointsRapidity = %d\n", numLatticePointsRapidity);
	printf("numProperTimePoints      = %d\n", numProperTimePoints);
	printf("latticeSpacingX          = %.3f\n", latticeSpacingX);
	printf("latticeSpacingY          = %.3f\n", latticeSpacingY);
	printf("latticeSpacingRapidity   = %.3f\n", latticeSpacingRapidity);
	printf("latticeSpacingProperTime = %.3f\n", latticeSpacingProperTime);
	printf("\n");


	// set the lattice struct
	struct LatticeParameters * lattice = (struct LatticeParameters *) params;

	lattice->numLatticePointsX = numLatticePointsX;
	lattice->numLatticePointsY = numLatticePointsY;
	lattice->numLatticePointsRapidity = numLatticePointsRapidity;
	lattice->numProperTimePoints = numProperTimePoints;

	lattice->numComputationalLatticePointsX = numLatticePointsX + N_GHOST_CELLS;
	lattice->numComputationalLatticePointsY = numLatticePointsY + N_GHOST_CELLS;
	lattice->numComputationalLatticePointsRapidity = numLatticePointsRapidity + N_GHOST_CELLS;

	lattice->latticeSpacingX = latticeSpacingX;
	lattice->latticeSpacingY = latticeSpacingY;
	lattice->latticeSpacingRapidity = latticeSpacingRapidity;
	lattice->latticeSpacingProperTime = latticeSpacingProperTime;
}


void loadInitialConditionParameters(config_t *cfg, void * params)
{
	// read the file
	char fname[255];
	sprintf(fname, "parameters/ic.properties");

	if(!config_read_file(cfg, fname))
	{
		fprintf(stderr, "No configuration file  %s found for initial condition parameters - %s.\n", fname, config_error_text(cfg));
		exit(-1);
	}

	int initialConditionType;
	int numberOfNucleonsPerNuclei;

	double initialCentralTemperatureGeV;
	double scatteringCrossSectionNN;
	double impactParameter;
	double fractionOfBinaryCollisions;
	double rapidityVariance;
	double rapidityMean;


	// get initial condition parameters
	getIntegerProperty(cfg, "initialConditionType", &initialConditionType);
	getIntegerProperty(cfg, "numberOfNucleonsPerNuclei", &numberOfNucleonsPerNuclei);

	getDoubleProperty(cfg, "initialCentralTemperatureGeV", &initialCentralTemperatureGeV);
	getDoubleProperty(cfg, "scatteringCrossSectionNN", &scatteringCrossSectionNN);
	getDoubleProperty(cfg, "impactParameter", &impactParameter);
	getDoubleProperty(cfg, "fractionOfBinaryCollisions", &fractionOfBinaryCollisions);
	getDoubleProperty(cfg, "rapidityVariance", &rapidityVariance);
	getDoubleProperty(cfg, "rapidityMean", &rapidityMean);


	printf("Initial condition parameters:");
	printf("\n-----------------------------\n");
	printf("initialConditionType         = %d\n", initialConditionType);
	printf("numberOfNucleonsPerNuclei    = %d\n", numberOfNucleonsPerNuclei);
	printf("initialCentralTemperatureGeV = %.3f\n", initialCentralTemperatureGeV);
	printf("scatteringCrossSectionNN     = %.2f\n", scatteringCrossSectionNN);
	printf("impactParameter              = %.2f\n", impactParameter);
	printf("fractionOfBinaryCollisions   = %.2f\n", fractionOfBinaryCollisions);
	printf("rapidityVariance             = %.3f\n", rapidityVariance);
	printf("rapidityMean                 = %.2f\n", rapidityMean);
	printf("\n");


	// this is a pointer so the address values are set (still not sure how params works...)

	// set the initial condition struct
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) params;

	initCond->initialConditionType = initialConditionType;
	initCond->numberOfNucleonsPerNuclei = numberOfNucleonsPerNuclei;

	initCond->initialCentralTemperatureGeV = initialCentralTemperatureGeV;
	initCond->scatteringCrossSectionNN = scatteringCrossSectionNN;
	initCond->impactParameter = impactParameter;
	initCond->fractionOfBinaryCollisions = fractionOfBinaryCollisions;
	initCond->rapidityVariance = rapidityVariance;
	initCond->rapidityMean = rapidityMean;
}


void loadHydroParameters(config_t *cfg, void * params)
{
	// read the file
	char fname[255];
	sprintf(fname, "parameters/hydro.properties");

	if(!config_read_file(cfg, fname))
	{
		fprintf(stderr, "No configuration file  %s found for hydrodynamic parameters - %s.\n", fname, config_error_text(cfg));
		exit(-1);
	}

	double initialProperTimePoint;
	double shearViscosityToEntropyDensity;
	double freezeoutTemperatureGeV;

	// get hydro parameters
	getDoubleProperty(cfg, "initialProperTimePoint", &initialProperTimePoint);
	getDoubleProperty(cfg, "shearViscosityToEntropyDensity", &shearViscosityToEntropyDensity);
	getDoubleProperty(cfg, "freezeoutTemperatureGeV", &freezeoutTemperatureGeV);

	printf("Hydro parameters:");
	printf("\n-----------------\n");
	printf("initialProperTimePoint         = %.3f\n", initialProperTimePoint);
	printf("shearViscosityToEntropyDensity = %.3f\n", shearViscosityToEntropyDensity);
	printf("freezeoutTemperatureGeV        = %.3f\n", freezeoutTemperatureGeV);
	printf("\n");


	// set hydro struct
	struct HydroParameters * hydro = (struct HydroParameters *) params;

	hydro->initialProperTimePoint = initialProperTimePoint;
	hydro->shearViscosityToEntropyDensity = shearViscosityToEntropyDensity;
	hydro->freezeoutTemperatureGeV = freezeoutTemperatureGeV;
}



