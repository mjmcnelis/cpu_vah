/*
 * LatticeParameters.h
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <libconfig.h>

#define N_GHOST_CELLS_M 2						// ghost cells on left side
#define N_GHOST_CELLS_P 2						// ghost cells on right side
#define N_GHOST_CELLS 4							// total 

void getIntegerProperty(config_t * cfg, const char * propName, int * propValue);
void getDoubleProperty(config_t * cfg, const char * propName, double * propValue);

struct LatticeParameters
{
	int numLatticePointsX;
	int numLatticePointsY;
	int numLatticePointsRapidity;

	int numComputationalLatticePointsX;			// includes ghost cells
	int numComputationalLatticePointsY;
	int numComputationalLatticePointsRapidity;

	int numProperTimePoints;

	double latticeSpacingX;
	double latticeSpacingY;
	double latticeSpacingRapidity;

	double latticeSpacingProperTime;
};

struct InitialConditionParameters
{
	int initialConditionType;

	int numberOfNucleonsPerNuclei;			// this is an A = B collision 

	double initialCentralTemperatureGeV;	// initial central temperature for a central A = B collision (b = 0)		
	double scatteringCrossSectionNN;
	double impactParameter;
	double fractionOfBinaryCollisions;

	// longitudinal energy density profile parameters
	double rapidityVariance; 	// \sigma^{2}_{\eta}
	double rapidityMean; 		// flat region around \ets_s = 0
};

struct HydroParameters
{
	double initialProperTimePoint;
	double shearViscosityToEntropyDensity;
	double freezeoutTemperatureGeV;
};

void loadLatticeParameters(config_t *cfg, void * params);
void loadInitialConditionParameters(config_t *cfg, void * params);
void loadHydroParameters(config_t *cfg, void * params);

#endif 
