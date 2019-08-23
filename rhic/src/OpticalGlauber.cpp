/*
 * GlauberModel.c
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#include "../include/OpticalGlauber.h"
#include "../include/Parameters.h"

#include <gsl/gsl_integration.h>

const double n0 = 0.17; 				// nuclear density in fm^(-3)
const double d = 0.54; 					// nuclear thickness in fm
const double fermiToMilliBarns = 0.1;	// conversion

struct nuclearThicknessFunctionParams
{
	double x;		// (x,y) coordinates and nuclei
	double y;
	double A;
};

double woodsSaxonDistribution(double r, double A)
{
	double Rn = 1.12 * pow(A, 1./3.)  -  0.86 * pow(A, -1./3.);		// nuclear radius in fm

	return n0 / (1.0 + exp((r - Rn) / d));
}

double nuclearThicknessFunctionIntegrand(double s, void * params)
{
	// declare a struct (don't understand what the * params input means)
	struct nuclearThicknessFunctionParams * ta_params = (struct nuclearThicknessFunctionParams *) params;

	double x = ta_params->x;
	double y = ta_params->y;
	double A = ta_params->A;

	double z = (1.0 - s) / s;		// variable substitution dz = ds / s2
	double s2 = s * s;

	return woodsSaxonDistribution(sqrt(x * x  +  y * y  +  z * z), A) / s2;
}

double nuclearThicknessFunction(double x, double y, double A)
{
	double result, error;
	double epsabs = 0.0;
	double epsrel = 1.0e-8;
	int n = 1000;				// number of integration points

	struct nuclearThicknessFunctionParams params = {x, y, A};

	gsl_function F;
	F.function = &nuclearThicknessFunctionIntegrand;	// this is the function set to the integrand
	F.params = &params;

	// gsl integrate wrt s
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(n);
	gsl_integration_qags(&F, 0.0, 1.0, epsabs, epsrel, n, w, &result, &error);
	gsl_integration_workspace_free(w);

	return 2.0 * result;
}

double binaryCollisionPairs(double TA, double TB, double snn)
{
	return snn * TA * TB;
}

double woundedNucleons(double TA, double TB, double A, double B, double snn)
{
	double part_A = 1.0  -  snn * TA * fermiToMilliBarns / A;
	double part_B = 1.0  -  snn * TB * fermiToMilliBarns / B;

	return TA * (1.0 - pow(1.0 - part_B, B))  +  TB * (1.0 - pow(1.0 - part_A, A));
}


// Mixture of wounded nucleon and binary collision energy density profile
void Optical_Glauber_energy_density_transverse_profile(double * const __restrict__ energyDensityTransverse, int nx, int ny, double dx, double dy, initial_condition_parameters initial)
{
	// nuclear collision parameters
	double A = initial.numberOfNucleonsPerNuclei;		// assume A = B
	double B = A;
	double b = initial.impactParameter;
	double snn = initial.scatteringCrossSectionNN;
	double alpha = initial.fractionOfBinaryCollisions;


	// normalization factors
	double TA_norm = nuclearThicknessFunction(0.0, 0.0, A);
	double TB_norm = nuclearThicknessFunction(0.0, 0.0, B);
	double binary_norm = 1.0 / binaryCollisionPairs(TA_norm, TB_norm, snn);
	double wounded_norm = 1.0 / woundedNucleons(TA_norm, TB_norm, A, B, snn);

	double b_over2 = b / 2.0;

	// normalized transverse energy density profile
	for(int i = 0; i < nx; i++)
	{
		double x = (i  -  (nx - 1.0)/2.0) * dx;

		for(int j = 0; j < ny; j++)
		{
			double y = (j  -  (ny - 1.0)/2.0) * dy;

			double TA = nuclearThicknessFunction(x - b_over2, y, A);
			double TB = nuclearThicknessFunction(x + b_over2, y, B);

			// binary collision and wounded nucleon energy density profiles
			double energy_density_binary = alpha * binary_norm * binaryCollisionPairs(TA, TB, snn);
			double energy_density_wounded = (1.0 - alpha) * wounded_norm * woundedNucleons(TA, TB, A, B, snn);

			energyDensityTransverse[i  +  j * nx] = (energy_density_binary + energy_density_wounded);
		}
	}
}



