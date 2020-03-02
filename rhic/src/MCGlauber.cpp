#include <stdlib.h>
#include <iostream>
#include <random>
//#include <gsl/gsl_integration.h>
#include "../include/OpticalGlauber.h"
#include "../include/MCGlauber.h"
#include "../include/Parameters.h"
using namespace std;


//const double SIG0 = 0.46;			// gaussian width
const double SIG0 = 0.8;			// nucleon width
const double two_pi = 2.0 * M_PI;

// generates A samples corresponding to Woods-Saxon nucleus with A nucleons
// results are loaded into x and y arrays (z is not used in current context)
// note: assumes that the the random number generator has already been seeded
void sampleWoodsSaxon(int A, double * const __restrict__ x, double * const __restrict__ y)
{
	const double M = 1.01 * 4.5310551374155095;	// max value of weight for Pb-Pb (with a fudge factor)

	const double rmax = 20;		// max radius (technically infinity)

	int nucleon = 0;			// number of sampled nucleons

	// sample nucleon positions
	while(nucleon < A)
	{
		double r = rmax * ((double)rand()) / ((double)RAND_MAX);

		double weight = r * r * woodsSaxonDistribution(r, A) / M;
		double propose = ((double)rand()) / ((double)RAND_MAX);

		if(fabs(weight - 0.5) > 0.5)
		{
			printf("sampleWoodsSaxon error: weight = %lf\n", weight);
			exit(-1);
		}

		if(propose < weight)
		{
			double costheta = 2.0 * (((double)rand()) / ((double)RAND_MAX) - 0.5);
			double phi = two_pi * ((double)rand()) / ((double)RAND_MAX);

			double sintheta = sqrt(1.0 - costheta * costheta);

			x[nucleon] = r * sintheta * cos(phi);	// set sampled nucleon positions
			y[nucleon] = r * sintheta * sin(phi);

			nucleon++;
		}
	}
}

int numberWoundedNucleons(int A, double b, double * const __restrict__ x, double * const __restrict__ y, double snn) {
	double x1[A], y1[A], x2[A], y2[A];
	int l1[A], l2[A];
	sampleWoodsSaxon(A,x1,y1);
	sampleWoodsSaxon(A,x2,y2);
	for (int i=0; i<A; i++) {
		x1[i] -= b/2;
		x2[i] += b/2;
	}
	double dn = sqrt(0.1*snn/M_PI);
	for(int i = 0; i < A; ++i) {
		l1[i] = 0;
		l2[i] = 0;
	}
	for(int i = 0; i < A; ++i) {
		for (int j = 0; j < A; ++j) {
			double dist = pow(x1[i]-x2[j],2);
			dist += pow(y1[i]-y2[j],2);
			dist = sqrt(dist);
			if (dist < dn) {
				l1[i] += 1;
				l2[j] += 1;
			}
		}
	}
	int n = 0;
	for(int i = 0; i < A; ++i) {
		if (l1[i] > 0) {
			x[n] = x1[i];
			y[n] = y1[i];
			n++;
		}
		if (l2[i] > 0) {
			x[n] = x2[i];
			y[n] = y2[i];
			n++;
		}
	}
	return n;
}

void MC_Glauber_energy_density_transverse_profile(double * const __restrict__ energyDensityTransverse, int nx, int ny, double dx, double dy, initial_condition_parameters initial)
{
	int A = initial.numberOfNucleonsPerNuclei;
	double b = initial.impactParameter;
	double snn = initial.scatteringCrossSectionNN;

	double xp[2 * A], yp[2 * A];
   	srand(1328398221);

	int wounded_nucleons = numberWoundedNucleons(A, b, xp, yp, snn);
	printf(": %d wounded nucleons ", wounded_nucleons);

	for(int i = 0; i < nx; i++)
	{
		double x = (i - (nx - 1.0)/2.0) * dx;

   		for(int j = 0; j < ny; j++)
   	 	{
   	 		double y = (j - (ny - 1.0)/2.0) * dy;

   	 		double eT = 0.0;
	      	for(int n = 0; n < wounded_nucleons; ++n)
	      	{
	         	double dx = x - xp[n];
	            double dy = y - yp[n];

	            eT += exp( - (dx * dx  +  dy * dy) / (2.0 * SIG0 * SIG0));  // gaussion bump in energy density
	        }

	        energyDensityTransverse[i + j * nx] = eT;
		}
	}
}


