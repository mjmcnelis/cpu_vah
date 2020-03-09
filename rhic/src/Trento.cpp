#include <stdlib.h>
#include <iostream>
#include <random>
#include <chrono>
#include <vector>
#include "../include/MCGlauber.h"
#include "../include/Parameters.h"
using namespace std;


//const double SIG0 = 0.46;			// gaussian width
const double w = 0.5;				// nucleon width
const double sgg = 20.4973;			// gluon-gluon cross section? (difference from snn?)


inline double canonical(default_random_engine & generator)
{
	// random number between [0,1)
	return generate_canonical<double, numeric_limits<double>::digits>(generator);
}


double woods_saxon(double r, double A)
{
	double n0 = 0.17; 												// nuclear density [fm^-3]
	double Rn = 1.12 * pow(A, 1./3.)  -  0.86 * pow(A, -1./3.);		// nuclear radius in fm
	double d = 0.54; 												// surface thickness in fm (Pb only?)

	return n0 / (1. + exp((r - Rn) / d));
}


double Tp(double distance2, double w2)
{
	return exp(- distance2 / (2. * w2)) / (2. * M_PI * w2);			// nucleon thickness function (right normalization?)
}


double Tpp_overlap(double distance2, double w2)
{
	return exp(- distance2 / (4. * w2)) / (4. * M_PI * w2);			// nucleon-nucleon overlap thickness function
}


void sample_nucleon_transverse_positions_in_nucleus(int A, double * x, double * y, default_random_engine & generator)
{
	// sample the transverse positions of nucleons in nucleus A with probability r^2.woods_saxon(r, A).dr.dphi.dcostheta

  	uniform_real_distribution<double> phi_distribution(0., 2. * M_PI);
	uniform_real_distribution<double> costheta_distribution(-1., nextafter(1., numeric_limits<double>::max()));

	// r2_woods_saxon_max is A dependent; pass separately (special cases 208, 197, etc)

	double r2_woods_saxon_max = 4.5310551374155095;	    			// max value of r^2.woods_saxon(r, A) (for A = 208)

	int n = 0;

	while(n < A)													// keep sampling until all nucleons positions sampled
	{		
		double r = 20. * canonical(generator);						// sample r uniformly between 0 and 20
		double phi = phi_distribution(generator);					// sample angles
	    double costheta = costheta_distribution(generator);
		double sintheta = sqrt(fabs(1. - costheta * costheta));

		double weight = r * r * woods_saxon(r, A) / r2_woods_saxon_max;

		if(fabs(weight - 0.5) > 0.5)
		{
			printf("sample_nucleon_transverse_positions_in_nucleus error: weight = %lf\n", weight);
			exit(-1);
		}

		if(canonical(generator) < weight)			// let's not enforce |x1 - x2| >= d_min yet
		{
	     	x[n] = r * sintheta * cos(phi);			// set nucleon's transverse position
			y[n] = r * sintheta * sin(phi);
			n++;
		}
	}
}



void wounded_nucleons(int A, int B, vector<double> * xA_wounded, vector<double> * yA_wounded, vector<double> * xB_wounded, vector<double> * yB_wounded, double b, double w, double sigma_gg, default_random_engine & generator)
{
	// computes the transverse positions of wounded nucleons in A and B

	
	double xA[A];		// sampled transverse positions of nucleons in A
	double yA[A];		
	sample_nucleon_transverse_positions_in_nucleus(A, xA, yA, generator);

	double xB[B];		// sampled transverse positions of nucleons in A
	double yB[B];
	sample_nucleon_transverse_positions_in_nucleus(B, xB, yB, generator);

	int wounded_nucleonsA[A];			// wounded nucleon tags
	int wounded_nucleonsB[B];			// 0 = not wounded, 1 = wounded

	for(int n = 0; n < A; n++)			// shift x positions by -/+ b/2
	{
		xA[n] -= b/2.;
		wounded_nucleonsA[n] = 0;
	}

	for(int n = 0; n < B; n++)
	{
		xB[n] += b/2.;
		wounded_nucleonsB[n] = 0;
	}

	
	for(int i = 0; i < A; i++)									// loop over all possible collision pairs and sample wounded nucleons						
	{
		for(int j = 0; j < B; j++)
		{
			double dx = xA[i] - xB[j];							// compute transverse distance between pairs
			double dy = yB[i] - yB[j];
			double distance_sq = dx * dx  +  dy * dy;

			double probability_collision = 1.  -  exp(- sigma_gg * Tpp_overlap(distance_sq, w * w));
			
			if(canonical(generator) < probability_collision)	// if pair collides, tag the wounded nucleons
			{
				wounded_nucleonsA[i] = 1;	
				wounded_nucleonsB[j] = 1;
			}

		} 
	} 

	
	for(int n = 0; n < A; n++)									// get transverse positions of wounded nucleons in A and B
	{
		if(wounded_nucleonsA[n])
		{
			xA_wounded->push_back(xA[n]);
			yA_wounded->push_back(yA[n]);
		}
	}

	for(int n = 0; n < B; n++)
	{
		if(wounded_nucleonsB[n])
		{
			xB_wounded->push_back(xB[n]);
			yB_wounded->push_back(yB[n]);
		}
	}
}


void trento_transverse_energy_profile(double * const __restrict__ energy_transverse, lattice_parameters lattice, initial_condition_parameters initial)
{
	// set the seed
	long unsigned seed = chrono::system_clock::now().time_since_epoch().count();
   	default_random_engine generator(seed);


	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	double dx = lattice.lattice_spacing_x;
	double dy = lattice.lattice_spacing_y;

	int A = initial.numberOfNucleonsPerNuclei;					// assumes A = B nuclear collisions
	int B = A;

	double b = initial.impactParameter;
	double snn = initial.scatteringCrossSectionNN;


	vector<double> * xA_wounded = new vector<double>();			// (x,y) positions of wounded nucleons in A and B
	vector<double> * yA_wounded = new vector<double>();
	vector<double> * xB_wounded = new vector<double>();
	vector<double> * yB_wounded = new vector<double>();  

	wounded_nucleons(A, B, xA_wounded, yA_wounded, xB_wounded, yB_wounded, b, w, sgg, generator); 	

	int nA_wounded = xA_wounded->size();
	int nB_wounded = xB_wounded->size();

	printf("%d wounded nucleons in A\n\n", nA_wounded);
	printf("%d wounded nucleons in B\n\n", nB_wounded);

	for(int i = 0; i < nx; i++)
	{
   		for(int j = 0; j < ny; j++)
   	 	{
   	 		double x = (i - (nx - 1.)/2.) * dx;			// fluid cell position
   	 		double y = (j - (ny - 1.)/2.) * dy;

   	 		double TA = 0;
   	 		double TB = 0;

			for(int n = 0; n < nA_wounded; n++)			// nuclear thickness function A 
			{
				double xA = (*xA_wounded)[n];
				double yA = (*yA_wounded)[n];

				double dx = x - xA;
				double dy = y - yA;
				double distance_sq = dx * dx  +  dy * dy;

				TA += Tp(distance_sq, w * w);
			}

			for(int n = 0; n < nB_wounded; n++)			// nuclear thickness function B = A
			{
				double xB = (*xB_wounded)[n];
				double yB = (*yB_wounded)[n];

				double dx = x - xB;
				double dy = y - yB;
				double distance_sq = dx * dx  +  dy * dy;
				
				TB += Tp(distance_sq, w * w);
			}

	        //double TR = sqrt(TA * TB);				// reduced nuclear thickness function (p = 0)
	        double TR = (TA + TB) / 2.;					// p = 1

	        double norm = 1;

	        energy_transverse[i + j * nx] = norm * TR;
		}
	}


	// I'll take care of the normalization = N / t0 separately
}






