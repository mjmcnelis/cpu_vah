#include <stdlib.h>
#include <iostream>
#include <random>
#include <chrono>
#include <vector>
#include "../include/Parameters.h"
#include "../include/DynamicalVariables.h"
#include "../include/Macros.h"
#include "../include/Hydrodynamics.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include "../include/OpenMP.h"

using namespace std;


inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


inline double Theta(double x)
{
	if(x > 0)
	{
		return 1.;
	}
	return 0;
}


inline double canonical(default_random_engine & generator)
{
	// random number between [0,1)
	return generate_canonical<double, numeric_limits<double>::digits>(generator);
}


double woods_saxon(double r, double A)
{
	double n0 = 0.17; 												// nuclear density [fm^-3]
	double R = 1.12 * pow(A, 1./3.)  -  0.86 * pow(A, -1./3.);		// nuclear radius in fm
	double d = 0.54; 												// surface thickness in fm (Pb only)

	return n0 / (1. + exp((r - R) / d));
}


inline double Tp(double distance2, double w2)
{
	return exp(- distance2 / (2. * w2)) / (2. * M_PI * w2);			// nucleon thickness function (right normalization?)
}


inline double Tpp_overlap(double distance2, double w2)
{
	return exp(- distance2 / (4. * w2)) / (4. * M_PI * w2);			// nucleon-nucleon overlap thickness function
}


double compute_sigma_gg(double w)
{
	// interpolation of sgg(w) data (note: this calculation is limited to sigma_NN = 6.4 fm^2 in Pb+Pb sqrt(s) = 2.76 TeV)

  	FILE * sgg_file;
  	sgg_file = fopen("tables/sgg.dat", "r");
  	if(sgg_file == NULL)
  	{
  		printf("Error: couldn't open sgg.dat file\n");
  	}

  	int points;
  	fscanf(sgg_file, "%d", &points);

  	double * width = (double *)calloc(points, sizeof(double));
  	double * sgg_values = (double *)calloc(points, sizeof(double));

	for(int i = 0; i < points; i++)
	{
		fscanf(sgg_file, "%lf\t%lf", &width[i], &sgg_values[i]);
	}

	fclose(sgg_file);

  	if(!(w >= width[0] && w <= width[points - 1]))
  	{
  		printf("compute_sigma_gg error: nucleon width w = %lf outside range [%lf,%lf]\n", w, width[0], width[points - 1]);
  		exit(-1);
  	}

	gsl_spline * sigma_gg_spline;
	gsl_interp_accel * accel = gsl_interp_accel_alloc();

    sigma_gg_spline = gsl_spline_alloc(gsl_interp_cspline, points);
    gsl_spline_init(sigma_gg_spline, width, sgg_values, points);

 	double sigma_gg = gsl_spline_eval(sigma_gg_spline, w, accel);

 	printf("sigma_gg = %lf\n\n", sigma_gg);

  	gsl_interp_accel_free(accel);
    gsl_spline_free(sigma_gg_spline);
    free(width);
	free(sgg_values);

	return sigma_gg;
}


void sample_nucleon_transverse_positions_in_nucleus(int A, double * xA, double * yA, double d_min, default_random_engine & generator)
{
	// sample the transverse positions of nucleons in nucleus A with probability r^2.woods_saxon(r, A).dr.dphi.dcostheta
	// also enforce minimum nucleon-nucleon separation

  	uniform_real_distribution<double> phi_distribution(0., 2. * M_PI);
	uniform_real_distribution<double> costheta_distribution(-1., nextafter(1., numeric_limits<double>::max()));

	double zA[A];												// longitudinal position (z not eta)

	double r_max = 20.;											// technically infinity but large enough coverage
	double r2_woods_saxon_max = 4.48620;						// max value of r^2.woods_saxon(r, 208) for Pb only

	int n = 0;													// number of sampled nucleons

	while(n < A)												// accept-reject sampling until all nucleons positions filled
	{
		double r = r_max * canonical(generator);				// sample r uniformly between 0 and r_max
		double phi = phi_distribution(generator);				// sample angles
	    double costheta = costheta_distribution(generator);
		double sintheta = sqrt(fabs(1. - costheta * costheta));

		double x = r * sintheta * cos(phi);						// cartesian position of proposed nucleon
		double y = r * sintheta * sin(phi);
		double z = r * costheta;

		double separation_weight = 1.;

		if(n > 0 && d_min > 0.)									// find distance btw proposed nucleon and nearest (already) sampled nucleon
		{
			double d = 1./0.;

			for(int i = 0; i < n; i++)
			{
				double xi = xA[i];
				double yi = yA[i];
				double zi = zA[i];

				double d_pair = sqrt((x - xi) * (x - xi)  + (y - yi) * (y - yi)  +  (z - zi) * (z - zi));

				d = fmin(d, d_pair);
			}

			separation_weight = Theta(d - d_min);
		}

		double weight = r * r * woods_saxon(r, A) / r2_woods_saxon_max * separation_weight;

		if(fabs(weight - 0.5) > 0.5)
		{
			printf("sample_nucleon_transverse_positions_in_nucleus warning: weight = %lf out of bounds\n", weight);
		}

		if(canonical(generator) < weight)
		{
	     	xA[n] = x;											// set nucleon's transverse position
			yA[n] = y;
			zA[n] = z;											// this is only needed for computing nucleon-nucleon separation
			n++;
		}
	}
}



void wounded_nucleons(int A, int B, vector<double> * xA_wounded, vector<double> * yA_wounded, vector<double> * xB_wounded, vector<double> * yB_wounded, double b, double w, double d_min, double sigma_gg, default_random_engine & generator)
{
	// computes the transverse positions of wounded nucleons in A and B

	double xA[A];												// sample nucleon transverse positions in nucleus A
	double yA[A];
	sample_nucleon_transverse_positions_in_nucleus(A, xA, yA, d_min, generator);

	double xB[B];												// sample nucleon transverse positions in nucleus B
	double yB[B];
	sample_nucleon_transverse_positions_in_nucleus(B, xB, yB, d_min, generator);

	int wounded_nucleonsA[A];									// wounded nucleon tags (default values are zero)
	int wounded_nucleonsB[B];									// 0 = not wounded, 1 = wounded

	for(int n = 0; n < A; n++)									// also shift x coordinates +/- b/2
	{
		xA[n] -= b/2.;
		wounded_nucleonsA[n] = 0;
	}

	for(int n = 0; n < B; n++)
	{
		xB[n] += b/2.;
		wounded_nucleonsB[n] = 0;
	}

	for(int i = 0; i < A; i++)									// loop over all possible collision pairs and tag wounded nucleons
	{
		for(int j = 0; j < B; j++)
		{
			double dx = xA[i] - xB[j];							// compute transverse distance between colliding pairs
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

	for(int n = 0; n < A; n++)									// collect transverse positions of wounded nucleons in A and B
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


void trento_transverse_energy_density_profile(double * const __restrict__ energy_density_transverse, lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro)
{
	// set seed and initialize sampler

	long unsigned seed = chrono::system_clock::now().time_since_epoch().count();

	if(initial.trento_fixed_seed)
	{
		seed = abs(initial.trento_fixed_seed);
	}

	printf("Trento seed = %lu\n", seed);

   	default_random_engine generator(seed);

	int nx = lattice.lattice_points_x;							// transverse grid
	int ny = lattice.lattice_points_y;
	double dx = lattice.lattice_spacing_x;
	double dy = lattice.lattice_spacing_y;

	int A = initial.nucleus_A;									// number of nucleons in nuclei A, B
	int B = initial.nucleus_B;

	double b = initial.impact_parameter;						// impact parameter [fm]
	double N = initial.trento_normalization_GeV;				// normalization factor [GeV]
	double w = initial.trento_nucleon_width;					// nucleon width [fm]
	double sigma_gg = compute_sigma_gg(w);						// compute parton cross section in trento model [fm^2]

	double d_min = initial.trento_min_nucleon_distance;			// minimum nucleon-nucleon separation in nucleus [fm]
	double p = initial.trento_geometric_parameter;				// geometric parameter p
	double sigma_k = initial.trento_gamma_standard_deviation;	// gamma distribution standard deviation

	double t0 = hydro.tau_initial;								// initial time

	int event_averaging = initial.trento_average_over_events;	// option to do event averaging for smoother energy density profile
	int number_of_events = 1;

	if(event_averaging)
	{
		number_of_events = initial.trento_number_of_average_events;
	}

	double alpha = 1. / (sigma_k * sigma_k);					// alpha = k (shape parameter)
	double beta = sigma_k * sigma_k;							// beta = 1 / alpha
	gamma_distribution<double> gamma(alpha, beta);				// gamma distribution (mean = 1, std = sigma_k)

	vector<double> * xA_wounded = new vector<double>();			// (x,y) positions of wounded nucleons in A and B
	vector<double> * yA_wounded = new vector<double>();
	vector<double> * xB_wounded = new vector<double>();
	vector<double> * yB_wounded = new vector<double>();

	for(int event = 0; event < number_of_events; event++)
	{
		wounded_nucleons(A, B, xA_wounded, yA_wounded, xB_wounded, yB_wounded, b, w, d_min, sigma_gg, generator);

		int nA_wounded = xA_wounded->size();
		int nB_wounded = xB_wounded->size();

		double gamma_A[nA_wounded];									// gamma distribution samples
		double gamma_B[nB_wounded];

		for(int n = 0; n < nA_wounded; n++)
		{
			if(sigma_k > 0.)
			{
				gamma_A[n] = gamma(generator);
			}
			else
			{
				gamma_A[n] = 1.;
			}
		}
		for(int n = 0; n < nB_wounded; n++)
		{
			if(sigma_k > 0.)
			{
				gamma_B[n] = gamma(generator);
			}
			else
			{
				gamma_B[n] = 1.;
			}
		}

		#pragma omp parallel for collapse(2)
		for(int j = 0; j < ny; j++)									// compute the transverse energy density profile
		{
	   		for(int i = 0; i < nx; i++)
	   	 	{
	   	 		double x = (i - (nx - 1.)/2.) * dx;					// fluid cell position x
	   	 		double y = (j - (ny - 1.)/2.) * dy;					// fluid cell position y

	   	 		double TA = 0;

				for(int n = 0; n < nA_wounded; n++)					// nuclear thickness function A [fm^-2]
				{
					double xA = (*xA_wounded)[n];
					double yA = (*yA_wounded)[n];

					double dx = x - xA;
					double dy = y - yA;
					double distance_sq = dx * dx  +  dy * dy;

					TA += gamma_A[n] * Tp(distance_sq, w * w);
				}

				double TB = 0;

				for(int n = 0; n < nB_wounded; n++)					// nuclear thickness function B = A [fm^-2]
				{
					double xB = (*xB_wounded)[n];
					double yB = (*yB_wounded)[n];

					double dx = x - xB;
					double dy = y - yB;
					double distance_sq = dx * dx  +  dy * dy;

					TB += gamma_B[n] * Tp(distance_sq, w * w);
				}

		        double TR = pow((pow(TA, p) + pow(TB, p)) / 2., 1./p);	// reduced nuclear thickness function

		        if(fabs(p) < 1.e-4)
		        {
		        	TR = sqrt(TA * TB);
		        }

		        energy_density_transverse[i + j * nx] += N * TR / (t0 * hbarc * (double)number_of_events);

			} // i
		} // j

		xA_wounded->clear();
		yA_wounded->clear();
		xB_wounded->clear();
		yB_wounded->clear();

		cout << "\r" << event + 1 << " / " << number_of_events << " events completed" << flush;

	} // events

	printf("\n\n");
}


void longitudinal_energy_density_extension(double * const __restrict__ eL, lattice_parameters lattice, initial_condition_parameters initial)
{
	int nz    = lattice.lattice_points_eta;
	double dz = lattice.lattice_spacing_eta;

	double eta_flat     = initial.rapidity_flat;
	double eta_variance = initial.rapidity_variance;

	// energy density profile along eta direction is a smooth plateu that exponentially decays when |eta| > eta_flat/2
	for(int k = 0; k < nz; k++)
	{
		double eta = (k - (nz - 1.)/2.) * dz;                   // fluid cell position eta
		double eta_scale = fabs(eta)  -  eta_flat / 2.;

		eL[k] = exp(- 0.5 * eta_scale * eta_scale / eta_variance * Theta(eta_scale));
	}
}


void set_trento_energy_density_and_flow_profile(lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro)
{
	precision e_min = hydro.energy_min;
	precision t0 = hydro.tau_initial;

	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	double eT[nx * ny];											// transverse energy density profile [fm^-4]
	double eL[nz];												// longitudinal extension [unitless]

	for(int i = 0; i < nx; i++)									// zero energy density profile
	{
		for(int j = 0; j < ny; j++)
		{
			eT[i + j * nx] = 0;
		}
	}

	for(int k = 0; k < nz; k++)
	{
		eL[k] = 0;
	}

#ifdef OPENMP
	double t1 = omp_get_wtime();
#endif

	trento_transverse_energy_density_profile(eT, lattice, initial, hydro);

	longitudinal_energy_density_extension(eL, lattice, initial);

	FILE *energy;														// make energy density block file (only if openmp turned off)

	if(initial.trento_average_over_events)
	{
		energy = fopen("tables/energy_density/new/e_block.dat", "w");	// block file of energy denisty [GeV/fm^-3]
		fprintf(energy, "%d\n%d\n%d\n", nx, ny, nz);					// write grid points at header ()
	}

	for(int k = 2; k < nz + 2; k++)
	{
        for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				precision e_s = eT[i - 2  +  (j - 2) * nx] * eL[k - 2];

				if(initial.trento_average_over_events)
				{
					fprintf(energy, "%.6g\t", e_s);
				}

				e[s] = energy_density_cutoff(e_min, e_s);

				u[s].ux = 0.0;		// zero initial velocity
				u[s].uy = 0.0;
			#ifndef BOOST_INVARIANT
				u[s].un = 0.0;
			#endif

				up[s].ux = 0.0;		// also set up = u
				up[s].uy = 0.0;
			#ifndef BOOST_INVARIANT
				up[s].un = 0.0;
			#endif
			}

			if(initial.trento_average_over_events)
			{
				if(j < ny + 1)
				{
					fprintf(energy, "\n");
				}
			}
		}

		if(initial.trento_average_over_events)
		{
			if(k < nz + 1)
			{
				fprintf(energy, "\n");
			}
		}

	}

	if(initial.trento_average_over_events)
	{
		fclose(energy);
	}

#ifdef OPENMP
   double t2 = omp_get_wtime();
   printf("Trento initial conditions took %lf seconds\n", t2 - t1);
#endif
}





