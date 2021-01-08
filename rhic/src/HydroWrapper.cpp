#include <iostream>
#include <sstream>
#include <string>
#include "../include/HydroWrapper.h"
#include "../include/EnergyVector.h"
#include "../include/Parameters.h"
#include "../include/Print.h"
#include "../include/Hydrodynamics.h"
#include "../include/Output.h"
#include "../include/OpenMP.h"

using namespace std;

// this is a C++ wrapper for OSU hydro codes (cpu-vh, gpu-vh, cpu_vah, etc)
// originally written by Derek Everett 2018


HYDRO::HYDRO()
{
#ifndef JETSCAPE
	freezeout_surface_file.open("output/surface.dat");			// for writing freezeout surface to file
#endif
}



HYDRO::~HYDRO()
{

}



void HYDRO::read_trento_energy_density_profile_from_memory(std::vector<double> trento_energy)
{
    // read initial energy density profile from TRENTo (for JETSCAPE, initial_condition_type = 5)
    // note: units should be [GeV/fm^3], later converted to [fm^-4]

	trento_energy_density_profile = trento_energy;
}



void HYDRO::store_freezeout_surface(freezeout_surface surface)
{
#ifdef JETSCAPE
	printf("\nStoring freezeout surface in memory for particlization module...\n\n");
#else
	printf("\nWriting freezeout surface to file...\n\n");
#endif

	long total_cells = surface.tau.size();

	if(total_cells == 0)
	{
		printf("HYDRO::store_freezeout_surface flag: freezeout surface is empty\n");
		return;
	}
	else
	{
		printf("Number of freezeout cells = %ld\n", total_cells);

		for(long i = 0; i < total_cells; i++)
		{
			// freezeout surface vectors for JETSCAPE are in physical units (i.e. energy density GeV/fm^3, etc)
		#ifdef JETSCAPE
			tau.push_back((double)surface.tau[i]);
			x.push_back(  (double)surface.x[i]);
			y.push_back(  (double)surface.y[i]);
			eta.push_back((double)surface.eta[i]);

			dsigma_tau.push_back((double)surface.dsigma_tau[i]);
			dsigma_x.push_back(  (double)surface.dsigma_x[i]);
			dsigma_y.push_back(  (double)surface.dsigma_y[i]);
			dsigma_eta.push_back((double)surface.dsigma_eta[i]);

			ux.push_back((double)surface.ux[i]);
			uy.push_back((double)surface.uy[i]);
			un.push_back((double)surface.un[i]);

			E.push_back((double)surface.E[i]);
			T.push_back((double)surface.T[i]);
			P.push_back((double)surface.P[i]);

			pixx.push_back((double)surface.pixx[i]);
			pixy.push_back((double)surface.pixy[i]);
			pixn.push_back((double)surface.pixn[i]);
			piyy.push_back((double)surface.piyy[i]);
			piyn.push_back((double)surface.piyn[i]);

			Pi.push_back((double)surface.Pi[i]);


			// freezeout surface file in hbarc = 1 units (i.e. energy density in fm^-4, etc)
			// these units will be undone in iS3D, which reads in this file as input
		#else
			freezeout_surface_file
			<< surface.tau[i]          << " " << surface.x[i]            << " " << surface.y[i]            << " " << surface.eta[i]        << " "
			<< surface.dsigma_tau[i]   << " " << surface.dsigma_x[i]     << " " << surface.dsigma_y[i]     << " " << surface.dsigma_eta[i] << " "
			<< surface.ux[i]           << " " << surface.uy[i]           << " " << surface.un[i]           << " "
			<< surface.E[i] / hbarc    << " " << surface.T[i] / hbarc    << " " << surface.P[i] / hbarc    << " "
			<< surface.pixx[i] / hbarc << " " << surface.pixy[i] / hbarc << " " << surface.pixn[i] / hbarc << " "
			<< surface.piyy[i] / hbarc << " " << surface.piyn[i] / hbarc << " " << surface.Pi[i] / hbarc   << endl;

		#endif

			// sims script: should surface.dat be empty if event fails?
		}
	}
#ifndef JETSCAPE
	freezeout_surface_file.close();
#endif
}



void HYDRO::free_freezeout_surface()
{
	printf("\nFreeing freezeout surface in hydro module...\n\n");

	tau.clear();
	x.clear();
	y.clear();
	eta.clear();

	dsigma_tau.clear();
	dsigma_x.clear();
	dsigma_y.clear();
	dsigma_eta.clear();

	ux.clear();
	uy.clear();
	un.clear();

	E.clear();
	T.clear();
	P.clear();

	pixx.clear();
	pixy.clear();
	pixn.clear();
	piyy.clear();
	piyn.clear();

	Pi.clear();
}



void HYDRO::start_hydro(int argc, char **argv)
{
#ifdef _OPENMP
	printf("Running CPU-VAH with OpenMP acceleration: number of threads = %d\n\n", omp_get_max_threads());
#endif

	bool sample_parameters = false;     // default: use model parameters in parameters/
	int sample = 0;
	random_model_parameters random;

#ifdef RANDOM_MODEL_PARAMETERS
	if(argc == 2)                       // read sample model parameters
	{
		sample_parameters = true;

		istringstream iss(argv[1]);

        if(!(iss >> sample))
        {
  			printf("HYDRO::start_hydro error: parameter sample index %s invalid\n", argv[1]);
  			exit(-1);
  		}
		else if(!iss.eof())
		{
			printf("HYDRO::start_hydro error: trailing characeters after sample index %s.\n", argv[1]);
			exit(-1);
		}
		if(sample <= 0)
		{
			printf("HYDRO::start_hydro error: parameter sample index %d must be greater than zero\n", sample);
			exit(-1);
		}

		random = load_random_model_parameters(sample);
	}
	else if(argc > 2)
	{
		printf("HYDRO::start_hydro error: too many command line arguments (max is 1)\n");
		exit(-1);
	}
#endif

	hydro_parameters hydro = load_hydro_parameters(sample_parameters, random);
	initial_condition_parameters initial = load_initial_condition_parameters(sample_parameters, random);
	lattice_parameters lattice = load_lattice_parameters(hydro, initial, sample_parameters, sample);

	if(hydro.run_hydro)                 // main hydro simulation (freezeout surface empty by default)
	{
		print_hydro_mode(hydro);
		store_freezeout_surface(run_hydro(lattice, initial, hydro, sample, trento_energy_density_profile));
	}

	printf("\nEnd of hydro simulation\n");
}



void HYDRO::start_hydro_no_arguments()
{
#ifdef _OPENMP
	printf("Running CPU-VAH with OpenMP acceleration: number of threads = %d\n\n", omp_get_max_threads());
#endif

	bool sample_parameters = false;     // will read default parameters in parameters/
	int sample = 0;                     // can't use auto_grid (need different way to read in sample model parameters)
	random_model_parameters random;

	hydro_parameters hydro = load_hydro_parameters(sample_parameters, random);
	initial_condition_parameters initial = load_initial_condition_parameters(sample_parameters, random);
	lattice_parameters lattice = load_lattice_parameters(hydro, initial, sample_parameters, sample);

	if(hydro.run_hydro)                 // main hydro simulation (freezeout surface empty by default)
	{
		print_hydro_mode(hydro);
		store_freezeout_surface(run_hydro(lattice, initial, hydro, sample, trento_energy_density_profile));
	}

	printf("\nEnd of hydro simulation\n");
}




