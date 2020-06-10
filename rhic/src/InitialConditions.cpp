#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include "../include/Macros.h"
#include "../include/Hydrodynamics.h"
#include "../include/DynamicalVariables.h"
#include "../include/Precision.h"
#include "../include/Parameters.h"
#include "../include/Trento.h"
#include "../include/ViscousBjorken.h"
#include "../include/AnisoBjorken.h"
#include "../include/IdealGubser.h"
#include "../include/ViscousGubser.h"
#include "../include/AnisoGubser.h"
#include "../include/EquationOfState.h"
#include "../include/Projections.h"
#include "../include/AnisoVariables.h"
#include "../include/Viscosities.h"
using namespace std;

#define THETA_FUNCTION(X) ((double)X < (double)0 ? (double)0 : (double)1)


inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


void set_initial_T_taumu_variables(double t, int nx, int ny, int nz, hydro_parameters hydro)
{
	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				precision e_s = e[s];

			#ifdef E_CHECK
				q[s].e_check = e_s;
			#endif

				precision ux = u[s].ux;
				precision uy = u[s].uy;
			#ifndef BOOST_INVARIANT
				precision un = u[s].un;
			#else
				precision un = 0;
			#endif

				precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t * t * un * un);

			#ifdef ANISO_HYDRO
			#ifndef BOOST_INVARIANT
				precision utperp = sqrt(1.  +  ux * ux  +  uy * uy);
				precision zt = t * un / utperp;
				precision zn = ut / t / utperp;
			#else
				precision zt = 0;
				precision zn = 1. / t;
			#endif
				precision pl = q[s].pl;
				precision pt = q[s].pt;
			#else
				equation_of_state_new eos(e_s, hydro.conformal_eos_prefactor);
				precision p = eos.equilibrium_pressure();
			#endif

			#ifdef PIMUNU
				precision pitt = q[s].pitt;
				precision pitx = q[s].pitx;
				precision pity = q[s].pity;
			#ifndef BOOST_INVARIANT
				precision pitn = q[s].pitn;
			#else
				precision pitn = 0;
			#endif
			#else
				precision pitt = 0;
				precision pitx = 0;
				precision pity = 0;
				precision pitn = 0;
			#endif

			#ifdef WTZMU
				precision WtTz = q[s].WtTz;
				precision WxTz = q[s].WxTz;
				precision WyTz = q[s].WyTz;
				precision WnTz = q[s].WnTz;
			#else
				precision WtTz = 0;
				precision WxTz = 0;
				precision WyTz = 0;
				precision WnTz = 0;
			#endif

			#ifdef PI
				precision Pi = q[s].Pi;
			#else
				precision Pi = 0;
			#endif

			#ifdef ANISO_HYDRO
				q[s].ttt = (e_s + pt) * ut * ut  -   pt  +  (pl - pt) * zt * zt  +  2. * WtTz * zt  +  pitt;
				q[s].ttx = (e_s + pt) * ut * ux  +  WxTz * zt  +  pitx;
				q[s].tty = (e_s + pt) * ut * uy  +  WyTz * zt  +  pity;
			#ifndef BOOST_INVARIANT
				q[s].ttn = (e_s + pt) * ut * un  +  (pl - pt) * zt * zn  +  WtTz * zn  +  WnTz * zt  +  pitn;
			#endif
			#else
				q[s].ttt = (e_s + p + Pi) * ut * ut  -   p  -  Pi  +  pitt;
				q[s].ttx = (e_s + p + Pi) * ut * ux  +  pitx;
				q[s].tty = (e_s + p + Pi) * ut * uy  +  pity;
			#ifndef BOOST_INVARIANT
				q[s].ttn = (e_s + p + Pi) * ut * un  +  pitn;
			#endif
			#endif
			}
		}
	}
}


void set_anisotropic_initial_condition(int nx, int ny, int nz, hydro_parameters hydro)
{
	precision plpt_ratio = hydro.plpt_ratio_initial;
	precision conformal_prefactor = hydro.conformal_eos_prefactor;
	precision t = hydro.tau_initial;

	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				precision e_s = e[s];

				equation_of_state_new eos(e_s, conformal_prefactor);
				precision p = eos.equilibrium_pressure();

				precision pl = 3. * p * plpt_ratio / (2. + plpt_ratio);
				precision pt = 3. * p / (2. + plpt_ratio);

			#ifdef ANISO_HYDRO						// anisotropic hydro initial conditions
				q[s].pl = pl;
				q[s].pt = pt;

			#ifdef LATTICE_QCD						// initialize mean field and anisotropic variables
				precision T = eos.T;
				precision b = eos.equilibrium_mean_field();
				precision mass = T * eos.z_quasi();

				precision lambda0 = T;				// initial guess for anisotropic variables
				precision aT0 = 1.;
				precision aL0 = 1.;

				aniso_variables X = find_anisotropic_variables(e_s, pl, pt, b, mass, lambda0, aT0, aL0);

				if(X.did_not_find_solution)
				{
					aniso_regulation[s] += 1;
				}

				q[s].b = b;
				lambda[s] = X.lambda;
				aT[s] = X.aT;
				aL[s] = X.aL;
			#endif

			#ifdef PIMUNU
		  		q[s].pitt = 0;
		  		q[s].pitx = 0;
		  		q[s].pity = 0;
		  		q[s].pixx = 0;
		  		q[s].pixy = 0;
		  		q[s].piyy = 0;
		  	#ifndef BOOST_INVARIANT
		  		q[s].pitn = 0;
		  		q[s].pixn = 0;
		  		q[s].piyn = 0;
		  		q[s].pinn = 0;
		  	#endif
			#endif

		  	#ifdef WTZMU
		  		q[s].WtTz = 0;
		  		q[s].WxTz = 0;
		  		q[s].WyTz = 0;
		  		q[s].WnTz = 0;
			#endif

			#else 									// viscous hydro initial conditions

		 	#ifdef PIMUNU
		  		precision ux = u[s].ux;
		  		precision uy = u[s].uy;
		  		precision un = 0;
		  	#ifndef BOOST_INVARIANT
		  		un = u[s].un;
		  	#endif
		  		precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t * t * un * un);
		  		precision utperp = sqrt(1.  +  ux * ux  +  uy * uy);
		  		precision zt = t * un / utperp;
		  		precision zn = ut / t / utperp;

		  		spatial_projection Delta(ut, ux, uy, un, t * t);

		  		// pi^\munu = (pl - pt)/3 . (\Delta^\munu + 3.z^\mu.z^\nu)
		  		q[s].pitt = (pl - pt)/3. * (Delta.Dtt  +  3. * zt * zt);
		  		q[s].pitx = (pl - pt)/3. * Delta.Dtx;
		  		q[s].pity = (pl - pt)/3. * Delta.Dty;
		  		q[s].pixx = (pl - pt)/3. * Delta.Dxx;
		  		q[s].pixy = (pl - pt)/3. * Delta.Dxy;
		  		q[s].piyy = (pl - pt)/3. * Delta.Dyy;
		  		q[s].pinn = (pl - pt)/3. * (Delta.Dnn  +  3. * zn * zn);

		  	#ifndef BOOST_INVARIANT
		  		q[s].pitn = (pl - pt)/3. * (Delta.Dtn  +  3. * zt * zn);
		  		q[s].pixn = (pl - pt)/3. * Delta.Dxn;
		  		q[s].piyn = (pl - pt)/3. * Delta.Dyn;
		  	#endif
		  	#endif

		  	#ifdef PI
		  		q[s].Pi = (pl + 2.*pt)/3.  -  p;
		  	#endif

		  	#endif
			}
		}
	}
}


void longitudinal_energy_density_profile(double * const __restrict__ eL, int nz, double dz, initial_condition_parameters initial)
{
	double etaFlat = initial.rapidity_mean;
	double etaVariance = initial.rapidity_variance;

	// profile along eta direction is a smooth plateu that exponentially decays when |eta| > etaFlat
	for(int k = 0; k < nz; k++)
	{
		double eta = (k - (nz - 1.)/2.) * dz;				// physical eta points
		double etaScaled = fabs(eta)  -  0.5 * etaFlat;

		eL[k] = exp(- 0.5 * etaScaled * etaScaled / etaVariance * THETA_FUNCTION(etaScaled));
	}
}

/*
// Optical or MC glauber initial energy density profile
void set_Glauber_energy_density_and_flow_profile(int nx, int ny, int nz, double dx, double dy, double dz, initial_condition_parameters initial, hydro_parameters hydro)
{
	int initial_condition_type = initial.initial_condition_type;

	double T0 = initial.initialCentralTemperatureGeV;									// central temperature (GeV)
	//double e0 = equilibrium_energy_density(T0 / hbarc, hydro.conformal_eos_prefactor);	// energy density scale factor
	double e0 = equilibrium_energy_density_new(T0 / hbarc, hydro.conformal_eos_prefactor);	// energy density scale factor


	precision e_min = hydro.energy_min;

	double eT[nx * ny];		// normalized transverse profile
	double eL[nz];			// normalized longitudinal profile

	// compute normalized transverse and longitudinal profiles
	if(initial_condition_type == 4)
	{
		Optical_Glauber_energy_density_transverse_profile(eT, nx, ny, dx, dy, initial);
	}
	else if(initial_condition_type == 5)
	{
		MC_Glauber_energy_density_transverse_profile(eT, nx, ny, dx, dy, initial);
	}

	longitudinal_energy_density_profile(eL, nz, dz, initial);

	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				precision e_s = e0 * eT[i - 2 + (j - 2) * nx] * eL[k - 2];

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
		}
	}
}
*/

void read_block_energy_density_file(int nx, int ny, int nz, hydro_parameters hydro)
{
	// load block energy density file to energy density e[s]
	// and set u, up components to 0 (assumes no flow profile)

	// note: e_block.dat in block format (nested matrix is nz x (ny x nx))
	// and contains a header with grid points (nx,ny,nz)

	// e.g. 2d block format:
	// nx
	// ny
	// nz
	// e(1,1)	...	e(nx,1)
	// ...		...	...
	// e(1,ny)	... e(nx,ny)

	// e.g. 3d block format:
	// nx
	// ny
	// nz
	// e(1,1,1)		...		e(nx,1,1)
	// ...			...		...
	// e(1,ny,1)	... 	e(nx,ny,1)
	// e(1,1,2)		...		e(nx,1,2)
	// ...			...		...
	// e(1,ny,2)	... 	e(nx,ny,2)
	// ...
	// e(1,1,nz)	...		e(nx,1,nz)
	// ...			...		...
	// e(1,ny,nz)	... 	e(nx,ny,nz)


	// see how block file is constructed in set_trento_energy_density_and_flow_profile() in Trento.cpp


	FILE * e_block;
  	e_block = fopen("tables/e_block.dat", "r");

  	if(e_block == NULL)
  	{
  		printf("Error: couldn't open e_block.dat file\n");
  	}

  	int nx_block;
  	int ny_block;
  	int nz_block;

  	fscanf(e_block, "%d\n%d\n%d", &nx_block, &ny_block, &nz_block);

	if(nx != nx_block || ny != ny_block || nz != nz_block)
	{
		printf("read_block_energy_density_file error: hydro grid and block file dimensions are inconsistent\n");
		exit(-1);
	}

	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				precision e_s;

				fscanf(e_block, "%lf\t", &e_s);

				e[s] = energy_density_cutoff(hydro.energy_min, e_s);

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
		}
	}
	fclose(e_block);
}



void set_trento_energy_density_profile(int nx, int ny, int nz, hydro_parameters hydro, std::vector<double> trento)
{
	if(trento.size() == 0)
	{
		printf("set_trento_energy_density_profile error: trento energy density profile is empty (initial condition type only compatible with JETSCAPE)\n");
		exit(-1);
	}

	if((nx * ny * nz) != trento.size())
	{
		printf("set_trento_energy_density_profile error: physical grid points and trento energy density vector size are inconsistent\n");
		exit(-1);
	}

  	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);
		        int st = linear_column_index(i - 2, j - 2, k - 2, nx + 4, ny + 4);   // TRENTo vector has no ghost cells

		        e[s] = energy_density_cutoff(hydro.energy_min, trento[st] / hbarc);  // convert units to [fm^-4]

		        u[s].ux = 0;		// zero initial velocity
				u[s].uy = 0;
			#ifndef BOOST_INVARIANT
				u[s].un = 0;
			#endif

				up[s].ux = 0;		// also set up = u
				up[s].uy = 0;
			#ifndef BOOST_INVARIANT
				up[s].un = 0;
			#endif
			}
		}
	}
}


void set_initial_conditions(precision t, lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro, std::vector<double> trento)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_x;
	int nz = lattice.lattice_points_eta;

	double dx = lattice.lattice_spacing_x;
	double dy = lattice.lattice_spacing_y;
	double dz = lattice.lattice_spacing_eta;

	precision dt = lattice.fixed_time_step;

	if(lattice.adaptive_time_step) dt = lattice.min_time_step;

	printf("\nInitial conditions = ");
	switch(initial.initial_condition_type)
	{
		case 1:		// viscous hydro or anisotropic hydro bjorken
		{
			printf("Bjorken\n\n");

		#ifndef BOOST_INVARIANT
			printf("Bjorken initial condition error: define BOOST_INVARIANT in /rhic/include/Macros.h\n");
			exit(-1);
		#endif

		#ifdef ANISO_HYDRO
			printf("\nRunning semi-analytic anisotropic Bjorken solution...\n");
			run_semi_analytic_aniso_bjorken(lattice, initial, hydro);
			set_aniso_bjorken_initial_condition(nx, ny, nz, initial, hydro);
		#else
			printf("\nRunning semi-analytic viscous Bjorken solution...\n");
			run_semi_analytic_viscous_bjorken(lattice, initial, hydro);
			set_viscous_bjorken_initial_condition(nx, ny, nz, initial, hydro);
		#endif

			set_initial_T_taumu_variables(t, nx, ny, nz, hydro);

			break;
		}
		case 2:		// viscous hydro gubser
		{
		#ifdef PIMUNU
			printf("Viscous Gubser\n\n");
		#else
			printf("Ideal Gubser\n\n");
		#endif

		#ifndef BOOST_INVARIANT
			printf("Gubser initial condition error: define BOOST_INVARIANT in /rhic/include/Marcos.h..\n");
			exit(-1);
		#endif
		#ifdef ANISO_HYDRO
			printf("Gubser initial condition error: undefine ANISO_HYDRO in /rhic/include/Marcos.h..\n");
			exit(-1);
		#endif
		#ifndef CONFORMAL_EOS
			printf("\nGubser initial condition error: define CONFORMAL_EOS in /rhic/include/Marcos.h, exiting...\n");
			exit(-1);
		#endif

		#ifdef PIMUNU
			printf("Running semi-analytic viscous Gubser solution...\n\n");
			double T0_hat = run_semi_analytic_viscous_gubser(lattice, initial, hydro);
			set_viscous_gubser_initial_condition(T0_hat, nx, ny, nz, dt, dx, dy, dz, hydro, initial);
		#else
			printf("Running analytic ideal Gubser solution...\n\n");
			run_analytic_ideal_gubser(lattice, initial, hydro);
			set_ideal_gubser_initial_conditions(lattice, dt, initial, hydro);
		#endif

			set_initial_T_taumu_variables(t, nx, ny, nz, hydro);

			break;
		}
		case 3:		// anisotropic hydro gubser
		{
			printf("Anisotropic Gubser (residual shear stress initialized to zero)\n\n");

		#ifndef ANISO_HYDRO
			printf("Aniso Gubser error: ANISO_HYDRO not defined in /rhic/include/Marcos.h, exiting...\n");
			exit(-1);
		#endif

		#ifndef BOOST_INVARIANT
			printf("Aniso Gubser error: BOOST_INVARIANT not defined in /rhic/include/Marcos.h, exiting...\n");
			exit(-1);
		#endif

		#ifndef CONFORMAL_EOS
			printf("\nGubser initial condition error: CONFORMAL_EOS not defined in /rhic/include/Marcos.h, exiting...\n");
			exit(-1);
		#endif

			if(hydro.temperature_etas)
			{
				printf("Gubser initial condition error: temperature dependent eta/s breaks conformal invariance\n");
				exit(-1);
			}

			printf("Running semi-analytic anisotropic Gubser solution...\n\n");

			double T0_hat = run_semi_analytic_aniso_gubser(lattice, initial, hydro);

			set_aniso_gubser_energy_density_and_flow_profile(T0_hat, nx, ny, nz, dt, dx, dy, dz, hydro, initial);
			set_initial_T_taumu_variables(t, nx, ny, nz, hydro);

			break;
		}
		case 4:		// trento
		{
			printf("Trento (fluid velocity initialized to zero)\n\n");
			set_trento_energy_density_and_flow_profile(lattice, initial, hydro);
			set_anisotropic_initial_condition(nx, ny, nz, hydro);
			set_initial_T_taumu_variables(t, nx, ny, nz, hydro);

			break;
		}
		case 5:		// read custom energy density block file
		{
			printf("Reading custom energy density block file... (fluid velocity initialized to zero)\n\n");
			read_block_energy_density_file(nx, ny, nz, hydro);
			set_anisotropic_initial_condition(nx, ny, nz, hydro);
			set_initial_T_taumu_variables(t, nx, ny, nz, hydro);

			break;
		}
		case 6:		// read trento energy density profile from JETSCAPE C++ vector
		{
			printf("Reading trento energy density profile from JETSCAPE C++ vector... (fluid velocity initialized to zero)\n\n");
			set_trento_energy_density_profile(nx, ny, nz, hydro, trento);
			set_anisotropic_initial_condition(nx, ny, nz, hydro);
			set_initial_T_taumu_variables(t, nx, ny, nz, hydro);
			break;
		}
		default:
		{
			printf("\n\nset_initial_conditions error: initial condition type %d is not an option (see parameters/initial.properties)\n", initial.initial_condition_type);
			exit(-1);
		}
	}
}


