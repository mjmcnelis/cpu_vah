#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <ctime>
#include <math.h>
#include "../include/Hydrodynamics.h"
#include "../include/Print.h"
#include "../include/FileIO.h"
#include "../include/Parameters.h"
#include "../include/Precision.h"
#include "../include/Macros.h"
#include "../include/GhostCells.h"
#include "../include/InitialConditions.h"
#include "../include/KurganovTadmor.h"
#include "../include/AdaptiveTimeStep.h"
using namespace std;

bool hit_CFL_bound = false;

bool after_output = false;
precision dt_after_output;

const int freezeout_frequency = 10;
const precision dt_eps = 1.e-8;


inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


bool all_cells_below_freezeout_temperature(lattice_parameters lattice, hydro_parameters hydro)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision T_switch = hydro.freezeout_temperature_GeV;
	precision e_switch = equilibrium_energy_density(T_switch / hbarc, hydro.conformal_eos_prefactor);

	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				if(e[s] > e_switch) return false;
			}
		}
	}
	printf("\nAll cells below freezeout temperature\n\n");
	return true;
}


precision set_the_time_step(int n, precision t, precision dt_prev, precision t_next_output, lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro)
{
	precision dt_min = lattice.min_time_step;

	precision dt = lattice.fixed_time_step;		// default time step

	if(lattice.adaptive_time_step)				// adaptive time step
	{
		dt = dt_min;

		if(n > 0)
		{
			precision dt_CFL = compute_dt_CFL(t, lattice, hydro);

			precision dt_source = 1./0.;

			if(!hit_CFL_bound)		// assumes that dt_source remains > dt_CFL after hit bound
			{
				euler_step(t, q, qI, e, lambda, aT, aL, up, u, dt_prev, dt_prev, lattice, hydro, 0, 0);		// compute source function

				dt_source = compute_dt_source(t, Q, q, qI, dt_prev, lattice);

				if(dt_source >= dt_CFL) hit_CFL_bound = true;
			}

			dt = compute_adaptive_time_step(t, dt_CFL, dt_source, dt_min);
		}
	}
	if(hydro.run_hydro == 1 && initial.initialConditionType != 1)						// adjust dt further (for timed hydro outputs, except Bjorken)
	{
		if(after_output)
		{
			dt = dt_after_output;
			after_output = false;
		}
		else if(t + dt_eps < t_next_output)
		{
			if(t + dt > t_next_output)
			{
				dt = fmax(0.001 * dt_min, t_next_output - t);

				dt_after_output = dt_prev;

				after_output = true;
			}
		}
	}

	return dt;
}


void run_hydro(lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro)
{
	precision dt_start = lattice.fixed_time_step;		// starting time step

	if(lattice.adaptive_time_step) dt_start = lattice.min_time_step;

	precision dt = dt_start;							// initialize time step
	precision dt_prev = dt;								// previous time step

	print_parameters(lattice, hydro);
	allocate_memory(lattice);

	precision t = hydro.tau_initial;					// initial longitudinal proper time

	set_initial_conditions(t, lattice, initial, hydro);	// initial conditions for (q, e, u)

	set_ghost_cells(q, e, u, lattice);					// initialize ghost cells in (q, e, u)

	precision t_out = t;								// output times
	precision dt_out = lattice.output_interval;

#ifdef BOOST_INVARIANT
	printf("Running 2+1d hydro simulation...\n\n");
#else
	printf("Running 3+1d hydro simulation...\n\n");
#endif

	double steps = 0;
	clock_t start = clock();


	// fluid dynamic evolution
	//----------------------------------------------------------
	for(int n = 0; n <= lattice.max_time_steps; n++)
	{
		dt = set_the_time_step(n, t, dt_prev, t_out + dt_out, lattice, initial, hydro);

		if(hydro.run_hydro == 1)		// outputs hydro data at regular time intervals
		{
			if(n == 0 || initial.initialConditionType == 1)
			{
				print_hydro_center(n, t, lattice, hydro);
				output_dynamical_variables(t, dt_prev, lattice, initial, hydro);

				if(all_cells_below_freezeout_temperature(lattice, hydro)) break;
			}
			else if(fabs(t - t_out - dt_out) < dt_eps)
			{
				print_hydro_center(n, t, lattice, hydro);
				output_dynamical_variables(t, dt_prev, lattice, initial, hydro);

				if(all_cells_below_freezeout_temperature(lattice, hydro)) break;

				t_out += dt_out;
			}
		}
		else if(hydro.run_hydro == 2)	// finding the freezeout surface
		{
			if(n % freezeout_frequency == 0)
			{
				print_hydro_center(n, t, lattice, hydro);

				if(all_cells_below_freezeout_temperature(lattice, hydro)) break; 	// replace with freezeout finder
			}
		}

		int update = 1;

		evolve_hydro_one_time_step(n, t, dt, dt_prev, lattice, hydro, update, hit_CFL_bound);

		t += dt;

		dt_prev = dt;

		steps++;
	}
	//----------------------------------------------------------


	double duration = (clock() - start) / (double)CLOCKS_PER_SEC;
	print_run_time(duration, steps, lattice);

	free_memory();

	printf("\nFinished hydro\n");
}





