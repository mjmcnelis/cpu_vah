#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <math.h>
#include "../include/Hydrodynamics.h"
#include "../include/Print.h"
#include "../include/Output.h"
#include "../include/Parameters.h"
#include "../include/Precision.h"
#include "../include/Macros.h"
#include "../include/GhostCells.h"
#include "../include/InitialConditions.h"
#include "../include/KurganovTadmor.h"
#include "../include/AdaptiveTimeStep.h"
#include "../include/FreezeoutFinder.h"
#include "../include/OpenMP.h"

using namespace std;

bool hit_CFL_bound = false;
bool after_output = false;
precision dt_after_output;

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
	precision e_switch = equilibrium_energy_density_new(T_switch / hbarc, hydro.conformal_eos_prefactor);

    for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				if(e[s] >= e_switch)
				{
                    return false;
				}
			}

		}
    }
	return true;
}

long number_of_cells_above_freezeout_temperature(lattice_parameters lattice, hydro_parameters hydro)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision T_switch = hydro.freezeout_temperature_GeV;
	precision e_switch = equilibrium_energy_density_new(T_switch / hbarc, hydro.conformal_eos_prefactor);

	long cells = 0;

	// add omp later
	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				if(e[s] >= e_switch)
				{
					cells++;
				}
			}
		}
	}
	return cells;
}


precision set_the_time_step(int n, precision t, precision dt_prev, precision t_next_output, lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro)
{
	precision dt_fix = lattice.fixed_time_step;
	precision dt_min = lattice.min_time_step;

	precision dt = dt_fix;                                          // default time step is fixed

	if(lattice.adaptive_time_step)                                  // adaptive time step
	{
		if(n == 0)
		{
			dt = dt_min;
		}
		else
		{
			precision dt_CFL = compute_dt_CFL(t, lattice, hydro);   // less strict CFL condition dx / (8.ax)

			if(lattice.adaptive_time_step == 2)
			{
				dt_CFL = dt_fix;                                    // strict CFL condition <= dx / 8 (lattice.fixed_time_step <= dx / 8)
			}

			precision dt_source = 1./0.;

			if(!hit_CFL_bound)      // assumes that dt_source remains > dt_CFL after hit bound
			{                       // (may not be always true, but it's convenient to skip this step eventually)

				euler_step(t, q, qI, e, lambda, aT, aL, up, u, dt_prev, dt_prev, lattice, hydro, 0, 0);    // compute source function

				dt_source = compute_dt_source(t, Q, q, qI, dt_prev, lattice);

				if(dt_source >= dt_CFL)
				{
					printf("\nHit CFL bound at t = %lf\n\n", t);
					hit_CFL_bound = true;
				}
			}

			dt = compute_adaptive_time_step(t, dt_CFL, dt_source, dt_min);
		}
	}

	if(hydro.run_hydro == 1 && initial.initial_condition_type != 1) // adjust dt further (for timed hydro outputs, except Bjorken)
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


freezeout_surface run_hydro(lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro, int sample, std::vector<double> trento)
{
	precision dt_start = lattice.fixed_time_step;       // starting time step

	if(lattice.adaptive_time_step)
	{
		dt_start = lattice.min_time_step;
	}

	precision dt = dt_start;                            // initialize time step
	precision dt_prev = dt;                             // previous time step

	print_parameters(lattice, hydro);
	allocate_memory(lattice);

	long nx = lattice.lattice_points_x;
	long ny = lattice.lattice_points_y;
	long nz = lattice.lattice_points_eta;
	long grid_size = nx * ny * nz;

	precision t = hydro.tau_initial;                             // initial longitudinal proper time

	set_initial_conditions(t, lattice, initial, hydro, trento);  // initial conditions for (q, e, u)

	set_ghost_cells(q, e, u, lattice);                           // initialize ghost cells in (q, e, u)

	precision t_out = t;                                         // output times
	precision dt_out = lattice.output_interval;

	bool double_dt_out = true;									 // double output interval after three outputs
  	int number_outputs = 0;                                      // set to false if want regular outputs

#ifdef BOOST_INVARIANT
	printf("Running 2+1d hydro simulation...\n\n");
#else
	printf("Running 3+1d hydro simulation...\n\n");
#endif

	freezeout_finder fo_finder(lattice, hydro);         // freezeout finder class
	int freezeout_period = lattice.tau_coarse_factor;   // time steps between freezeout finder calls
	int grid_below_Tswitch = 0;                         // number of times freezeout finder searches a grid below Tswitch
	int freezeout_depth = 3;                            // max number of time steps freezeout finder goes below Tswitch

	int steps = 0;

#ifdef OPENMP
  	double t1 = omp_get_wtime();                        // is this the right omp time I want?
    printf("Staring omp time = %lf\n", t1);
#else
  	clock_t start = clock();
#endif

	// fluid dynamic evolution
	//----------------------------------------------------------
	for(int n = 0; n < lattice.max_time_steps; n++)
	{
		dt = set_the_time_step(n, t, dt_prev, t_out + dt_out, lattice, initial, hydro);

		if(hydro.run_hydro == 1)                                // outputs hydro data at regular time intervals (if hydro.output = 1)
		{
		#ifdef FREEZEOUT_SLICE
			output_freezeout_slice_x(t, lattice, hydro);		// output freezeout surface slices
		#ifndef BOOST_INVARIANT
			output_freezeout_slice_z(t, lattice, hydro);
		#endif
		#endif

			if(n == 0 || initial.initial_condition_type == 1)   // output first time or output bjorken at every time step
			{
			#ifdef PRINT_HYDRO
				long cells_above_Tswitch = number_of_cells_above_freezeout_temperature(lattice, hydro);
				print_hydro_center(n, t, lattice, hydro, cells_above_Tswitch);
			#endif

				if(hydro.output)
				{
					number_outputs++;
					output_dynamical_variables(t, dt_prev, lattice, initial, hydro);
				}

				if(all_cells_below_freezeout_temperature(lattice, hydro))
				{
					break;
				}
			}
			else if(fabs(t - t_out - dt_out) < dt_eps) 			// output hydrodynamic quantities at regular time intervals
			{
			#ifdef PRINT_HYDRO
				long cells_above_Tswitch = number_of_cells_above_freezeout_temperature(lattice, hydro);
				print_hydro_center(n, t, lattice, hydro, cells_above_Tswitch);
			#endif

				if(hydro.output)
				{
					number_outputs++;
					output_dynamical_variables(t, dt_prev, lattice, initial, hydro);
				}

				if(all_cells_below_freezeout_temperature(lattice, hydro))
				{
					break;
				}

				t_out += dt_out;
			}
		}
		else if(hydro.run_hydro == 2)                           // construct freezeout surface
		{
			if(n == 0)                                          // initialize hydro info in freezeout finder
			{
			#ifdef PRINT_HYDRO
				long cells_above_Tswitch = number_of_cells_above_freezeout_temperature(lattice, hydro);
				print_hydro_center(n, t, lattice, hydro, cells_above_Tswitch);
			#endif

				fo_finder.set_hydro_evolution(t, q, e, u);
			}
			else if(n % freezeout_period == 0)                  // find freezeout cells
			{
			#ifdef PRINT_HYDRO
				long cells_above_Tswitch = number_of_cells_above_freezeout_temperature(lattice, hydro);
				print_hydro_center(n, t, lattice, hydro, cells_above_Tswitch);
			#endif

				fo_finder.swap_and_set_hydro_evolution(q, e, u);

			#ifdef BOOST_INVARIANT
				fo_finder.find_2d_freezeout_cells(t, hydro);
			#else
				fo_finder.find_3d_freezeout_cells(t, hydro);
			#endif

				if(all_cells_below_freezeout_temperature(lattice, hydro))
				{
					grid_below_Tswitch++;
					printf("\nNumber of times grid went below freezeout temperature during freezeout finder call: %d\n\n", grid_below_Tswitch);
				}
			}

			if(grid_below_Tswitch >= freezeout_depth)
			{
				break;
			}
		}
		else if(hydro.run_hydro == 3)                           // print center every PRINT_PERIOD step
		{
			if(n % PRINT_PERIOD == 0)
			{
			#ifdef PRINT_HYDRO
				long cells_above_Tswitch = number_of_cells_above_freezeout_temperature(lattice, hydro);
				print_hydro_center(n, t, lattice, hydro, cells_above_Tswitch);
			#endif
                if(all_cells_below_freezeout_temperature(lattice, hydro))
                {
                    break;
                }
			}
		}

		int update = 1;

		evolve_hydro_one_time_step(n, t, dt, dt_prev, lattice, hydro, update, hit_CFL_bound);

		t += dt;

		dt_prev = dt;

		steps++;

		if(number_outputs == 3 && double_dt_out)
		{
			dt_out = 2. * lattice.output_interval;
		}
	}
	//----------------------------------------------------------

	if(steps >= lattice.max_time_steps)
	{
		printf("run_hydro error: exceeded max number of time steps = %d (hydro simulation failed)\n", steps);
		exit(-1);
	}

#ifdef OPENMP
  	double t2 = omp_get_wtime();
    printf("End omp time = %lf\n", t2);
  	double duration = t2 - t1;
#else
  	double duration = (clock() - start) / (double)CLOCKS_PER_SEC;  // hydro evolution runtime in seconds (ignores initalization runtime)
#endif

	print_run_time(t, duration, (double)steps, lattice, sample);

	free_memory();

	if(hydro.run_hydro == 2)
	{
		fo_finder.free_finder_memory(sample);
		printf("\nFinished hydro evolution\n");
		return fo_finder.surface;                   // return freezeout surface
	}

	printf("\nFinished hydro evolution\n");

	freezeout_surface surface;                      // return an empty surface (default)

	return surface;
}




