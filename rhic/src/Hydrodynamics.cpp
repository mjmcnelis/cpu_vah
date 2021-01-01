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


precision set_time_step(int n, precision t, precision dt_prev, precision t_next_output, lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro)
{
	precision dt_min = lattice.min_time_step;

	precision dt = lattice.fixed_time_step;                         // default time step is fixed

	if(lattice.adaptive_time_step)                                  // compute adaptive time step
	{
		if(n == 0)
		{
			dt = dt_min;
		}
		else
		{
			precision dt_CFL = compute_dt_CFL(t, lattice, hydro);

			precision dt_source = 1./0.;

			if(!hit_CFL_bound)                                      // skip dt_source calculation after hit CFL bound
			{
				int update = 0;

				// get total source function qI <= E
				euler_step(t, q, qI, e, lambda, aT, aL, up, u, dt_prev, dt_prev, lattice, hydro, update, hit_CFL_bound);

				// here Q holds previous q (from swap_hydro_variables)
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

	// further adjust dt to output hydro evolution at specific times (except Bjorken)
	if(hydro.run_hydro == 1 && initial.initial_condition_type != 1)
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
	print_parameters(lattice, hydro);
	allocate_memory(lattice);                                    // grid allocation

	precision t = hydro.tau_initial;                             // starting longtudinal proper time
	precision dt = lattice.min_time_step;                        // initial time step

	if(!lattice.adaptive_time_step)
	{
		dt = lattice.fixed_time_step;
	}

	precision dt_prev = dt;                                      // previous time step

	set_initial_conditions(t, lattice, initial, hydro, trento);  // initialize (q, e, u, up)
	set_ghost_cells(q, e, u, lattice);                           // for (q, e, u)

	precision t_out = t;                                         // track hydro evolution output times
	precision dt_out = lattice.output_interval;
  	int number_outputs = 0;

	freezeout_finder finder(lattice, hydro);                     // initialize freezeout finder

	if(hydro.run_hydro == 2)
	{
		finder.load_initial_grid(t, q, e, u);
	}

	int grid_below_Tswitch = 0;                                  // number of times freezeout finder searches a grid below Tswitch
	int freezeout_depth = 3;                                     // stop hydro evolution once grid_below_Tswitch = freezeout_depth

	int steps = 0;

#ifdef _OPENMP
  	double t1 = omp_get_wtime();
#else
  	clock_t start = clock();
#endif
	// fluid dynamic evolution
	for(int n = 0; n < lattice.max_time_steps; n++)
	{
		dt = set_time_step(n, t, dt_prev, t_out + dt_out, lattice, initial, hydro);

		if(hydro.run_hydro == 1)                                    // output hydro evolution
		{
		#ifdef FREEZEOUT_SLICE
			output_freezeout_slice_x(t, lattice, hydro);        // output freezeout surface slices
		#ifndef BOOST_INVARIANT
			output_freezeout_slice_z(t, lattice, hydro);
		#endif
		#endif
			if(n == 0 || initial.initial_condition_type == 1)   // output first step (or every time step if Bjorken)
			{
				long cells_above_Tswitch = number_of_cells_above_freezeout_temperature(lattice, hydro);
				print_hydro_center(n, t, lattice, hydro, cells_above_Tswitch);

				number_outputs++;
				output_hydro_simulation(t, dt_prev, lattice, initial, hydro);

				if(all_cells_below_freezeout_temperature(lattice, hydro))
				{
					break;                              // stop hydro evolution
				}
			}
			else if(fabs(t - t_out - dt_out) < dt_eps)          // output at regular time intervals
			{
				long cells_above_Tswitch = number_of_cells_above_freezeout_temperature(lattice, hydro);
				print_hydro_center(n, t, lattice, hydro, cells_above_Tswitch);

				number_outputs++;
				output_hydro_simulation(t, dt_prev, lattice, initial, hydro);

				if(all_cells_below_freezeout_temperature(lattice, hydro))
				{
				#ifdef FREEZEOUT_SLICE
					if(t > 17.0)                        // for freezeout slice plot only
					{
						printf("Ending hydro simulation at t = %lf fm/c for freezeout slice plot\n", t);
						break;
					}
				#endif
					break;                              // stop hydro evolution
				}

				t_out += dt_out;
			}
		}
		else if(hydro.run_hydro == 2)                               // construct particlization hypersurface
		{
			if(n == 0)
			{
				long cells_above_Tswitch = number_of_cells_above_freezeout_temperature(lattice, hydro);
				print_hydro_center(n, t, lattice, hydro, cells_above_Tswitch);

			}
			else if(n % hydro.freezeout_finder_period == 0)         // search for freezeout cells
			{
				long cells_above_Tswitch = number_of_cells_above_freezeout_temperature(lattice, hydro);
				print_hydro_center(n, t, lattice, hydro, cells_above_Tswitch);

				finder.load_current_grid(q, e, u);

			#ifdef BOOST_INVARIANT
				finder.find_2d_freezeout_cells(t, hydro);
			#else
				finder.find_3d_freezeout_cells(t, hydro);
			#endif
				if(all_cells_below_freezeout_temperature(lattice, hydro))
				{
					grid_below_Tswitch++;
				}
			}

			if(grid_below_Tswitch >= freezeout_depth)
			{
				break;		// stop hydro evolution
			}
		}

		int update = 1;                                             // RK2 iteration:

		evolve_hydro_one_time_step(n, t, dt, dt_prev, lattice, hydro, update, hit_CFL_bound);

		t += dt;

		dt_prev = dt;

		steps++;

		if(number_outputs == 3)
		{
			dt_out = 2. * lattice.output_interval;		// double output interval after three outputs
		}
	}

	if(steps >= lattice.max_time_steps)
	{
		printf("run_hydro error: exceeded max number of time steps = %d (simulation failed to finish)\n", steps);
		exit(-1);
	}

#ifdef _OPENMP
  	double t2 = omp_get_wtime();
  	double duration = t2 - t1;
#else
  	double duration = (clock() - start) / (double)CLOCKS_PER_SEC;
#endif
	print_run_time(t, duration, (double)steps, lattice, sample);        // output runtime benchmarks

	free_memory();                                                      // deallocate grid

	if(hydro.run_hydro == 2)
	{
		finder.free_finder_memory(sample);                          // deallocate grid stack in freezeout finder (not freezeout surface)
	}

	printf("\nFinished hydro evolution\n");

	return finder.surface;
}




