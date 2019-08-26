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
#include "../include/GhostCells.h"
#include "../include/InitialConditions.h"
#include "../include/KurganovTadmor.h"
#include "../include/AdaptiveTimeStep.h"
using namespace std;

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
	precision e_switch = equilibriumEnergyDensity(T_switch / hbarc, hydro.conformal_eos_prefactor);

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
	return true;
}



void test_output(lattice_parameters lattice, hydro_parameters hydro)
{
	precision dt_start = lattice.fixed_time_step;
	if(lattice.adaptive_time_step) dt_start = lattice.min_time_step;

	precision dt = dt_start;
	precision dt_out = lattice.output_interval;
	precision dt_min = lattice.min_time_step;

	precision t = hydro.tau_initial;
	precision t_out = t;

	cout << "\n0\t" << t << "\t" << dt << "\t" << t_out + dt_out << endl;

	print_hydro_center(0, t, lattice, hydro);

	for(int n = 0; n <= 200; n++)
	{	
		dt = dt_start;

		if(t + dt_eps < t_out + dt_out)
		{
			if(t + dt > t_out + dt_out)
			{
				dt = fmax(dt_min, t_out + dt_out - t);		// should be incorporated in next time step
			}
		}
		
		cout << n << "\t" << t << "\t" << dt << "\t" << t_out + dt_out << endl;

		if(fabs(t - t_out - dt_out) < dt_eps)
		{
			//print_hydro_center(n, t, lattice, hydro);
			t_out += dt_out;
		}
		
		t += dt;
		
	}
	exit(-1);
}


precision set_time_step(int n, precision t_current, precision dt_prev, precision t_next_output, lattice_parameters lattice)
{
	if(n == 0) return dt_prev;


}





void run_hydro(lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro)
{
	precision dt_start = lattice.fixed_time_step;		// starting time step

	if(lattice.adaptive_time_step) dt_start = lattice.min_time_step;

	precision dt = dt_start;							// initialize time step
	precision dt_prev = dt;								// previous time step

	print_parameters(lattice, hydro);
	allocate_memory(lattice);

	//test_output(lattice, hydro);

	precision t = hydro.tau_initial;					// initial longitudinal proper time
	set_initial_conditions(t, lattice, initial, hydro);	// initial conditions for (q, e, u)
	set_ghost_cells(q, e, u, lattice);					// initialize ghost cells in (q, e, u)


	precision t_out = t;
	precision dt_out = lattice.output_interval;

	double steps = 0;
	clock_t start = clock();

	print_hydro_center(0, t, lattice, hydro);

	for(int n = 0; n <= lattice.max_time_steps; n++)	// fluid dynamic evolution
	{
		dt = set_time_step(n, t, dt_prev, t_out + dt_out, lattice);

		if(fabs(t - t_out - dt_out) < dt_eps)
		{
			print_hydro_center(n, t, lattice, hydro);

			if(hydro.run_hydro == 1) output_dynamical_variables(t, dt_prev, lattice, initial, hydro);


			if(all_cells_below_freezeout_temperature(lattice, hydro)) 	// replace with freezeout finder
			{
				printf("\nAll cells below freezeout temperature\n\n"); 
				break;
			}


			t_out += dt_out;
		}

		evolve_hydro_one_time_step(t, dt, dt_prev, lattice, hydro);

		t += dt;
		dt_prev = dt;

		// if(lattice.adaptive_time_step)
		// {
		// 	precision dt_min = lattice.min_time_step;
		// 	// compute the smallest time scale in the fluid
		// 	//hydro_time_scales dt_hydro = compute_hydro_time_scales(t, q, e, u, up, nx, ny, nz, dt, dt_prev, dx, dy, dz, etabar_const, dt_min);
		// 	//dt = set_time_step(dt_hydro, dt_min);
		// }
		steps++;
	}

	double duration = (clock() - start) / (double)CLOCKS_PER_SEC;
	print_run_time(duration, steps, lattice);

	free_memory();

	printf("\nFinished hydro\n");
}





