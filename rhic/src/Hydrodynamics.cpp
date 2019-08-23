#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <ctime>
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

#define FREQ 10


int central_index(int nx, int ny, int nz, int ncx, int ncy, int ncz)
{
	// (not sure what this means)
	int ictr = (nx % 2 == 0) ? ncx/2 : (ncx-1)/2;
	int jctr = (ny % 2 == 0) ? ncy/2 : (ncy-1)/2;
	int kctr = (nz % 2 == 0) ? ncz/2 : (ncz-1)/2;

	return ictr  +  ncx * (jctr  +  ncy * kctr);
}


void run_hydro(lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro)
{
	// hydro parameters
	precision t0 = hydro.tau_initial;
	precision etabar_const = hydro.shear_viscosity;
	precision T_switch = hydro.freezeout_temperature_GeV;


	// physical grid
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;
	int nt = lattice.max_number_of_time_steps;


	// computational grid = physical + ghost + white
	int ncx = nx + 4;
	int ncy = ny + 4;
	int ncz = nz + 4;


	// lattice spacing
	precision dx = lattice.lattice_spacing_x;
	precision dy = lattice.lattice_spacing_y;
	precision dz = lattice.lattice_spacing_eta;

	// option to turn on adaptive time step
	int adaptive_time_step = lattice.adaptive_time_step;
	precision dt_min = lattice.min_time_step;

	precision dt = lattice.fixed_time_step;	// default time step (fixed)
	if(adaptive_time_step)
	{
		dt = dt_min;
	}
	precision dt_prev = dt;		// previous time step for computing time derivatives


	// freezeout energy density
	precision e_switch = equilibriumEnergyDensity(T_switch / hbarc);


	print_parameters(lattice, t0, T_switch, etabar_const);


	// allocate memory for computational grid points
	allocate_memory(ncx * ncy * ncz);


	// fluid dynamic initialization
	precision t = t0;
	set_initial_conditions(t, lattice, initial, hydro);		// initial conditions for (q, u)
	set_ghost_cells(q, e, u, nx, ny, nz);					// initialize ghost cells in q

	printf("\n");

	int s = central_index(nx, ny, nz, ncx, ncy, ncz);

	double steps = 0;
	clock_t start = clock();

	// fluid dynamic evolution
	for(int n = 0; n <= nt; n++)
	{
		if(n % FREQ == 0)	// output variables to file for testing
		{
			precision e_s = e[s] * hbarc;

			print_hydro_center(t, e_s, s);

			//output_dynamical_variables(t, nx, ny, nz, dt, dx, dy, dz, initial, etabar_const);

			if(e[s] < e_switch) 	// replace with freezeout finder not finding any cells
			{
				printf("\nReached freezeout temperature at the center.\n\n");
				break;	// need to change it so that all cells below freezeout temperature
			}
		}

		evolve_hydro_one_time_step(t, dt, dt_prev, nx, ny, nz, dx, dy, dz, etabar_const);

		t += dt;

		dt_prev = dt;

		if(adaptive_time_step)
		{
			// compute the smallest time scale in the fluid
			hydro_time_scales dt_hydro = compute_hydro_time_scales(t, q, e, u, up, nx, ny, nz, dt, dt_prev, dx, dy, dz, etabar_const, dt_min);
			dt = set_time_step(dt_hydro, dt_min);
		}

		steps += 1.0;
	}

	double duration     = (clock() - start) / (double)CLOCKS_PER_SEC;
	double spatial_grid = nx * ny * nz;

	cout << "Total time             = " << setprecision(4) << duration << " s\n";
	cout << "Average time/step      = " << setprecision(3) << 1000. * duration / steps << " ms\n";
	cout << "Average time/cell/step = " << setprecision(3) << 1000. * duration / (spatial_grid * steps) << " ms\n";

	free_memory();
}



