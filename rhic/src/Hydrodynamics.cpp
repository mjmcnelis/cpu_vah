/*
 * HydroPlugin.c
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <ctime>
#include "../include/Hydrodynamics.h"
#include "../include/DynamicalVariables.h"
#include "../include/Parameters.h"
#include "../include/FileIO.h"
#include "../include/InitialConditions.h"
#include "../include/KurganovTadmor.h"
#include "../include/EquationOfState.h"

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


void run_hydro(void * latticeParams, void * initCondParams, void * hydroParams)
{
	const double hbarc = 0.197326938;

	// Parameters
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
	struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;

	// System configuration
	int nx = lattice->numLatticePointsX;									// physical grid
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;
	int nt = lattice->numProperTimePoints;

	int ncx = lattice->numComputationalLatticePointsX;						// (physical + ghost + white) grid
	int ncy = lattice->numComputationalLatticePointsY;
	int ncz = lattice->numComputationalLatticePointsRapidity;

	PRECISION dx = lattice->latticeSpacingX;								// lattice spacing
	PRECISION dy = lattice->latticeSpacingY;
	PRECISION dz = lattice->latticeSpacingRapidity;
	double dt = lattice->latticeSpacingProperTime;							// time resolution

	PRECISION etabar = (PRECISION)(hydro->shearViscosityToEntropyDensity);	// shear viscosity


	double t0 = hydro->initialProperTimePoint;								// initial longitudinal proper time
																			// if use F.S. need t_fs instead


	const double T_switch = (hydro->freezeoutTemperatureGeV) / hbarc;		// switching temperature [fm^-1]
	const double e_switch = equilibriumEnergyDensity(T_switch);				// switching energy density [fm^-4]


	printf("Time resolution    = %.3f fm/c\n", dt);
	printf("Hydro grid         = %d x %d x %d\n", nx, ny, nz);
	printf("Spatial resolution = (%.3f fm, %.3f fm, %.3f)\n", dx, dy, dz);
	printf("Spatial size       = %.3f fm  x  %.3f fm  x  %.3f\n", 0.5 * (nx - 1.0) * dx, 0.5 * (ny - 1.0) * dy, 0.5 * (nz - 1.0) * dz);

	printf("\nFreezeout temperature = %.3f GeV\t(e = %.3f GeV/fm^3)\n", T_switch * hbarc, e_switch * hbarc);
	printf("Shear viscosity       = %.3f\n", etabar);
#ifdef CONFORMAL_EOS
	printf("Equation of state     = Conformal\n");
#else
	printf("Equation of state     = QCD\n");
#endif

	// allocate memory for computational grid points
	allocate_memory(ncx * ncy * ncz);

	// fluid dynamic initialization
	double t = t0;
	set_initial_conditions(t, latticeParams, initCondParams, hydroParams);	// generate initial conditions (Tmunu, e, p, u, up, pl, pimunu, Wmu)
	set_ghost_cells(q, e, p, u, nx, ny, nz, ncx, ncy);						// initialize ghost cells (all current variables)

	printf("\n");

	int sctr = central_index(nx, ny, nz, ncx, ncy, ncz);	// central index

	//cout << sctr << endl;

	double steps = 0.0;
	clock_t start = clock();

	// fluid dynamic evolution
	for(int n = 0; n <= nt; n++)
	{
		// output variables to file occasionally (for testing)
		if(n % FREQ == 0)
		{
			printf("t = %.3f fm/c\t\te = %.3f GeV/fm^3\tp = %.3f GeV/fm^3\tpl = %.3f GeV/fm^3\tT = %.3f GeV\n", t, e[sctr] * hbarc, p[sctr] * hbarc, q->pl[sctr] * hbarc, effectiveTemperature(e[sctr]) * hbarc);

			output_dynamical_variables(t, nx, ny, nz, dx, dy, dz);

			if(e[sctr] < e_switch) 	// replace with freezeout finder not finding any cells
			{
				printf("\nReached freezeout temperature at the center.\n\n");
				break;	// need to change it so that all cells below freezeout temperature
			}
		}
		rungeKutta2(t, dt, q, Q, nx, ny, nz, ncx, ncy, dx, dy, dz, etabar);		// evolve one time step
		steps += 1.0;
		t += dt;
	}

	// time benchmarks
	double duration = (clock() - start) / (double)CLOCKS_PER_SEC;
	double spatial_grid = (double)(nx * ny * nz);

	cout << "Total time             = " << setprecision(4) << duration << " s\n";
	cout << "Average time/cell/step = " << setprecision(3) << 1000.0 * duration / (spatial_grid * steps) << " ms\n";

	free_memory();
}



