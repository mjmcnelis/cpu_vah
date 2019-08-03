#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <ctime>
#include "../include/Hydrodynamics.h"
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"
#include "../include/FluxTerms.h"
#include "../include/GhostCells.h"
#include "../include/Parameters.h"
#include "../include/FileIO.h"
#include "../include/InitialConditions.h"
#include "../include/KurganovTadmor.h"
#include "../include/EquationOfState.h"
using namespace std;

#define FREQ 100

const double hbarc = 0.197326938;


int central_index(int nx, int ny, int nz, int ncx, int ncy, int ncz)
{
	// (not sure what this means)
	int ictr = (nx % 2 == 0) ? ncx/2 : (ncx-1)/2;
	int jctr = (ny % 2 == 0) ? ncy/2 : (ncy-1)/2;
	int kctr = (nz % 2 == 0) ? ncz/2 : (ncz-1)/2;

	return ictr  +  ncx * (jctr  +  ncy * kctr);
}

void print_hydro_center(double t, double e, int s)
{
	precision p = equilibriumPressure(e / hbarc) * hbarc;
	precision T = effectiveTemperature(e / hbarc) * hbarc;
#ifdef ANISO_HYDRO
	precision pl  = q[s].pl * hbarc;
	#if (PT_MATCHING == 1)
		precision pt = q[s].pt * hbarc;
	#else
		precision pt = (e - pl) / 2.;
	#endif
	printf("t = %.3f fm/c\t\te = %.3f GeV/fm^3\tpeq = %.3f GeV/fm^3\tpl = %.3f GeV/fm^3\tpt = %.3f GeV/fm^3\tT = %.3f GeV\n", t, e, p, pl, pt, T);
#else
	printf("t = %.3f fm/c\t\te = %.3f GeV/fm^3\tpeq = %.3f GeV/fm^3\tT = %.3f GeV\n", t, e, p, T);
#endif
}


void print_parameters(int nx, int ny, int nz, double dt, double dx, double dy, double dz, double t0, double T_switch, double etabar)
{
// lattice parameters
	printf("Time resolution     = %.4f fm/c\n", dt);
	printf("Spatial grid points = %d x %d x %d\n", nx, ny, nz);
	printf("Spatial resolution  = [%.3f fm, %.3f fm, %.3f]\n", dx, dy, dz);
	printf("Spatial dimensions  = %.3f fm  x  %.3f fm  x  %.3f\n", (nx - 1.) * dx, (ny - 1.) * dy, (nz - 1.) * dz);
// hydro parameters
	printf("\nHydro time 	      = %.3f fm/c\n", t0);
	printf("Shear viscosity       = %.3f\n", etabar);
	printf("Freezeout temperature = %.3f GeV\n", T_switch * hbarc);
	printf("Flux limiter          = %.2f\n", THETA);
// equation of state
#ifdef CONFORMAL_EOS
	printf("\nEquation of state = Conformal\n\n");
#else
	printf("\nEquation of state = QCD\n\n");
#endif
// viscous pressures
#ifdef ANISO_HYDRO
#ifdef PIMUNU
	printf("Transverse shear stress 	= On\n");
#else
	printf("Transverse shear stress 	= Off\n");
#endif
#ifdef WTZMU
	printf("Longitudinal momentum diffusion = On\n");
#else
	printf("Longitudinal momentum diffusion = Off\n");
#endif
#else
#ifdef PIMUNU
	printf("Shear stress  = On\n");
#else
	printf("Shear stress  = Off\n");
#endif
#ifdef PI
	printf("Bulk pressure = On\n");
#else
	printf("Bulk pressure = Off\n");
#endif
#endif
}


void run_hydro(void * latticeParams, void * initCondParams, void * hydroParams)
{
	// parameters
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
	struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;

	// system configuration
	int initialType = initCond->initialConditionType;	// initial condition type

	int nx = lattice->numLatticePointsX;				// physical grid
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;
	int nt = lattice->numProperTimePoints;

	int ncx = nx + 4;									// computational grid = physical + ghost + white
	int ncy = ny + 4;
	int ncz = nz + 4;

	precision dx = lattice->latticeSpacingX;			// lattice spacing
	precision dy = lattice->latticeSpacingY;
	precision dz = lattice->latticeSpacingRapidity;
	double dt = lattice->latticeSpacingProperTime;		// time resolution

	precision etabar = hydro->shear_viscosity;			// shear viscosity

	double T0 = initCond->initialCentralTemperatureGeV;
	double e0 = equilibriumEnergyDensity(T0 / hbarc);
	double t0 = hydro->tau_initial;						// initial longitudinal proper time (if use F.S. need t_fs instead)
	double T_switch = (hydro->freezeoutTemperatureGeV) / hbarc;		// switching temperature [fm^-1]
	double e_switch = equilibriumEnergyDensity(T_switch);			// switching energy density [fm^-4]

	print_parameters(nx, ny, nz, dt, dx, dy, dz, t0, T_switch, etabar);

	allocate_memory(ncx * ncy * ncz);					// allocate memory for computational grid points

	// fluid dynamic initialization
	double t = t0;
	set_initial_conditions(t, latticeParams, initCondParams, hydroParams);	// generate initial (Tmunu, e, p, u, up, pl, pimunu, Wmu)
	set_ghost_cells(q, e, u, nx, ny, nz);									// initialize ghost cells (all current variables)

	printf("\n");

	int s = central_index(nx, ny, nz, ncx, ncy, ncz);	// central index

	double steps = 0;
	clock_t start = clock();

	// fluid dynamic evolution
	for(int n = 0; n <= nt; n++)
	{
		if(n % FREQ == 0)	// output variables to file for testing
		{
			precision e_s = e[s] * hbarc;

			print_hydro_center(t, e_s, s);

			output_dynamical_variables(t, nx, ny, nz, dt, dx, dy, dz, initialType, e0, etabar);

			if(e[s] < e_switch) 	// replace with freezeout finder not finding any cells
			{
				printf("\nReached freezeout temperature at the center.\n\n");
				break;	// need to change it so that all cells below freezeout temperature
			}
		}

		evolve_hydro_one_time_step(t, dt, nx, ny, nz, dx, dy, dz, etabar);
		steps += 1.;
		t += dt;
	}

	double duration     = (clock() - start) / (double)CLOCKS_PER_SEC;
	double spatial_grid = (double)(nx * ny * nz);

	cout << "Total time             = " << setprecision(4) << duration << " s\n";
	cout << "Average time/step      = " << setprecision(3) << 1000. * duration / steps << " ms\n";
	cout << "Average time/cell/step = " << setprecision(3) << 1000. * duration / (spatial_grid * steps) << " ms\n";

	free_memory();
}



