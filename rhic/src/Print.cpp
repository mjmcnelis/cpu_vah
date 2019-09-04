#include <stdlib.h>
#include <stdio.h>
#include <string>
#include "../include/Macros.h"
#include "../include/Hydrodynamics.h"
#include "../include/FileIO.h"
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"
#include "../include/Parameters.h"
using namespace std;


void print_line()
{
	printf("--------------------------------------------------------------------------------");
	printf("--------------------------------------------------------------------------------\n");
}


void print_hydro_mode(hydro_parameters hydro)
{
	string mode = "Running";
	if(hydro.run_hydro == 1) mode = "Testing";
	
#ifdef ANISO_HYDRO
	printf("\n:::::::::::::::::::::::::::::::::::::::::::\n");
	printf(":::  %s viscous anisotropic hydro  :::\n", mode.c_str());
	printf(":::::::::::::::::::::::::::::::::::::::::::\n\n");
#else	
	printf("\n:::::::::::::::::::::::::::::::::::::::::::\n");
	printf(":::   %s second order viscous hydro   :::\n", mode.c_str());
	printf(":::::::::::::::::::::::::::::::::::::::::::\n\n");	
#endif
}


void print_run_time(double duration, double steps, lattice_parameters lattice)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	printf("Total time               = %.3g s\n", duration);
	printf("Number of time steps     = %d\n", (int)steps);
	printf("Average time/step        = %.3g ms\n", 1000. * duration / steps);
	printf("Average time/cell/step   = %.3g ms\n", 1000. * duration / (nx * ny * nz * steps));
}


void print_hydro_center(int n, double t, lattice_parameters lattice, hydro_parameters hydro)
{	
	int s = central_index(lattice);

	if(n == 0)
	{
		printf("\tn\t|\tt (fm/c)\t|\tT (GeV)\t\t|\te (GeV/fm^3)\t|\tpeq (GeV/fm^3)\t|\tpl (GeV/fm^3)\t|\tpt (GeV/fm^3)\n");
		print_line();
	}
	precision e_s = e[s] * hbarc;
	precision p = equilibriumPressure(e[s]) * hbarc;
	precision T = effectiveTemperature(e[s], hydro.conformal_eos_prefactor) * hbarc;
#ifdef ANISO_HYDRO
	precision pl  = q[s].pl * hbarc;
	#if (PT_MATCHING == 1)
		precision pt = q[s].pt * hbarc;
	#else
		precision pt = (e_s - pl) / 2.;
	#endif
	printf("\t%d\t|\t%.3f\t\t|\t%.3f\t\t|\t%.3f\t\t|\t%.3f\t\t|\t%.3f\t\t|\t%.3f\n", n, t, T, e_s, p, pl, pt);
#else
	printf("%d\tt = %.3f fm/c\t\te = %.3f GeV/fm^3\tpeq = %.3f GeV/fm^3\tT = %.3f GeV\n", n, t, e_s, p, T);
#endif
}


void print_parameters(lattice_parameters lattice, hydro_parameters hydro)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision dx = lattice.lattice_spacing_x;
	precision dy = lattice.lattice_spacing_y;
	precision dz = lattice.lattice_spacing_eta;

	precision dt = lattice.fixed_time_step;
	string time_step = "(fixed)";

	if(lattice.adaptive_time_step) 
	{
		dt = lattice.min_time_step;
		time_step = "(adaptive)";
	}
	
	// lattice parameters
	printf("Time resolution     = %.3g fm/c %s\n", dt, time_step.c_str());
	printf("Spatial grid points = %d x %d x %d\n", nx, ny, nz);
	printf("Spatial resolution  = [%.3f fm, %.3f fm, %.3f]\n", dx, dy, dz);
	printf("Spatial dimensions  = %.3f fm  x  %.3f fm  x  %.3f\n\n", (nx - 1.) * dx, (ny - 1.) * dy, (nz - 1.) * dz);

	// hydro parameters
	printf("Initial time 	       = %.3f fm/c\n", 	hydro.tau_initial);

	if(hydro.temperature_etas)
	{
		printf("Shear viscosity        = Temperature dependent\n");
	}
	else
	{
		printf("Shear viscosity        = %.3f (fixed)\n", hydro.constant_etas);
	}

	printf("Freezeout temperature  = %.3f GeV\n",	hydro.freezeout_temperature_GeV);
	printf("Flux limiter           = %.2f\n",		hydro.flux_limiter);
	printf("Minimum energy density = %.2e\n",		hydro.energy_min);

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

