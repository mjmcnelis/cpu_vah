#include <stdlib.h>
#include <stdio.h>
#include "../include/Hydrodynamics.h"
#include "../include/FluxTerms.h"
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"
#include "../include/Parameters.h"


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


void print_parameters(void * lattice, double t0, double T_switch, double etabar)
{
	struct LatticeParameters * lattice_parameters = (struct LatticeParameters *) lattice;

	int nx = lattice_parameters->lattice_points_x;
	int ny = lattice_parameters->lattice_points_y;
	int nz = lattice_parameters->lattice_points_eta;

	precision dx = lattice_parameters->lattice_spacing_x;
	precision dy = lattice_parameters->lattice_spacing_y;
	precision dz = lattice_parameters->lattice_spacing_eta;

	int adaptive_time_step = lattice_parameters->adaptive_time_step;

	precision dt = lattice_parameters->fixed_time_step;
	if(adaptive_time_step) 
	{
		dt = lattice_parameters->min_time_step;
	}

// lattice parameters
	if(adaptive_time_step)
	{
		printf("Time resolution     = %.3g (adaptive) \n", dt);
	}
	else
	{
		printf("Time resolution     = %.3g fm/c (fixed)\n", dt);
	}
	printf("Spatial grid points = %d x %d x %d\n", nx, ny, nz);
	printf("Spatial resolution  = [%.3f fm, %.3f fm, %.3f]\n", dx, dy, dz);
	printf("Spatial dimensions  = %.3f fm  x  %.3f fm  x  %.3f\n", (nx - 1.) * dx, (ny - 1.) * dy, (nz - 1.) * dz);
// hydro parameters
	printf("\nInitial time 	       = %.3f fm/c\n", t0);
	printf("Shear viscosity        = %.3f\n", etabar);
	printf("Freezeout temperature  = %.3f GeV\n", T_switch);
	printf("Flux limiter           = %.2f\n", THETA);
	printf("Minimum energy density = %.2e\n", E_MIN);
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


