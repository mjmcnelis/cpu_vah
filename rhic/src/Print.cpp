#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <math.h>
#include "../include/Macros.h"
#include "../include/Hydrodynamics.h"
#include "../include/FileIO.h"
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"
#include "../include/Parameters.h"
using namespace std;

inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}

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
	printf("\n::::::::::::::::::::::::::::::::::::::::::::::\n");
	printf(":::   %s second order viscous hydro   :::\n", mode.c_str());
	printf("::::::::::::::::::::::::::::::::::::::::::::::\n\n");
#endif
}


void print_run_time(double duration, double steps, lattice_parameters lattice)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	printf("Total time               = %.4g s\n", duration);
	printf("Number of time steps     = %d\n", (int)steps);
	printf("Average time/step        = %.4g s\n", duration / steps);
	printf("Average time/cell/step   = %.4g ms\n", 1000. * duration / (nx * ny * nz * steps));
}

precision compute_max_ut(const fluid_velocity * const __restrict__ u, precision t, lattice_parameters lattice)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision max = 0.0;

	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				precision ux = u[s].ux;
				precision uy = u[s].uy;
			#ifndef BOOST_INVARIANT
				precision un = u[s].un;
			#else
				precision un = 0;
			#endif
				precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t * t * un * un);

				if(ut > max)
				{
					max = ut;
				}
			}
		}
	}
	return max;
}


void print_hydro_center(int n, double t, lattice_parameters lattice, hydro_parameters hydro, long cells_above_Tswitch)
{
	int s = central_index(lattice);

	if(n == 0)
	{
	#ifdef ANISO_HYDRO
		printf("\tn\t|\tt\t|\tT\t|\te\t|\tp\t|\tpl\t|\tpt\t|\tut_max\t|\tT > Tsw\t|\n");
	#else
		printf("\tn\t|\tt\t|\tT\t|\te\t|\tp\t|\tut_max\t|\tT > Tsw\t|\n");
	#endif
		print_line();
	}
	precision e_s = e[s] * hbarc;

	equation_of_state eos(e[s]);
	precision p = eos.equilibrium_pressure() * hbarc;
	precision T = eos.effective_temperature(hydro.conformal_eos_prefactor) * hbarc;

	precision gamma_max = compute_max_ut(u, t, lattice);

	long nx = lattice.lattice_points_x;
	long ny = lattice.lattice_points_y;
	long nz = lattice.lattice_points_eta;

	precision percent_above_Tsw = 100. * (double)cells_above_Tswitch / (nx * ny * nz);
#ifdef ANISO_HYDRO
	precision pl = q[s].pl * hbarc;
	precision pt = q[s].pt * hbarc;

	// T = GeV
	// e, pl, pt = GeV/fm^3

	printf("\t%d\t|\t%.4g\t|\t%.4g\t|\t%.5g\t|\t%.5g\t|\t%.5g\t|\t%.5g\t|\t%.5g\t|\t%.3f%%\t|\n", n, t, T, e_s, p, pl, pt, gamma_max, percent_above_Tsw);
#else
	printf("\t%d\t|\t%.4g\t|\t%.4g\t|\t%.5g\t|\t%.5g\t|\t%.5g\t|\t%.3f%%\t|\n", n, t, T, e_s, p, gamma_max, percent_above_Tsw);
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
	printf("Initial time 	    = %.3g fm/c\n", hydro.tau_initial);
	printf("Initial pl/pt ratio = %.3g\n\n", 	hydro.plpt_ratio_initial);

	if(hydro.temperature_etas)
	{
		printf("Shear viscosity model:\teta/s = max(%.3f, %.3f  +  %.3f(T / GeV - 0.154))\n\n", hydro.etas_min, hydro.etas_min, hydro.etas_slope);
	}
	else
	{
		printf("Shear viscosity model:\teta/s = %.3f (fixed)\n\n", hydro.constant_etas);
	}

	printf("Bulk viscosity normalization    = %.3g\n", 			hydro.zetas_normalization_factor);
	printf("Bulk viscosity peak temperature = %.3g MeV\n\n", 	hydro.zetas_peak_temperature_GeV * 1000.);

	if(hydro.kinetic_theory_model == 0)
	{
		printf("Kinetic theory model = Small mass expansion\n\n");
	}
	else
	{
		printf("Kinetic theory model = Quasiparticle\n\n");
	}

	printf("Freezeout temperature  = %.3g MeV\n",	hydro.freezeout_temperature_GeV * 1000.);
	printf("Flux limiter           = %.2f\n",		hydro.flux_limiter);
	printf("Minimum energy density = %.2e\n",		hydro.energy_min);
	printf("Minimum pressure       = %.2e\n",		hydro.pressure_min);

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
#ifdef B_FIELD
	printf("Mean field                      = On\n");
#else
	printf("Mean field                      = Off\n");
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


