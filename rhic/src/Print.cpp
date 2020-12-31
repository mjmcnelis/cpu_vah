#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <math.h>
#include "../include/Macros.h"
#include "../include/Hydrodynamics.h"
#include "../include/Output.h"
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"
#include "../include/Parameters.h"
#include "../include/Print.h"
#include "../include/OpenMP.h"
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
	string dimension = "3+1d";
#ifdef BOOST_INVARIANT
	dimension = "2+1d";
#endif

#ifdef ANISO_HYDRO
	printf("\n:::::::::::::::::::::::::::::::::::::::::::\n");
	printf(":::   Running %s anisotropic hydro    :::\n", dimension.c_str());
	printf(":::::::::::::::::::::::::::::::::::::::::::\n\n");
#else
	printf("\n:::::::::::::::::::::::::::::::::::::::::::::::::::\n");
	printf(":::   Running %s second order viscous hydro   :::\n", dimension.c_str());
	printf(":::::::::::::::::::::::::::::::::::::::::::::::::::\n\n");
#endif
}


void print_run_time(double t, double duration, double steps, lattice_parameters lattice, int sample)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

    int threads = 1;
#ifdef _OPENMP
    threads = omp_get_max_threads();
#endif

	printf("\nLifetime               = %.5g fm/c\n", t);
	printf("Run time               = %.4g s\n", duration);
	printf("Number of time steps   = %d\n", (int)steps);
	printf("Average time/step      = %.4g s\n", duration / steps);
	printf("Average time/cell/step = %.4g ms\n", 1000. * duration / (nx * ny * nz * steps));
    printf("Number of threads      = %d\n", threads);

#ifdef BENCHMARKS
	FILE * benchmarks;

	if(sample == 0)
	{
		benchmarks = fopen("output/benchmarks/benchmarks.dat", "a");
	}
	else
	{
		char fname[255];
		sprintf(fname, "output/benchmarks/benchmarks_%d.dat", sample);
		benchmarks = fopen(fname, "a");
	}

	fprintf(benchmarks, "%.5g\t%.5g\t%d\t%.5g\t%.5g\t%d\n", t, duration, (int)steps, duration / steps, 1000. * duration / (nx * ny * nz * steps), threads);
	fclose(benchmarks);
#endif
}

hydro_max compute_hydro_max(const precision * const __restrict__ e, const fluid_velocity * const __restrict__ u, precision t, lattice_parameters lattice, hydro_parameters hydro)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision ut_max = 0.;
	precision T_max = 0.;
	precision e_max = 0.;

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

				if(ut > ut_max)
				{
					ut_max = ut;
				}

				precision e_s = e[s];

				equation_of_state_new eos(e_s, hydro.conformal_eos_prefactor);
				precision T = eos.T;

				if(e_s > e_max)
				{
					e_max = e_s;
				}

				if(T > T_max)
				{
					T_max = T;
				}
			}
		}
	}

	hydro_max hydro_max_values;
	hydro_max_values.ut_max = ut_max;
	hydro_max_values.e_max = e_max;
	hydro_max_values.T_max = T_max;

	return hydro_max_values;
}


int compute_total_regulations(const int * const __restrict__ regulation, lattice_parameters lattice)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	int total_regulations = 0;

	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				total_regulations += regulation[s];
			}
		}
	}

	return total_regulations;
}




void print_hydro_center(int n, double t, lattice_parameters lattice, hydro_parameters hydro, long cells_above_Tswitch)
{
	int s = central_index(lattice);

	if(n == 0)
	{
	#ifdef ANISO_HYDRO

		// X_reg = percentage of anisotropic variable (X) regulations in the grid

		printf("\tn\t|\tt\t|\tT\t|\te\t|\tp\t|\tpl\t|\tpt\t|\temax\t|\tXreg\t|\te > esw\t|\n");
	#else
		printf("\tn\t|\tt\t|\tT\t|\te\t|\tp\t|\te_max\t|\tT_max\t|\tut_max\t|\tvisc_reg\t|\tT > Tsw\t|\n");
	#endif
		print_line();
	}
	precision e_s = e[s] * hbarc;

	equation_of_state_new eos(e[s], hydro.conformal_eos_prefactor);
	precision p = eos.equilibrium_pressure() * hbarc;
	precision T = eos.T * hbarc;

	hydro_max hydro_max_values = compute_hydro_max(e, u, t, lattice, hydro);
	precision gamma_max = hydro_max_values.ut_max;
	precision e_max = hbarc * hydro_max_values.e_max;
	precision T_max = hbarc * hydro_max_values.T_max;

	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision percent_above_Tsw = 100. * (double)cells_above_Tswitch / (double)(nx * ny * nz);
#ifdef ANISO_HYDRO
	precision pl = q[s].pl * hbarc;
	precision pt = q[s].pt * hbarc;

	// T = GeV
	// e, pl, pt = GeV/fm^3

	double percent_aniso_reg = 0;
#ifdef LATTICE_QCD
	percent_aniso_reg = 100. * (double)compute_total_regulations(aniso_regulation, lattice) / (double)(nx * ny * nz);
#endif

	printf("\t%d\t|\t%.4g\t|\t%.4g\t|\t%.4g\t|\t%.4g\t|\t%.4g\t|\t%.4g\t|\t%.4g\t|\t%.1f%%\t|\t%.1f%%\t|\n", n, t, T, e_s, p, pl, pt, e_max, percent_aniso_reg, percent_above_Tsw);
#else

#ifdef MONITOR_REGULATIONS
	double percent_viscous_reg = 100. * (double)compute_total_regulations(viscous_regulation, lattice) / (double)(nx * ny * nz);
#else
	double percent_viscous_reg = sqrt(-1.);		// print nan
#endif

	printf("\t%d\t|\t%.4g\t|\t%.4g\t|\t%.4g\t|\t%.4g\t|\t%.4g\t|\t%.3g\t|\t%.3g\t|\t%.1f%%\t|\t%.1f%%\t|\n", n, t, T, e_s, p, e_max, T_max, gamma_max, percent_viscous_reg, percent_above_Tsw);
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
	printf("Initial time  = %.3g fm/c\n", hydro.tau_initial);
	printf("Initial pl/pt = %.3g\n\n", 	hydro.plpt_ratio_initial);

	if(hydro.temperature_etas)
	{
		printf("etas(T)  = %.3f  +  (T - %.3f GeV) . (%.2f GeV^-1 . Theta(%.3f GeV - T)  +  %.2f GeV^-1 . Theta(T - %.3f GeV))\n\n", hydro.etas_etask, hydro.etas_Tk_GeV, hydro.etas_aL, hydro.etas_Tk_GeV, hydro.etas_aH, hydro.etas_Tk_GeV);
	}
	else
	{
		printf("eta/s(T)  = %.3f (fixed)\n\n", hydro.constant_etas);
	}

#ifdef CONFORMAL_EOS
	printf("zetas = 0\n\n");
#else
	char sign = '+';
	if(hydro.zetas_skew < 0)
	{
		sign = '-';
	}
	printf("zetas(T) = %.3f Lambda^2 / (Lambda^2  +  (T - %.3f GeV)^2)\t\tLambda = %.3f GeV . (1 %c %.3f sign(T - %.3f GeV))\n\n", hydro.zetas_normalization_factor, hydro.zetas_peak_temperature_GeV, hydro.zetas_width_GeV, sign, fabs(hydro.zetas_skew), hydro.zetas_peak_temperature_GeV);
#endif

	if(hydro.kinetic_theory_model == 0)
	{
		printf("Kinetic theory model = Small mass expansion\n\n");
	}
	else
	{
		printf("Kinetic theory model = Quasiparticle\n\n");
	}

	precision T_switch = hydro.freezeout_temperature_GeV;
	precision e_switch = hbarc * equilibrium_energy_density_new(T_switch / hbarc, hydro.conformal_eos_prefactor);

	printf("Freezeout temperature     = %.3g MeV\n", 		T_switch * 1000.);
	printf("Freezeout energy density  = %.3g GeV/fm^3\n",	e_switch);
	printf("Flux limiter              = %.2f\n",			hydro.flux_limiter);
	printf("Minimum energy density    = %.2e fm^-4\n",		hydro.energy_min);
	printf("Minimum pressure          = %.2e fm^-4\n",		hydro.pressure_min);

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

#ifdef VORTICITY
	printf("\nVorticity terms = On\n");
#else
	printf("\nVorticity terms = Off\n");
#endif
}
