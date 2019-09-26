#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../include/Parameters.h"
#include "../include/Precision.h"
#include "../include/Macros.h"
#include "../include/EquationOfState.h"
#include "../include/DynamicalVariables.h"
#include "../include/TransportCoefficients.h"
#include "../include/Projections.h"
#include "../include/NeighborCells.h"
#include "../include/Hydrodynamics.h"
#include "../include/AnisoBjorken.h"
#include "../include/AnisoGubser.h"

inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


int central_index(lattice_parameters lattice)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	int ncx = nx + 4;
	int ncy = ny + 4;
	int ncz = nz + 4;

	// (not sure what this means)
	int i = (nx % 2 == 0) ? ncx/2 : (ncx-1)/2;
	int j = (ny % 2 == 0) ? ncy/2 : (ncy-1)/2;
	int k = (nz % 2 == 0) ? ncz/2 : (ncz-1)/2;

	return linear_column_index(i, j, k, ncx, ncy);
}


double compute_conformal_aL(double pl, double e)
{
	precision x   = pl / e;		// x = pl / e
	precision x2  = x   * x;
	precision x3  = x2  * x;
	precision x4  = x3  * x;
	precision x5  = x4  * x;
	precision x6  = x5  * x;
	precision x7  = x6  * x;
	precision x8  = x7  * x;
	precision x9  = x8  * x;
	precision x10 = x9  * x;
	precision x11 = x10 * x;
	precision x12 = x11 * x;
	precision x13 = x12 * x;
	precision x14 = x13 * x;

	precision aL = (5.6098342562962155e-24 + 1.0056714201158781e-17*x + 8.574287549260127e-13*x2 + 8.639689853874967e-9*x3 + 0.000014337184308704522*x4 +
     0.0047402683487226555*x5 + 0.3461801244895056*x6 + 5.3061287395562*x7 + 3.7804213528647956*x8 - 55.646719325650224*x9 +
     71.68906037132133*x10 + 0.6485422288016947*x11 - 52.86438720903515*x12 + 32.635674688615836*x13 - 5.899614102635062*x14)/
   (1.2460117685059638e-20 + 3.9506205613753145e-15*x + 1.090135069930889e-10*x2 + 4.2931027828550746e-7*x3 + 0.00030704101799886117*x4 +
     0.04575504592470687*x5 + 1.4479634250149949*x6 + 6.077429142899631*x7 - 29.171395065126873*x8 + 13.501854646832847*x9 +
     65.98203155631907*x10 - 111.65365949648432*x11 + 71.83676912638525*x12 - 19.66184593458614*x13 + 1.5947903161928916*x14);

   	return fmax(0.001, fmin(aL, 20.0));
}


inline precision central_derivative(const precision * const __restrict__ f, int n, precision dx)
{
	// f[n] = fm  |	 f[n+1] = fp  (appears counterintuitive it's f's array structure)
	return (f[n + 1] - f[n]) / (2. * dx);
}


void output_residual_shear_validity(const hydro_variables * const __restrict__ q, const fluid_velocity * const __restrict__ u, const precision * const e, double t, lattice_parameters lattice)
{
#ifdef PIMUNU
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision dx = lattice.lattice_spacing_x;
	precision dy = lattice.lattice_spacing_y;
	precision dn = lattice.lattice_spacing_eta;

	FILE *RpiT_inverse;
	char fname1[255];
	sprintf(fname1, "output/RpiT_inv_%.3f.dat", t);
	RpiT_inverse = fopen(fname1, "w");

	precision t2 = t * t;
	precision t4 = t2 * t2;

	for(int k = 2; k < nz + 2; k++)
	{
		double z = (k - 2. - (nz - 1.)/2.) * dn;

		for(int j = 2; j < ny + 2; j++)
		{
			double y = (j - 2. - (ny - 1.)/2.) * dy;

			for(int i = 2; i < nx + 2; i++)
			{
				double x = (i - 2. - (nx - 1.)/2.) * dx;

				int s = linear_column_index(i, j, k, nx + 4, ny + 4);
				
				precision e_s = e[s];
				precision pt;

				precision pitt = q[s].pitt;			// get shear stress
				precision pitx = q[s].pitx;
				precision pity = q[s].pity;
				precision pitn = 0;
				precision pixx = q[s].pixx;
				precision pixy = q[s].pixy;
				precision pixn = 0;
				precision piyy = q[s].piyy;
				precision piyn = 0;
				precision pinn = 0;

			#ifndef BOOST_INVARIANT
				pitn = q[s].pitn;
				pixn = q[s].pixn;
				piyn = q[s].piyn;
				pinn = q[s].pinn;
			#else
			#ifndef ANISO_HYDRO
				pinn = q[s].pinn;
			#endif
			#endif

			#ifdef ANISO_HYDRO						// get transverse pressure

			#if (PT_MATCHING == 1)					
				pt = q[s].pt;
			#else
				pt = (e_s - q[s].pl) / 2.;
			#endif

			#else									// get transverse pressure and double-transverse project viscous hydro shear stress
				precision P = equilibriumPressure(e_s);
				
			#ifdef PI
				P += q[s].Pi;
			#endif

				precision ux = u[s].ux;
				precision uy = u[s].uy;
				precision un = 0;

			#ifndef BOOST_INVARIANT
				un = u[s].un;
			#endif
				precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t2 * un * un);
				precision utperp = sqrt(1.  +  ux * ux  +  uy * uy);

				precision zt = t * un / utperp;
				precision zn = ut / t / utperp;

				transverse_projection Xi(ut, ux, uy, un, zt, zn, t2);	
				double_transverse_projection Xi_2(Xi, t2, t4);

				Xi_2.double_transverse_project_tensor(pitt, pitx, pity, pitn, pixx, pixy, pixn, piyy, piyn, pinn);

				pt = P  -  (zt * zt * pitt  +  t4 * zn * zn * pinn  +  2. * t2 * zt * zn * pitn) / 2.;
			#endif

				precision pi_mag = sqrt(fabs(pitt * pitt  +  pixx * pixx  +  piyy * piyy  +  t4 * pinn * pinn  -  2. * (pitx * pitx  +  pity * pity  -  pixy * pixy  +  t2 * (pitn * pitn  -  pixn * pixn  -  piyn * piyn))));

				fprintf(RpiT_inverse, "%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, pi_mag / (sqrt(2.) * pt));
			}
		}
	}
	fclose(RpiT_inverse);
#endif
}


void output_viscous_bjorken(const hydro_variables * const __restrict__ q, const precision * const e, precision e0, precision t, lattice_parameters lattice, hydro_parameters hydro)
{
#ifndef ANISO_HYDRO
	FILE *energy, *plptratio;

	energy = fopen("output/eratio.dat", "a");
	plptratio = fopen("output/plptratio.dat", "a");

	int s = central_index(lattice);

	precision e_s = e[s];
	precision p = equilibriumPressure(e_s);

#ifdef PIMUNU
	precision pinn = q[s].pinn;
#else
	precision pinn = 0;
#endif

#ifdef PI
	precision Pi = q[s].Pi;
#else
	precision Pi = 0;
#endif

	precision pl = p  +  Pi  +  t * t * pinn;
	precision pt = p  +  Pi  -  t * t * pinn / 2.;

	fprintf(energy, 	"%.8f\t%.8e\n", t, e_s / e0);
	fprintf(plptratio, 	"%.8f\t%.8e\n", t, pl / pt);

	fclose(energy);
	fclose(plptratio);
#endif
}


void output_aniso_bjorken(const hydro_variables * const __restrict__ q, const precision * const e, double e0, double t, lattice_parameters lattice, hydro_parameters hydro)
{
#ifdef ANISO_HYDRO
	FILE *energy, *plptratio;

	energy = fopen("output/eratio.dat", "a");
	plptratio = fopen("output/plptratio.dat", "a");

	int s = central_index(lattice);

	precision e_s = e[s];
	precision pl = q[s].pl;
#if (PT_MATCHING == 1)
	precision pt = q[s].pt;
#else
	precision pt = (e_s - pl) / 2.;
#endif

	fprintf(energy, 	"%.8f\t%.8e\n", t, e_s / e0);
	fprintf(plptratio, 	"%.8f\t%.8e\n", t, pl / pt);

	fclose(energy);
	fclose(plptratio);
#endif
}


void output_gubser(const hydro_variables * const __restrict__ q, const fluid_velocity  * const __restrict__ u, const precision * const e, double t, lattice_parameters lattice)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision dx = lattice.lattice_spacing_x;
	precision dy = lattice.lattice_spacing_y;
	precision dz = lattice.lattice_spacing_eta;

	FILE *energy, *plptratio, *uxplot, *urplot;
	char fname1[255], fname2[255], fname3[255], fname4[255];

	sprintf(fname1, "output/e_%.3f.dat", t);
	sprintf(fname2, "output/plpt_%.3f.dat", t);
	sprintf(fname3, "output/ux_%.3f.dat", t);
	sprintf(fname4, "output/ur_%.3f.dat", t);

	energy      = fopen(fname1, "w");
	plptratio 	= fopen(fname2, "w");
	uxplot    	= fopen(fname3, "w");
	urplot		= fopen(fname4, "w");

	precision t2 = t * t;
	precision t4 = t2 * t2;

	for(int k = 2; k < nz + 2; k++)
	{
		double z = (k - 2. - (nz - 1.)/2.) * dz;

		for(int j = 2; j < ny + 2; j++)
		{
			double y = (j - 2. - (ny - 1.)/2.) * dy;

			for(int i = 2; i < nx + 2; i++)
			{
				double x = (i - 2. - (nx - 1.)/2.) * dx;

				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				precision ux = u[s].ux;
				precision uy = u[s].uy;
				precision ur = sqrt(ux * ux  +  uy * uy);
				precision e_s = e[s];

			#ifdef ANISO_HYDRO
				precision pl = q[s].pl;
			#else

			#ifndef BOOST_INVARIANT
				precision un = u[s].un;
				precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t2 * un * un);
				precision utperp = sqrt(1.  +  ux * ux  +  uy * uy);
				precision zt = t * un / utperp;
				precision zn = ut / t / utperp;
			#else
				precision zt = 0;
				precision zn = 1. / t;
			#endif

			#ifdef PIMUNU
				precision pitt = q[s].pitt;
				precision pinn = q[s].pinn;

			#ifndef BOOST_INVARIANT
				precision pitn = q[s].pitn;
			#else
				precision pitn = 0;
			#endif

			#else
				precision pitt = 0;
				precision pitn = 0;
				precision pinn = 0;
			#endif

				// pl = zmu.z.nu.Tmunu
				precision pl = e_s/3.  +  zt * zt * pitt  +  t4 * zn * zn * pinn  +  2. * t2 * zt * zn * pitn;
			#endif

				precision pt = (e_s - pl) / 2.;

				fprintf(energy, 	"%.5f\t%.3f\t%.3f\t%.8f\n", x, y, z, e_s * hbarc);	// GeV / fm^3
				fprintf(plptratio,	"%.5f\t%.3f\t%.3f\t%.8f\n", x, y, z, pl / pt);
				fprintf(uxplot, 	"%.5f\t%.3f\t%.3f\t%.8f\n", x, y, z, ux);
				fprintf(urplot, 	"%.5f\t%.3f\t%.3f\t%.8f\n", x, y, z, ur);
			}
		}
	}
	fclose(energy);
	fclose(plptratio);
	fclose(uxplot);
	fclose(urplot);
}


void output_optical_glauber(const hydro_variables * const __restrict__ q, const fluid_velocity  * const __restrict__ u, const precision * const e, double t, lattice_parameters lattice)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision dx = lattice.lattice_spacing_x;
	precision dy = lattice.lattice_spacing_y;
	precision dz = lattice.lattice_spacing_eta;

	FILE *energy, *plptratio, *uxplot, *uyplot;
	char fname1[255], fname2[255], fname3[255], fname4[255];

	sprintf(fname1, "output/e_%.3f.dat", t);
	sprintf(fname2, "output/plpt_%.3f.dat", t);
	sprintf(fname3, "output/ux_%.3f.dat", t);
	sprintf(fname4, "output/uy_%.3f.dat", t);

	energy      = fopen(fname1, "w");
	plptratio 	= fopen(fname2, "w");
	uxplot    	= fopen(fname3, "w");
	uyplot    	= fopen(fname4, "w");

	precision t2 = t * t;
	precision t4 = t2 * t2;

	for(int k = 2; k < nz + 2; k++)
	{
		double z = (k - 2. - (nz - 1.)/2.) * dz;

		for(int j = 2; j < ny + 2; j++)
		{
			double y = (j - 2. - (ny - 1.)/2.) * dy;

			for(int i = 2; i < nx + 2; i++)
			{
				double x = (i - 2. - (nx - 1.)/2.) * dx;

				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				precision ux = u[s].ux;
				precision uy = u[s].uy;
				precision e_s = e[s];

			#ifdef ANISO_HYDRO
				precision pl = q[s].pl;
			#if (PT_MATCHING == 1)
				precision pt = q[s].pt;
			#else
				precision pt = (e_s - pl) / 2.;
			#endif
			#else
				precision p = equilibriumPressure(e_s);

			#ifndef BOOST_INVARIANT
				precision un = u[s].un;
				precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t2 * un * un);
				precision utperp = sqrt(1.  +  ux * ux  +  uy * uy);
				precision zt = t * un / utperp;
				precision zn = ut / t / utperp;
			#else
				precision zt = 0;
				precision zn = 1. / t;
			#endif

			#ifdef PIMUNU
				precision pitt = q[s].pitt;
			#ifndef BOOST_INVARIANT
				precision pitn = q[s].pitn;
			#else
				precision pitn = 0;
			#endif
				precision pinn = q[s].pinn;
			#else
				precision pitt = 0, pitn = 0, pinn = 0;
			#endif
			#ifdef PI
				precision Pi = q[s].Pi;
			#else
				precision Pi = 0;
			#endif
				// pl = zmu.z.nu.Tmunu
				precision pl = p  +  Pi  +  zt * zt * pitt  +  t4 * zn * zn * pinn  +  2. * t2 * zt * zn * pitn;
				// pt = -Ximunu.Tmunu / 2 
				precision pt = p  +  Pi  -  (zt * zt * pitt  +  t4 * zn * zn * pinn  +  2. * t2 * zt * zn * pitn) / 2.;
			#endif

				fprintf(energy, 	"%.3f\t%.3f\t%.3f\t%.8e\n", x, y, z, e_s * hbarc);
				fprintf(plptratio,	"%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, pl / pt);
				fprintf(uxplot, 	"%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, ux);
				fprintf(uyplot, 	"%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, uy);
			}
		}
	}
	fclose(energy);
	fclose(plptratio);
	fclose(uxplot);
	fclose(uyplot);
}


void output_dynamical_variables(double t, double dt, lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro)
{
	int initial_type = initial.initialConditionType;

	if(initial_type == 1)
	{
		precision T0 = initial.initialCentralTemperatureGeV;
		precision e0 = equilibriumEnergyDensity(T0 / hbarc, hydro.conformal_eos_prefactor);

	#ifdef ANISO_HYDRO
		output_aniso_bjorken(q, e, e0, t, lattice, hydro);
	#else
		output_viscous_bjorken(q, e, e0, t, lattice, hydro);
	#endif
	}
	else if(initial_type == 2 || initial_type == 3)
	{
		output_gubser(q, u, e, t, lattice);
		output_residual_shear_validity(q, u, e, t, lattice);
	}
	else if(initial_type == 4)
	{
		output_optical_glauber(q, u, e, t, lattice);
		output_residual_shear_validity(q, u, e, t, lattice);
	}
}


void output_semi_analytic_solution_if_any(lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro)
{
	switch(initial.initialConditionType)
	{
		case 1:		// Bjorken
		{
			printf("\nRunning semi-analytic anisotropic Bjorken solution...\n");
			run_semi_analytic_aniso_bjorken(lattice, initial, hydro);
			break;
		}
		case 3:		// anisotropic Gubser
		{
			printf("\nRunning semi-analytic anisotropic Gubser solution!...\n");
			double T0_hat = run_semi_analytic_aniso_gubser(lattice, initial, hydro);
			break;
		}
		default:
		{
			printf("\nNo semi-analytic solution to run...\n");
			break;
		}
	}
}













