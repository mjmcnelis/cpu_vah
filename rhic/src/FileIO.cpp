#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../include/Parameters.h"
#include "../include/Precision.h"
#include "../include/Macros.h"
#include "../include/EquationOfState.h"
#include "../include/DynamicalVariables.h"
#include "../include/Projections.h"
#include "../include/NeighborCells.h"
#include "../include/Hydrodynamics.h"
#include "../include/AnisoBjorken.h"
#include "../include/ViscousBjorken.h"
#include "../include/AnisoGubser.h"
#include "../include/ViscousGubser.h"
#include "../include/Viscosities.h"

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


void output_mean_field(const hydro_variables * const __restrict__ q, const fluid_velocity * const __restrict__ u, const precision * const e, double t, lattice_parameters lattice, hydro_parameters hydro)
{
#ifdef LATTICE_QCD
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision dx = lattice.lattice_spacing_x;
	precision dy = lattice.lattice_spacing_y;
	precision dn = lattice.lattice_spacing_eta;

	FILE *db_peq;			// from time evolution of relaxation equation
	FILE *dbasy_peq;		// asymptotic approximation (ignore residual shear contribution)
	FILE *db2_peq; 			// 2nd order approximation

	char fname1[255];
	char fname2[255];
	char fname3[255];

	sprintf(fname1, "output/db_peq_%.3f.dat", t);
	sprintf(fname2, "output/dbasy_peq_%.3f.dat", t);
	sprintf(fname3, "output/db2_peq_%.3f.dat", t);

#ifdef ANISO_HYDRO
	db_peq 		= fopen(fname1, "w");
	dbasy_peq 	= fopen(fname2, "w");
#else
	db2_peq 	= fopen(fname3, "w");
#endif

	precision t2 = t * t;
	precision t4 = t2 * t2;

	int stride_y = nx + 4;					// strides for neighbor cells along x, y, n (stride_x = 1)
	int stride_z = (nx + 4) * (ny + 4);		// stride formulas based from linear_column_index()

	precision ui1[6];		// fluid velocity of neighbor cells along x [i-1, i+1]
	precision uj1[6];		// fluid velocity of neighbor cells along y [j-1, j+1]
	precision uk1[6];		// fluid velocity of neighbor cells along n [k-1, k+1]

	precision vxi[4];		// vx of neighbor cells along x [i-2, i-1, i+1, i+2]
	precision vyj[4];		// vy of neighbor cells along y [j-2, j-1, j+1, j+2]
	precision vnk[4];		// vn of neighbor cells along n [k-2, k-1, k+1, k+2]


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

				equation_of_state_new eos(e_s, hydro.conformal_eos_prefactor);
				precision p = eos.equilibrium_pressure();
				precision cs2 = eos.speed_of_sound_squared();
				precision beq = eos.equilibrium_mean_field();
				precision T = eos.T;
				precision mass = T * eos.z_quasi();
				precision mdmde = eos.mdmde_quasi();

			#ifdef ANISO_HYDRO
				precision pl = q[s].pl;
				precision pt = q[s].pt;
				precision b = q[s].b;
				precision Pi = (pl + 2.*pt) / 3. - p;
			#else
				precision Pi = q[s].Pi;
			#endif

				precision ux = u[s].ux;
				precision uy = u[s].uy;

			#ifndef BOOST_INVARIANT
				precision un = u[s].un;
			#else
				precision un = 0;
			#endif

				precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t2 * un * un);
				precision utperp = sqrt(1.  +  ux * ux  +  uy * uy);

				precision zt = t * un / utperp;
				precision zn = ut / t / utperp;

				int simm = s - 2;			// neighbor cell indices (x)
				int sim  = s - 1;
				int sip  = s + 1;
				int sipp = s + 2;

				int sjmm = s - 2*stride_y;	// neighbor cell indices (y)
				int sjm  = s - stride_y;
				int sjp  = s + stride_y;
				int sjpp = s + 2*stride_y;

				int skmm = s - 2*stride_z;	// neighbor cell indices (n)
				int skm  = s - stride_z;
				int skp  = s + stride_z;
				int skpp = s + 2*stride_z;

				get_fluid_velocity_neighbor_cells(u[simm], u[sim], u[sip], u[sipp], u[sjmm], u[sjm], u[sjp], u[sjpp], u[skmm], u[skm], u[skp], u[skpp], ui1, uj1, uk1, vxi, vyj, vnk, t2);


			#ifdef ANISO_HYDRO
				fprintf(db_peq, "%.2f\t%.2f\t%.2f\t%.6f\n", x, y, z, (b - beq) / p);
			#else

			#endif
			}
		}
	}
#ifdef ANISO_HYDRO
	fclose(db_peq);
#else

#endif
#endif
}


void output_residual_shear_validity(const hydro_variables * const __restrict__ q, const fluid_velocity * const __restrict__ u, const precision * const e, double t, lattice_parameters lattice, hydro_parameters hydro)
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
	sprintf(fname1, "output/piperp_pt_%.3f.dat", t);
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
				precision pt = q[s].pt;
			#else 									// get transverse pressure and double-transverse project viscous hydro shear stress

				equation_of_state_new eos(e_s, hydro.conformal_eos_prefactor);
				precision P = eos.equilibrium_pressure();

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

				precision pt = P  -  (zt * zt * pitt  +  t4 * zn * zn * pinn  +  2. * t2 * zt * zn * pitn) / 2.;
			#endif

				precision pi_mag = sqrt(fabs(pitt * pitt  +  pixx * pixx  +  piyy * piyy  +  t4 * pinn * pinn  -  2. * (pitx * pitx  +  pity * pity  -  pixy * pixy  +  t2 * (pitn * pitn  -  pixn * pixn  -  piyn * piyn))));

				fprintf(RpiT_inverse, "%.2f\t%.2f\t%.2f\t%.6f\n", x, y, z, pi_mag / (sqrt(2.) * pt));
			}
		}
	}
	fclose(RpiT_inverse);
#endif
}


void output_viscous_bjorken(const hydro_variables * const __restrict__ q, const precision * const e, precision e0, precision t, lattice_parameters lattice, hydro_parameters hydro)
{
#ifndef ANISO_HYDRO
	FILE *energy, *plptratio, *shear, *bulk, *db;

	energy = fopen("output/e_e0.dat", "a");
	plptratio = fopen("output/pl_pt.dat", "a");
	shear = fopen("output/shear_peq.dat", "a");
	bulk = fopen("output/bulk_peq.dat", "a");
	db = fopen("output/db_peq.dat", "a");

	int s = central_index(lattice);

	precision e_s = e[s];

	equation_of_state_new eos(e_s, hydro.conformal_eos_prefactor);
	precision p = eos.equilibrium_pressure();
	precision T = eos.T;
	precision mdmde = eos.mdmde_quasi();
	precision m = T * eos.z_quasi();
	precision zeta = ((e_s + p) / T) * zeta_over_s(T, hydro);
	precision taubulk = zeta / eos.beta_bulk();

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

	fprintf(energy,		"%.8f\t%.8e\n", t, e_s / e0);
	fprintf(plptratio,	"%.8f\t%.8e\n", t, pl / pt);
	fprintf(shear,		"%.8f\t%.8e\n", t, 2. * (pl - pt) / (3.*p));
	fprintf(bulk,		"%.8f\t%.8e\n", t, Pi / p);
	fprintf(db,			"%.8f\t%.8e\n", t, -3. * taubulk * mdmde * (e_s + p) * Pi / (t * m * m * p));

	fclose(energy);
	fclose(plptratio);
	fclose(shear);
	fclose(bulk);
	fclose(db);

#ifdef E_CHECK
	FILE *energy_check;
	energy_check = fopen("output/e_check_e0.dat", "a");
	precision e_check = q[s].e_check;

	fprintf(energy_check, "%.8f\t%.8e\n", t, e_check / e0);
	fclose(energy_check);
#endif

#endif
}


void output_aniso_bjorken(const hydro_variables * const __restrict__ q, const precision * const e, double e0, double t, lattice_parameters lattice, hydro_parameters hydro)
{
#ifdef ANISO_HYDRO
	FILE *energy, *plptratio;
	FILE *bulk, *shear;

	energy = fopen("output/e_e0.dat", "a");
	plptratio = fopen("output/pl_pt.dat", "a");
	bulk = fopen("output/bulk_peq.dat", "a");
	shear = fopen("output/shear_peq.dat", "a");

	int s = central_index(lattice);

	precision e_s = e[s];

	equation_of_state_new eos(e_s, hydro.conformal_eos_prefactor);
	precision p = eos.equilibrium_pressure();
	precision beq = eos.equilibrium_mean_field();

	precision pl = q[s].pl;
	precision pt = q[s].pt;

	fprintf(energy, 	"%.8f\t%.8e\n", t, e_s / e0);
	fprintf(plptratio, 	"%.8f\t%.8e\n", t, pl / pt);
	fprintf(bulk, 		"%.8f\t%.8e\n", t, (pl + 2.*pt) / (3.*p) - 1.);
	fprintf(shear, 		"%.8f\t%.8e\n", t, 2.*(pl - pt) / (3.*p));

	fclose(energy);
	fclose(plptratio);
	fclose(shear);
	fclose(bulk);

#ifdef B_FIELD
	FILE *bfield;
	bfield = fopen("output/db_peq.dat", "a");
	precision b = q[s].b;

	fprintf(bfield, "%.8f\t%.8e\n", t, (b - beq) / p);		// db / peq
	fclose(bfield);
#endif


#ifdef E_CHECK
	FILE *energy_check;
	energy_check = fopen("output/e_check_e0.dat", "a");
	precision e_check = q[s].e_check;

	fprintf(energy_check, "%.8f\t%.8e\n", t, e_check / e0);
	fclose(energy_check);
#endif

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

				fprintf(energy, 	"%.2f\t%.2f\t%.2f\t%.6f\n", x, y, z, e_s * hbarc);	// GeV / fm^3
				fprintf(plptratio,	"%.2f\t%.2f\t%.2f\t%.6f\n", x, y, z, pl / pt);
				fprintf(uxplot, 	"%.2f\t%.2f\t%.2f\t%.6f\n", x, y, z, ux);
				fprintf(urplot, 	"%.2f\t%.2f\t%.2f\t%.6f\n", x, y, z, ur);
			}
		}
	}
	fclose(energy);
	fclose(plptratio);
	fclose(uxplot);
	fclose(urplot);
}


void output_hydro(const hydro_variables * const __restrict__ q, const fluid_velocity  * const __restrict__ u, const precision * const e, double t, lattice_parameters lattice, hydro_parameters hydro)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision dx = lattice.lattice_spacing_x;
	precision dy = lattice.lattice_spacing_y;
	precision dz = lattice.lattice_spacing_eta;

	FILE *energy;			// energy density [GeV/fm^3]
	FILE *plptratio;		// plpt ratio
	FILE *uxplot;			// ux
	FILE *shear;			// pressure anisotropy
	FILE *bulk;				// bulk pressure

#ifdef E_CHECK
	FILE *energy_check;		// energy density evolved separately with KT algorithm [GeV/fm^3]
#endif

#ifdef B_FIELD
	FILE *mean_field;
#endif

	char fname1[255];
	char fname2[255];
	char fname3[255];
	char fname4[255];
	char fname5[255];
	char fname6[255];
	char fname7[255];

	sprintf(fname1, "output/e_%.3f.dat", t);
	sprintf(fname2, "output/plpt_%.3f.dat", t);
	sprintf(fname3, "output/ux_%.3f.dat", t);
	sprintf(fname4, "output/shear_peq_%.3f.dat", t);
	sprintf(fname5, "output/bulk_peq_%.3f.dat", t);
	sprintf(fname6, "output/e_check_%.3f.dat", t);
	sprintf(fname7, "output/db_peq_%.3f.dat", t);

	energy 		= fopen(fname1, "w");
	plptratio 	= fopen(fname2, "w");
	uxplot    	= fopen(fname3, "w");
	shear 	  	= fopen(fname4, "w");
	bulk 		= fopen(fname5, "w");
#ifdef E_CHECK
	energy_check= fopen(fname6, "w");
#endif

#ifdef B_FIELD
	mean_field 	= fopen(fname7, "w");
#endif

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

			#ifndef BOOST_INVARIANT
				precision un = u[s].un;
			#else
				precision un = 0;
			#endif

				precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t2 * un * un);

				equation_of_state_new eos(e_s, hydro.conformal_eos_prefactor);
				precision p = eos.equilibrium_pressure();

			#ifdef ANISO_HYDRO
				precision pl = q[s].pl;
				precision pt = q[s].pt;
			#else
				precision utperp = sqrt(1.  +  ux * ux  +  uy * uy);
				precision zt = t * un / utperp;
				precision zn = ut / t / utperp;

			#ifdef PIMUNU
				precision pitt = q[s].pitt;
			#ifndef BOOST_INVARIANT
				precision pitn = q[s].pitn;
			#else
				precision pitn = 0;
			#endif
				precision pinn = q[s].pinn;
			#else
				precision pitt = 0;
				precision pitn = 0;
				precision pinn = 0;
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

				fprintf(energy, 	"%.2f\t%.2f\t%.2f\t%.6e\n", x, y, z, e_s * hbarc);
				fprintf(plptratio,	"%.2f\t%.2f\t%.2f\t%.6f\n", x, y, z, pl / pt);
				fprintf(uxplot, 	"%.2f\t%.2f\t%.2f\t%.6f\n", x, y, z, ux);
				fprintf(shear,	 	"%.2f\t%.2f\t%.2f\t%.6f\n", x, y, z, 2. * (pl - pt) / (3. * p));
				fprintf(bulk,	 	"%.2f\t%.2f\t%.2f\t%.6f\n", x, y, z, (pl + 2.*pt) / (3. * p)  -  1.);
			#ifdef E_CHECK
				fprintf(energy_check, "%.2f\t%.2f\t%.2f\t%.6e\n", x, y, z, q[s].e_check * hbarc);
			#endif

			#ifdef B_FIELD
				precision db = q[s].b - eos.equilibrium_mean_field();
				fprintf(mean_field, "%.2f\t%.2f\t%.2f\t%.6e\n", x, y, z, db / p);
			#endif
			}
		}
	}
	fclose(energy);
	fclose(plptratio);
	fclose(uxplot);
	fclose(shear);
	fclose(bulk);

#ifdef E_CHECK
	fclose(energy_check);
#endif

#ifdef B_FIELD
	fclose(mean_field);
#endif
}

void output_Tmunu_violations(const float * const __restrict__ Tmunu_violations, double t, lattice_parameters lattice)
{
#ifdef MONITOR_TTAUMU
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision dx = lattice.lattice_spacing_x;
	precision dy = lattice.lattice_spacing_y;
	precision dz = lattice.lattice_spacing_eta;

	FILE *violations;
	char fname[255];
	sprintf(fname, "output/Tmunu_violations_%.3f.dat", t);

	violations = fopen(fname, "w");

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

				fprintf(violations, "%.3f\t%.3f\t%.3f\t%f\n", x, y, z, Tmunu_violations[s]);
			}
		}
	}
	fclose(violations);
#endif
}


void output_regulations(const int * const __restrict__ regulation, double t, lattice_parameters lattice)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision dx = lattice.lattice_spacing_x;
	precision dy = lattice.lattice_spacing_y;
	precision dz = lattice.lattice_spacing_eta;

	FILE *regulations;
	char fname[255];

#ifdef ANISO_HYDRO
#ifdef LATTICE_QCD
	sprintf(fname, "output/aniso_regulations_%.3f.dat", t);
#endif
#else
	sprintf(fname, "output/viscous_regulations_%.3f.dat", t);
#endif

	regulations = fopen(fname, "w");

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

			#ifdef ANISO_HYDRO
			#ifdef LATTICE_QCD
				fprintf(regulations, "%.3f\t%.3f\t%.3f\t%d\n", x, y, z, aniso_regulation[s]);
			#endif
			#else
				fprintf(regulations, "%.3f\t%.3f\t%.3f\t%d\n", x, y, z, viscous_regulation[s]);
			#endif
			}
		}
	}
	fclose(regulations);
}


void output_dynamical_variables(double t, double dt, lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro)
{
	int initial_type = initial.initial_condition_type;

	if(initial_type == 1)
	{
		precision T0 = initial.initialCentralTemperatureGeV;
		precision e0 = equilibrium_energy_density_new(T0 / hbarc, hydro.conformal_eos_prefactor);

	#ifdef ANISO_HYDRO
		output_aniso_bjorken(q, e, e0, t, lattice, hydro);
	#else
		output_viscous_bjorken(q, e, e0, t, lattice, hydro);
	#endif
	}
	else if(initial_type == 2 || initial_type == 3)
	{
		output_gubser(q, u, e, t, lattice);
		output_residual_shear_validity(q, u, e, t, lattice, hydro);
	}
	else
	{
		output_hydro(q, u, e, t, lattice, hydro);
		output_residual_shear_validity(q, u, e, t, lattice, hydro);

	#ifdef MONITOR_TTAUMU
		output_Tmunu_violations(Tmunu_violations, t, lattice);
	#endif

	#ifdef ANISO_HYDRO
	#ifdef LATTICE_QCD
		output_regulations(aniso_regulation, t, lattice);
	#endif
	#else
		output_regulations(viscous_regulation, t, lattice);
	#endif
	}
}


void output_semi_analytic_solution_if_any(lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro)
{
	switch(initial.initial_condition_type)
	{
		case 1:		// Bjorken
		{
		#ifdef ANISO_HYDRO
			printf("\nRunning semi-analytic anisotropic Bjorken solution...\n");
			run_semi_analytic_aniso_bjorken(lattice, initial, hydro);
		#else
			printf("\nRunning semi-analytic viscous Bjorken solution...\n");
			run_semi_analytic_viscous_bjorken(lattice, initial, hydro);
		#endif
			break;
		}
		case 3:		// anisotropic Gubser
		{

		#ifdef ANISO_HYDRO
			printf("\nRunning semi-analytic anisotropic Gubser solution...\n");
			double T0_hat = run_semi_analytic_aniso_gubser(lattice, initial, hydro);
		#else
			printf("\nRunning semi-analytic viscous Gubser solution...\n");
			double T0_hat = run_semi_analytic_viscous_gubser(lattice, initial, hydro);
		#endif
			break;
		}
		default:
		{
			printf("\nNo semi-analytic solution to run...\n");
			break;
		}
	}
}













