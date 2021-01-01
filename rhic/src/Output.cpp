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

inline precision central_derivative(const precision * const __restrict__ f, int n, precision dx)
{
	return (f[n + 1] - f[n]) / (2. * dx);		// f[n] = fm  |	 f[n+1] = fp
}


void output_mean_field(const hydro_variables * const __restrict__ q, const fluid_velocity * const __restrict__ u, const fluid_velocity  * const __restrict__ up, const precision * const e, double t, double dt_prev, lattice_parameters lattice, hydro_parameters hydro)
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

	fprintf(db_peq,    "%d\n%d\n%d\n", nx, ny, nz);
	fprintf(dbasy_peq, "%d\n%d\n%d\n", nx, ny, nz);
#endif
	db2_peq 	= fopen(fname3, "w");
	fprintf(db2_peq, "%d\n%d\n%d\n", nx, ny, nz);

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
				precision mdmde = eos.mdmde_quasi();		// why use m dmde instead of dmde?
				precision betabulk = eos.beta_bulk();
				precision zetas = zeta_over_s(T, hydro);
				precision zeta = zetas * (e_s + p) / T;
				precision taubulk = zeta / betabulk;

				precision ux = u[s].ux;
				precision uy = u[s].uy;
				precision un = 0;

				precision ux_p = up[s].ux;
				precision uy_p = up[s].uy;
				precision un_p = 0;

			#ifndef BOOST_INVARIANT
				un = u[s].un;
				un_p = up[s].un;
			#endif

				precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t2 * un * un);
				precision utperp = sqrt(1.  +  ux * ux  +  uy * uy);

				precision vx = ux / ut;
				precision vy = uy / ut;
				precision vn = un / ut;

				precision zt = t * un / utperp;
				precision zn = ut / t / utperp;

				precision zt2 = zt * zt;
				precision zn2 = zn * zn;
				precision ztzn = zt * zn;

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

				precision dux_dt = (ux - ux_p) / dt_prev;
				precision dux_dx = central_derivative(ui1, 0, dx);
				precision dux_dy = central_derivative(uj1, 0, dy);
				precision dux_dn = central_derivative(uk1, 0, dn);

				precision duy_dt = (uy - uy_p) / dt_prev;
				precision duy_dx = central_derivative(ui1, 2, dx);
				precision duy_dy = central_derivative(uj1, 2, dy);
				precision duy_dn = central_derivative(uk1, 2, dn);

				precision dun_dt = (un - un_p) / dt_prev;
				precision dun_dx = central_derivative(ui1, 4, dx);
				precision dun_dy = central_derivative(uj1, 4, dy);
				precision dun_dn = central_derivative(uk1, 4, dn);

				precision dut_dt = vx * dux_dt  +  vy * duy_dt  +  t2 * vn * dun_dt  +  t * vn * un;
				precision dut_dx = vx * dux_dx  +  vy * duy_dx  +  t2 * vn * dun_dx;
				precision dut_dy = vx * dux_dy  +  vy * duy_dy  +  t2 * vn * dun_dy;
				precision dut_dn = vx * dux_dn  +  vy * duy_dn  +  t2 * vn * dun_dn;

				precision theta = dut_dt  +  dux_dx  +  duy_dy  +  dun_dn  +  ut / t;
				precision thetaL = - zt2 * dut_dt  +  t2 * zn2 * dun_dn  +  ztzn * (t2 * dun_dt  -  dut_dn)  +  t * zn2 * ut;
				precision thetaT = theta  -  thetaL;

			#ifdef ANISO_HYDRO
				precision pl = q[s].pl;
				precision pt = q[s].pt;
				precision b = q[s].b;
				precision Pi = (pl + 2.*pt) / 3. - p;

				precision edot = -(e_s + pl) * thetaL  -  (e_s + pt) * thetaT;		// ignores residual shear

				precision mdot = mdmde * edot / mass;

				precision dbasy = 3. * taubulk * mdot * Pi / (mass - 4.*taubulk*mdot);

				fprintf(db_peq, 	"%.2f\t%.2f\t%.2f\t%.4e\n", x, y, z, (b - beq) / p);
				fprintf(dbasy_peq, 	"%.2f\t%.2f\t%.2f\t%.4e\n", x, y, z, dbasy / p);

			#else
			#ifdef PI
				precision Pi = q[s].Pi;
			#else
				precision Pi = 0;
			#endif
			#endif

				precision db2 = -3.* taubulk * mdmde * (e_s + p) * Pi * theta / (mass * mass);

				fprintf(db2_peq, "%.2f\t%.2f\t%.2f\t%.4e\n", x, y, z, db2 / p);
			}
		}
	}
#ifdef ANISO_HYDRO
	fclose(db_peq);
	fclose(dbasy_peq);
#endif
	fclose(db2_peq);
#endif
}


void output_residual_transverse_shear_validity(const hydro_variables * const __restrict__ q, const fluid_velocity * const __restrict__ u, const precision * const e, double t, lattice_parameters lattice, hydro_parameters hydro)
{
#ifdef PIMUNU
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision dx = lattice.lattice_spacing_x;
	precision dy = lattice.lattice_spacing_y;
	precision dn = lattice.lattice_spacing_eta;

	FILE *piperp_pt;
	char fname1[255];
	sprintf(fname1, "output/piperp_pt_%.3f.dat", t);
	piperp_pt = fopen(fname1, "w");

	fprintf(piperp_pt, "%d\n%d\n%d\n", nx, ny, nz);

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

				fprintf(piperp_pt, "%.2f\t%.2f\t%.2f\t%.4e\n", x, y, z, pi_mag / (sqrt(2.) * pt));
			}
		}
	}
	fclose(piperp_pt);
#endif
}


void output_residual_longitudinal_shear_validity(const hydro_variables * const __restrict__ q, const fluid_velocity * const __restrict__ u, const precision * const e, double t, lattice_parameters lattice, hydro_parameters hydro)
{
#ifdef BOOST_INVARIANT
	return;
#endif

#ifdef ANISO_HYDRO
#ifndef WTZMU
	return;
#endif
#else
#ifndef PIMUNU
	return;
#endif
#endif

	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision dx = lattice.lattice_spacing_x;
	precision dy = lattice.lattice_spacing_y;
	precision dn = lattice.lattice_spacing_eta;

	FILE *WTz_pl_2pt;
	char fname[255];
	sprintf(fname, "output/WTz_pl_2pt_%.3f.dat", t);
	WTz_pl_2pt = fopen(fname, "w");

	fprintf(WTz_pl_2pt, "%d\n%d\n%d\n", nx, ny, nz);

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

			#ifdef ANISO_HYDRO
				precision pl = q[s].pl;
				precision pt = q[s].pt;
			#ifdef WTZMU
				precision WtTz = q[s].WtTz;
				precision WxTz = q[s].WxTz;
				precision WyTz = q[s].WyTz;
				precision WnTz = q[s].WnTz;
			#else
				precision WtTz = 0;
				precision WxTz = 0;
				precision WyTz = 0;
				precision WnTz = 0;
			#endif

			#else

				precision e_s = e[s];
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

			#ifdef PIMUNU
				precision pitt = q[s].pitt;
				precision pitx = q[s].pitx;
				precision pity = q[s].pity;
				precision pitn = 0;
				precision pixx = q[s].pixx;
				precision pixy = q[s].pixy;
				precision pixn = 0;
				precision piyy = q[s].piyy;
				precision piyn = 0;
				precision pinn = q[s].pinn;
			#ifndef BOOST_INVARIANT
				pitn = q[s].pitn;
				pixn = q[s].pixn;
				piyn = q[s].piyn;
			#endif
			#else
				precision pitt = 0;
				precision pitx = 0;
				precision pity = 0;
				precision pitn = 0;
				precision pixx = 0;
				precision pixy = 0;
				precision pixn = 0;
				precision piyy = 0;
				precision piyn = 0;
				precision pinn = 0;
			#endif

				// pl = zmu.z.nu.Tmunu
				precision pl = P  +  zt * zt * pitt  +  t4 * zn * zn * pinn  +  2. * t2 * zt * zn * pitn;
				// pt = -Ximunu.Tmunu / 2
				precision pt = P  -  (zt * zt * pitt  +  t4 * zn * zn * pinn  +  2. * t2 * zt * zn * pitn) / 2.;


				// z_\nu . pi^{\nu\alpha}
				precision pizt = zt * pitt  -  t2 * zn * pitn;
				precision pizx = zt * pitx  -  t2 * zn * pixn;
				precision pizy = zt * pity  -  t2 * zn * piyn;
				precision pizn = zt * pitn  -  t2 * zn * pinn;

				transverse_projection Xi(ut, ux, uy, un, zt, zn, t2);

				precision Xitt = Xi.Xitt;
				precision Xitx = Xi.Xitx;
				precision Xity = Xi.Xity;
				precision Xitn = Xi.Xitn;
				precision Xixx = Xi.Xixx;
				precision Xixy = Xi.Xixy;
				precision Xixn = Xi.Xixn;
				precision Xiyy = Xi.Xiyy;
				precision Xiyn = Xi.Xiyn;
				precision Xinn = Xi.Xinn;

				// WTz^\mu = - \Xi^\mu_\alpha . z_\nu . pi^{\nu\alpha}
				precision WtTz = - Xitt * pizt  +  Xitx * pizx  +  Xity * pizy  +  t2 * Xitn * pizn;
				precision WxTz = - Xitx * pizt  +  Xixx * pizx  +  Xixy * pizy  +  t2 * Xixn * pizn;
				precision WyTz = - Xity * pizt  +  Xixy * pizx  +  Xiyy * pizy  +  t2 * Xiyn * pizn;
				precision WnTz = - Xitn * pizt  +  Xixn * pizx  +  Xiyn * pizy  +  t2 * Xinn * pizn;
			#endif

				precision W_mag = sqrt(2. * fabs(WtTz * WtTz  -  WxTz * WxTz  -  WyTz * WyTz  -  t2 * WnTz * WnTz));

				fprintf(WTz_pl_2pt, "%.2f\t%.2f\t%.2f\t%.4e\n", x, y, z, W_mag / sqrt(pl * pl  +  2. * pt * pt));
			}
		}
	}
	fclose(WTz_pl_2pt);
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

	fprintf(energy,    "%d\n%d\n%d\n", nx, ny, nz);
	fprintf(plptratio, "%d\n%d\n%d\n", nx, ny, nz);
	fprintf(uxplot,    "%d\n%d\n%d\n", nx, ny, nz);
	fprintf(urplot,    "%d\n%d\n%d\n", nx, ny, nz);

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

#ifndef BOOST_INVARIANT
	FILE *unplot;			// un
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
	sprintf(fname7, "output/un_%.3f.dat", t);

	energy       = fopen(fname1, "w");
	plptratio    = fopen(fname2, "w");
	uxplot       = fopen(fname3, "w");

	fprintf(energy,    "%d\n%d\n%d\n", nx, ny, nz);
	fprintf(plptratio, "%d\n%d\n%d\n", nx, ny, nz);
	fprintf(uxplot,    "%d\n%d\n%d\n", nx, ny, nz);

#ifndef CONFORMAL_EOS
	shear        = fopen(fname4, "w");
	bulk         = fopen(fname5, "w");

	fprintf(shear, "%d\n%d\n%d\n", nx, ny, nz);
	fprintf(bulk,  "%d\n%d\n%d\n", nx, ny, nz);
#endif
#ifdef E_CHECK
	energy_check = fopen(fname6, "w");

	fprintf(energy_check, "%d\n%d\n%d\n", nx, ny, nz);
#endif
#ifndef BOOST_INVARIANT
	unplot       = fopen(fname7, "w");

	fprintf(unplot, "%d\n%d\n%d\n", nx, ny, nz);
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

				fprintf(energy, 	"%.2f\t%.2f\t%.2f\t%.4e\n", x, y, z, e_s * hbarc);
				fprintf(plptratio,	"%.2f\t%.2f\t%.2f\t%.4e\n", x, y, z, pl / pt);
				fprintf(uxplot, 	"%.2f\t%.2f\t%.2f\t%.4e\n", x, y, z, ux);
			#ifndef BOOST_INVARIANT
				fprintf(unplot, 	"%.2f\t%.2f\t%.2f\t%.4e\n", x, y, z, un);
			#endif
			#ifndef CONFORMAL_EOS
				fprintf(shear,	 	"%.2f\t%.2f\t%.2f\t%.4e\n", x, y, z, 2. * (pl - pt) / (3. * p));
				fprintf(bulk,	 	"%.2f\t%.2f\t%.2f\t%.4e\n", x, y, z, (pl + 2.*pt) / (3. * p)  -  1.);
			#endif
			#ifdef E_CHECK
				fprintf(energy_check, "%.2f\t%.2f\t%.2f\t%.4e\n", x, y, z, q[s].e_check * hbarc);
			#endif
			}
		}
	}
	fclose(energy);
	fclose(plptratio);
	fclose(uxplot);
#ifndef BOOST_INVARIANT
	fclose(unplot);
#endif
#ifndef CONFORMAL_EOS
	fclose(shear);
	fclose(bulk);
#endif
#ifdef E_CHECK
	fclose(energy_check);
#endif
}


void output_freezeout_slice_x(double t, lattice_parameters lattice, hydro_parameters hydro)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;
	precision dx = lattice.lattice_spacing_x;

	precision T_switch = hydro.freezeout_temperature_GeV;
	precision e_switch = equilibrium_energy_density_new(T_switch / hbarc, hydro.conformal_eos_prefactor);

	FILE * surface_slice_x;	// tau-x slice of freezeout surface (only freezeout cell positions)
	FILE * Rinv_pi_x;		// tau-x slice of (residual) shear inverse Reynolds number
	FILE * Rinv_bulk_x;		// tau-x slice of bulk inverse Reynolds number
	FILE * Rinv_dB_x;		// tau-x slice of dB/Peq inverse Reynolds number
	FILE * reg_viscous_x;	// tau-x slice of residual shear or shear + bulk regulations
	FILE * reg_dB_x;		// tau-x slice of dB regulation;

	char fname1[255];
	char fname2[255];
	char fname3[255];
	char fname4[255];
	char fname5[255];
	char fname6[255];

	sprintf(fname1, "output/surface_slice_x.dat");
	sprintf(fname2, "output/Rinv_shear_x.dat");
	sprintf(fname3, "output/Rinv_bulk_x.dat");
	sprintf(fname4, "output/Rinv_dB_x.dat");
	sprintf(fname5, "output/regulate_viscous_x.dat");
	sprintf(fname6, "output/regulate_dB_x.dat");

	surface_slice_x  = fopen(fname1, "a");
	Rinv_pi_x        = fopen(fname2, "a");
	Rinv_bulk_x      = fopen(fname3, "a");
	Rinv_dB_x        = fopen(fname4, "a");
	reg_viscous_x    = fopen(fname5, "a");
	reg_dB_x         = fopen(fname6, "a");

	precision t2 = t * t;
	precision t4 = t2 * t2;

	double xL;	  // left-side position and energy density
	double eL;

	for(int i = 2; i < nx + 2; i++)
	{
		double x = (i - 2. - (nx - 1.)/2.) * dx;

		int j = 2  +  (ny - 1) / 2;
		int k = 2  +  (nz - 1) / 2;

		int s = linear_column_index(i, j, k, nx + 4, ny + 4);

		precision e_s = e[s];

		// search for freezeout cells along x-direction
		if(i > 2)
		{
			// check if crossed freezeout energy density (freezeout surface construction)
			if((eL >= e_switch && e_s <= e_switch) || (eL <= e_switch && e_s >= e_switch))
			{
				double x_cell = xL  +  (e_switch - eL) * dx / (e_s - eL);	// linear interpolation

				fprintf(surface_slice_x, "%.3f\t%.4e\n", x_cell, t);
			}
		}

    	xL = x;		// update left side
    	eL = e_s;

    	equation_of_state_new eos(e_s, hydro.conformal_eos_prefactor);
		precision p = eos.equilibrium_pressure();

    #ifdef B_FIELD
    	precision b = q[s].b;
    	precision beq = eos.equilibrium_mean_field();
    	precision db =  b - beq;
    #endif

    	precision pitt = 0;
		precision pitx = 0;
		precision pity = 0;
		precision pitn = 0;
		precision pixx = 0;
		precision pixy = 0;
		precision pixn = 0;
		precision piyy = 0;
		precision piyn = 0;
		precision pinn = 0;

	#ifdef PIMUNU
		pitt = q[s].pitt;
		pitx = q[s].pitx;
		pity = q[s].pity;
		pitn = 0;
		pixx = q[s].pixx;
		pixy = q[s].pixy;
		pixn = 0;
		piyy = q[s].piyy;
		piyn = 0;
		pinn = 0;

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
	#endif

		precision WtTz = 0;
		precision WxTz = 0;
		precision WyTz = 0;
		precision WnTz = 0;

	#ifdef WTZMU
		WtTz = q[s].WtTz;
		WxTz = q[s].WxTz;
		WyTz = q[s].WyTz;
		WnTz = q[s].WnTz;
	#endif

	#ifdef ANISO_HYDRO
		precision pi_mag = sqrt(fabs(-2. * (WtTz * WtTz  -  WxTz * WxTz  -  WyTz * WyTz  -  t2 * WnTz * WnTz)  +  pitt * pitt  +  pixx * pixx  +  piyy * piyy  +  t4 * pinn * pinn  -  2. * (pitx * pitx  +  pity * pity  -  pixy * pixy  +  t2 * (pitn * pitn  -  pixn * pixn  -  piyn * piyn))));
	#else
		precision pi_mag = sqrt(fabs(pitt * pitt  +  pixx * pixx  +  piyy * piyy  +  t4 * pinn * pinn  -  2. * (pitx * pitx  +  pity * pity  -  pixy * pixy  +  t2 * (pitn * pitn  -  pixn * pixn  -  piyn * piyn))));
	#endif


	#ifdef ANISO_HYDRO			// aniso hydro slice:
		double pl = q[s].pl;
		double pt = q[s].pt;

		// sqrt((piT.piT-2.WTz.WTz)/(pl.pl+2.pt.pt))
		fprintf(Rinv_pi_x, "%.3f\t%.4e\t%.4e\n", x, t, pi_mag / sqrt(pl*pl + 2.*pt*pt));

	#ifdef B_FIELD
		// |dB| / p
		fprintf(Rinv_dB_x, "%.3f\t%.4e\t%.4e\n", x, t, fabs(db / p));
	#endif

	#ifdef MONITOR_B
		// dB regulation
		fprintf(reg_dB_x, "%.3f\t%.4e\t%d\n", x, t, b_regulation[s]);
	#endif


	#else 						// viscous hydro slice:

		// sqrt(pi.pi/3) / p
		fprintf(Rinv_pi_x, "%.3f\t%.4e\t%.4e\n", x, t, pi_mag / (sqrt(3.) * p));

	#ifdef PI
		// |Pi| / p
		fprintf(Rinv_bulk_x, "%.3f\t%.4e\t%.4e\n", x, t, fabs(q[s].Pi) / p);
	#endif

	#endif


	#ifdef MONITOR_REGULATIONS
		// residual shear or shear + bulk regulation
		fprintf(reg_viscous_x, "%.3f\t%.4e\t%d\n", x, t, viscous_regulation[s]);
	#endif
	}

	fclose(surface_slice_x);
	fclose(Rinv_pi_x);
	fclose(Rinv_bulk_x);
	fclose(Rinv_dB_x);
	fclose(reg_viscous_x);
	fclose(reg_dB_x);
}


void output_freezeout_slice_z(double t, lattice_parameters lattice, hydro_parameters hydro)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;
	precision dz = lattice.lattice_spacing_eta;

	precision T_switch = hydro.freezeout_temperature_GeV;
	precision e_switch = equilibrium_energy_density_new(T_switch / hbarc, hydro.conformal_eos_prefactor);

	FILE * surface_slice_z;	// tau-eta slice of freezeout surface (only freezeout cell positions)
	FILE * Rinv_pi_z;		// tau-eta slice of (residual) shear inverse Reynolds number
	FILE * Rinv_bulk_z;		// tau-eta slice of bulk inverse Reynolds number
	FILE * Rinv_dB_z;		// tau-eta slice of dB/Peq inverse Reynolds number
	FILE * reg_viscous_z;	// tau-eta slice of residual shear or shear + bulk regulations
	FILE * reg_dB_z;		// tau-eta slice of dB regulation;

	char fname1[255];
	char fname2[255];
	char fname3[255];
	char fname4[255];
	char fname5[255];
	char fname6[255];

	sprintf(fname1, "output/surface_slice_z.dat");
	sprintf(fname2, "output/Rinv_shear_z.dat");
	sprintf(fname3, "output/Rinv_bulk_z.dat");
	sprintf(fname4, "output/Rinv_dB_z.dat");
	sprintf(fname5, "output/regulate_viscous_z.dat");
	sprintf(fname6, "output/regulate_dB_z.dat");

	surface_slice_z  = fopen(fname1, "a");
	Rinv_pi_z        = fopen(fname2, "a");
	Rinv_bulk_z      = fopen(fname3, "a");
	Rinv_dB_z        = fopen(fname4, "a");
	reg_viscous_z    = fopen(fname5, "a");
	reg_dB_z         = fopen(fname6, "a");

	precision t2 = t * t;
	precision t4 = t2 * t2;

	double zL;	  // left-side position and energy density
	double eL;

	for(int k = 2; k < nz + 2; k++)
	{
		double z = (k - 2. - (nz - 1.)/2.) * dz;

		int i = 2  +  (nx - 1) / 2;											// x = y = 0
		int j = 2  +  (ny - 1) / 2;
		int s = linear_column_index(i, j, k, nx + 4, ny + 4);

		precision e_s = e[s];

		// search for freezeout cells along eta-direction
		if(k > 2)
		{
			// check if crossed freezeout energy density
			if((eL >= e_switch && e_s <= e_switch) || (eL <= e_switch && e_s >= e_switch))
			{
				double z_cell = zL  +  (e_switch - eL) * dz / (e_s - eL);	// linear interpolation

				fprintf(surface_slice_z, "%.3f\t%.4e\n", z_cell, t);
			}
		}

    	zL = z;		// update left side
    	eL = e_s;

    	equation_of_state_new eos(e_s, hydro.conformal_eos_prefactor);
		precision p = eos.equilibrium_pressure();

    #ifdef B_FIELD
    	precision b = q[s].b;
    	precision beq = eos.equilibrium_mean_field();
    	precision db =  b - beq;
    #endif

    	precision pitt = 0;
		precision pitx = 0;
		precision pity = 0;
		precision pitn = 0;
		precision pixx = 0;
		precision pixy = 0;
		precision pixn = 0;
		precision piyy = 0;
		precision piyn = 0;
		precision pinn = 0;

	#ifdef PIMUNU
		pitt = q[s].pitt;
		pitx = q[s].pitx;
		pity = q[s].pity;
		pixx = q[s].pixx;
		pixy = q[s].pixy;
		piyy = q[s].piyy;

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
	#endif

		precision pi_mag = sqrt(fabs(pitt * pitt  +  pixx * pixx  +  piyy * piyy  +  t4 * pinn * pinn  -  2. * (pitx * pitx  +  pity * pity  -  pixy * pixy  +  t2 * (pitn * pitn  -  pixn * pixn  -  piyn * piyn))));


	#ifdef ANISO_HYDRO			// aniso hydro slice:

	#ifdef B_FIELD
		// |dB| / p
		fprintf(Rinv_dB_z, "%.3f\t%.4e\t%.4e\n", z, t, fabs(db / p));
	#endif

	#ifdef MONITOR_B
		// dB regulation
		fprintf(reg_dB_z, "%.3f\t%.4e\t%d\n", z, t, b_regulation[s]);
	#endif


	#else 						// viscous hydro slice:

		// sqrt(pi.pi/3) / p
		fprintf(Rinv_pi_z, "%.3f\t%.4e\t%.4e\n", z, t, pi_mag / (sqrt(3.) * p));

	#ifdef PI
		// |Pi| / p
		fprintf(Rinv_bulk_z, "%.3f\t%.4e\t%.4e\n", z, t, fabs(q[s].Pi) / p);
	#endif

	#endif


	#ifdef MONITOR_REGULATIONS
		// residual shear or shear + bulk regulation
		fprintf(reg_viscous_z, "%.3f\t%.4e\t%d\n", z, t, viscous_regulation[s]);
	#endif
	}


	// separate loop for residual shear inverse Reynolds number (shift x to ~ right-half of grid)
#ifdef ANISO_HYDRO
	for(int k = 2; k < nz + 2; k++)
	{
		double z = (k - 2. - (nz - 1.)/2.) * dz;

		int i = 2  +  (3 * (nx - 1)) / 4;	// right half
		int j = 2  +  (ny - 1) / 2;
		int s = linear_column_index(i, j, k, nx + 4, ny + 4);

    	precision pitt = 0;
		precision pitx = 0;
		precision pity = 0;
		precision pitn = 0;
		precision pixx = 0;
		precision pixy = 0;
		precision pixn = 0;
		precision piyy = 0;
		precision piyn = 0;
		precision pinn = 0;

	#ifdef PIMUNU
		pitt = q[s].pitt;
		pitx = q[s].pitx;
		pity = q[s].pity;
		pixx = q[s].pixx;
		pixy = q[s].pixy;
		piyy = q[s].piyy;

	#ifndef BOOST_INVARIANT
		pitn = q[s].pitn;
		pixn = q[s].pixn;
		piyn = q[s].piyn;
		pinn = q[s].pinn;
	#endif
	#endif

		precision WtTz = 0;
		precision WxTz = 0;
		precision WyTz = 0;
		precision WnTz = 0;

	#ifdef WTZMU
		WtTz = q[s].WtTz;
		WxTz = q[s].WxTz;
		WyTz = q[s].WyTz;
		WnTz = q[s].WnTz;
	#endif

		precision pi_mag = sqrt(fabs(-2. * (WtTz * WtTz  -  WxTz * WxTz  -  WyTz * WyTz  -  t2 * WnTz * WnTz)  +  pitt * pitt  +  pixx * pixx  +  piyy * piyy  +  t4 * pinn * pinn  -  2. * (pitx * pitx  +  pity * pity  -  pixy * pixy  +  t2 * (pitn * pitn  -  pixn * pixn  -  piyn * piyn))));

		precision pl = q[s].pl;
		precision pt = q[s].pt;

		// sqrt((piT.piT-2.WTz.WTz)/(pl.pl+2.pt.pt))
		fprintf(Rinv_pi_z, "%.3f\t%.4e\t%.4e\n", z, t, pi_mag / sqrt(pl*pl + 2.*pt*pt));
	}
#endif

	fclose(surface_slice_z);
	fclose(Rinv_pi_z);
	fclose(Rinv_bulk_z);
	fclose(Rinv_dB_z);
	fclose(reg_viscous_z);
	fclose(reg_dB_z);
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

	fprintf(violations, "%d\n%d\n%d\n", nx, ny, nz);

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

				fprintf(violations, "%.2f\t%.2f\t%.2f\t%.4e\n", x, y, z, hbarc * Tmunu_violations[s]);
			}
		}
	}
	fclose(violations);
#endif
}


void output_plpt_regulations(const int * const __restrict__ plpt_regulation, double t, lattice_parameters lattice)
{
#ifdef MONITOR_PLPT
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision dx = lattice.lattice_spacing_x;
	precision dy = lattice.lattice_spacing_y;
	precision dz = lattice.lattice_spacing_eta;

	FILE *regulation;
	char fname[255];
	sprintf(fname, "output/plpt_regulation_%.3f.dat", t);

	regulation = fopen(fname, "w");

	fprintf(regulation, "%d\n%d\n%d\n", nx, ny, nz);

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

				fprintf(regulation, "%.2f\t%.2f\t%.2f\t%d\n", x, y, z, plpt_regulation[s]);
			}
		}
	}
	fclose(regulation);
#endif
}


void output_b_regulations(const int * const __restrict__ b_regulation, double t, lattice_parameters lattice)
{
#ifdef MONITOR_B
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision dx = lattice.lattice_spacing_x;
	precision dy = lattice.lattice_spacing_y;
	precision dz = lattice.lattice_spacing_eta;

	FILE *regulation;
	char fname[255];
	sprintf(fname, "output/b_regulation_%.3f.dat", t);

	regulation = fopen(fname, "w");

	fprintf(regulation, "%d\n%d\n%d\n", nx, ny, nz);

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

				fprintf(regulation, "%.2f\t%.2f\t%.2f\t%d\n", x, y, z, b_regulation[s]);
			}
		}
	}
	fclose(regulation);
#endif
}


void output_viscous_regulations(const int * const __restrict__ viscous_regulation, double t, lattice_parameters lattice)
{
#ifdef MONITOR_REGULATIONS
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision dx = lattice.lattice_spacing_x;
	precision dy = lattice.lattice_spacing_y;
	precision dz = lattice.lattice_spacing_eta;

	FILE *reg_file;
	char fname[255];

	sprintf(fname, "output/viscous_regulations_%.3f.dat", t);

	reg_file = fopen(fname, "w");

	fprintf(reg_file, "%d\n%d\n%d\n", nx, ny, nz);

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

				fprintf(reg_file, "%.2f\t%.2f\t%.2f\t%d\n", x, y, z, viscous_regulation[s]);
			}
		}
	}
	fclose(reg_file);
#endif
}


void output_hydro_simulation(double t, double dt_prev, lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro)
{
	int initial_condition_type = initial.initial_condition_type;

	if(initial_condition_type == 1)
	{
		precision T0 = initial.initialCentralTemperatureGeV;
		precision e0 = equilibrium_energy_density_new(T0 / hbarc, hydro.conformal_eos_prefactor);

	#ifdef ANISO_HYDRO
		output_aniso_bjorken(q, e, e0, t, lattice, hydro);
	#else
		output_viscous_bjorken(q, e, e0, t, lattice, hydro);
	#endif
	}
	else if(initial_condition_type == 2)
	{
		output_gubser(q, u, e, t, lattice);
		output_residual_transverse_shear_validity(q, u, e, t, lattice, hydro);
	}
	else
	{
		output_hydro(q, u, e, t, lattice, hydro);
		output_residual_transverse_shear_validity(q, u, e, t, lattice, hydro);
		output_residual_longitudinal_shear_validity(q, u, e, t, lattice, hydro);

	#ifdef LATTICE_QCD
		output_mean_field(q, u, up, e, t, dt_prev, lattice, hydro);
	#endif

	#ifdef MONITOR_PLPT
		output_plpt_regulations(plpt_regulation, t, lattice);
	#endif

	#ifdef MONITOR_B
		output_b_regulations(b_regulation, t, lattice);
	#endif

	#ifdef MONITOR_TTAUMU
		output_Tmunu_violations(Tmunu_violations, t, lattice);
	#endif

	#ifdef MONITOR_REGULATIONS
		output_viscous_regulations(viscous_regulation, t, lattice);
	#endif
	}
}




