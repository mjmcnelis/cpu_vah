#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../include/FileIO.h"
#include "../include/Parameters.h"
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"
#include "../include/EquationOfState.h"
#include "../include/TransportCoefficients.h"
#include "../include/Projections.h"
#include "../include/NeighborCells.h"
#include "../include/Hydrodynamics.h"
#include "../include/AnisoBjorken.h"


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


void output_residual_gradients(const hydro_variables * const __restrict__ q, const fluid_velocity * const __restrict__ u, const fluid_velocity * const __restrict__ up, const precision * const e, double t, double dt, lattice_parameters lattice, hydro_parameters hydro)
{
#ifdef ANISO_HYDRO
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision dx = lattice.lattice_spacing_x;
	precision dy = lattice.lattice_spacing_y;
	precision dn = lattice.lattice_spacing_eta;

	FILE *Knpi;
	char fname1[255];
	sprintf(fname1, "output/Knpi_%.3f.dat", t);		// tau_pi . sqrt(\sigma_T . \sigma_T)
	Knpi = fopen(fname1, "w");

	precision t2 = t * t;
	precision t4 = t2 * t2;

	int stride_y = nx + 4;
	int stride_z = (nx + 4) * (ny + 4);

	precision ui1[6], uj1[6], uk1[6];
	precision vxi[4], vyj[4], vnk[4];

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

				int simm = s - 2;
				int sim  = s - 1;
				int sip  = s + 1;
				int sipp = s + 2;

				int sjmm = s - 2*stride_y;
				int sjm  = s - stride_y;
				int sjp  = s + stride_y;
				int sjpp = s + 2*stride_y;

				int skmm = s - 2*stride_z;
				int skm  = s - stride_z;
				int skp  = s + stride_z;
				int skpp = s + 2*stride_z;

				get_fluid_velocity_neighbor_cells(u[simm], u[sim], u[sip], u[sipp], u[sjmm], u[sjm], u[sjp], u[sjpp], u[skmm], u[skm], u[skp], u[skpp], ui1, uj1, uk1, vxi, vyj, vnk, t2);

				precision e_s = e[s];
				precision T = effectiveTemperature(e_s);

				precision etabar = eta_over_s(T, hydro);
				precision tau_pi = (5. * etabar) / T;

				precision ux = u[s].ux;
				precision uy = u[s].uy;
				precision un = u[s].un;
				precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t2 * un * un);

				precision ux_p = up[s].ux;
				precision uy_p = up[s].uy;
				precision un_p = up[s].un;

				precision utperp = sqrt(1.  +  ux * ux  +  uy * uy);
				precision zt = t * un / utperp;
				precision zn = ut / t / utperp;

				precision vx = ux / ut;
				precision vy = uy / ut;
				precision vn = un / ut;

				// compute \sigma_T^{\mu\nu}
				precision dux_dt = (ux - ux_p) / dt;
				precision dux_dx = central_derivative(ui1, 0, dx);
				precision dux_dy = central_derivative(uj1, 0, dy);
				precision dux_dn = central_derivative(uk1, 0, dn);

				precision duy_dt = (uy - uy_p) / dt;
				precision duy_dx = central_derivative(ui1, 2, dx);
				precision duy_dy = central_derivative(uj1, 2, dy);
				precision duy_dn = central_derivative(uk1, 2, dn);

				precision dun_dt = (un - un_p) / dt;
				precision dun_dx = central_derivative(ui1, 4, dx);
				precision dun_dy = central_derivative(uj1, 4, dy);
				precision dun_dn = central_derivative(uk1, 4, dn);

				precision dut_dt = vx * dux_dt  +  vy * duy_dt  +  t2 * vn * dun_dt  +  t * vn * un;
				precision dut_dx = vx * dux_dx  +  vy * duy_dx  +  t2 * vn * dun_dx;
				precision dut_dy = vx * dux_dy  +  vy * duy_dy  +  t2 * vn * dun_dy;
				precision dut_dn = vx * dux_dn  +  vy * duy_dn  +  t2 * vn * dun_dn;

				precision sTtt = dut_dt;
				precision sTtx = (dux_dt  -  dut_dx) / 2.;
				precision sTty = (duy_dt  -  dut_dy) / 2.;
				precision sTtn = (dun_dt  +  un / t  -  (dut_dn  +  t * un) / t2) / 2.;
				precision sTxx = - dux_dx;
				precision sTxy = - (dux_dy  +  duy_dx) / 2.;
				precision sTxn = - (dun_dx  +  dux_dn / t2) / 2.;
				precision sTyy = - duy_dy;
				precision sTyn = - (dun_dy  +  duy_dn / t2) / 2.;
				precision sTnn = - (dun_dn  +  ut / t) / t2;

				transverse_projection Xi(ut, ux, uy, un, zt, zn, t2);
				double_transverse_projection Xi_2(Xi, t2, t4);
				Xi_2.double_transverse_project_tensor(sTtt, sTtx, sTty, sTtn, sTxx, sTxy, sTxn, sTyy, sTyn, sTnn);

				precision sT_mag = sqrt(fabs(sTtt * sTtt  +  sTxx * sTxx  +  sTyy * sTyy  +  t4 * sTnn * sTnn  -  2. * (sTtx * sTtx  +  sTty * sTty  -  sTxy * sTxy  +  t2 * (sTtn * sTtn  -  sTxn * sTxn  -  sTyn * sTyn))));

				fprintf(Knpi, "%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, tau_pi * sT_mag);
			}
		}
	}
	fclose(Knpi);
#endif
}

void output_residual_shear_validity(const hydro_variables * const __restrict__ q, const fluid_velocity * const __restrict__ u, const precision * const e, double t, double dt, lattice_parameters lattice)
{
#ifdef ANISO_HYDRO
#ifdef PIMUNU
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision dx = lattice.lattice_spacing_x;
	precision dy = lattice.lattice_spacing_y;
	precision dz = lattice.lattice_spacing_eta;

	FILE *RpiInverse;
	FILE *piu_ortho0, *piu_ortho1, *piu_ortho2, *piu_ortho3;
	FILE *piz_ortho0, *piz_ortho1, *piz_ortho2, *piz_ortho3;
	FILE *pi_trace;

	char fname1[255];
	char fname2[255], fname3[255], fname4[255], fname5[255];
	char fname6[255], fname7[255], fname8[255], fname9[255];
	char fname10[255];

	sprintf(fname1, "output/RpiInv_%.3f.dat", t);
	sprintf(fname2, "output/piu_ortho0_%.3f.dat", t);
	sprintf(fname3, "output/piu_ortho1_%.3f.dat", t);
	sprintf(fname4, "output/piu_ortho2_%.3f.dat", t);
	sprintf(fname5, "output/piu_ortho3_%.3f.dat", t);
	sprintf(fname6, "output/piz_ortho0_%.3f.dat", t);
	sprintf(fname7, "output/piz_ortho1_%.3f.dat", t);
	sprintf(fname8, "output/piz_ortho2_%.3f.dat", t);
	sprintf(fname9, "output/piz_ortho3_%.3f.dat", t);
	sprintf(fname10, "output/pi_trace_%.3f.dat", t);

	RpiInverse = fopen(fname1, "w");
	piu_ortho0 = fopen(fname2, "w");
	piu_ortho1 = fopen(fname3, "w");
	piu_ortho2 = fopen(fname4, "w");
	piu_ortho3 = fopen(fname5, "w");
	piz_ortho0 = fopen(fname6, "w");
	piz_ortho1 = fopen(fname7, "w");
	piz_ortho2 = fopen(fname8, "w");
	piz_ortho3 = fopen(fname9, "w");
	pi_trace   = fopen(fname10, "w");

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

			#if (PT_MATCHING == 1)
				precision pt = q[s].pt;
			#else
				precision e_s = e[s];
				precision pl = q[s].pl;
				precision pt = (e_s - pl) / 2.;
			#endif

				precision pitt = q[s].pitt;
				precision pitx = q[s].pitx;
				precision pity = q[s].pity;
				precision pitn = q[s].pitn;
				precision pixx = q[s].pixx;
				precision pixy = q[s].pixy;
				precision pixn = q[s].pixn;
				precision piyy = q[s].piyy;
				precision piyn = q[s].piyn;
				precision pinn = q[s].pinn;

				precision ux = u[s].ux;
				precision uy = u[s].uy;
				precision un = u[s].un;
				precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t2 * un * un);

				precision utperp = sqrt(1.  +  ux * ux  +  uy * uy);
				precision zt = t * un / utperp;
				precision zn = ut / t / utperp;

				precision pi_mag = sqrt(fabs(pitt * pitt  +  pixx * pixx  +  piyy * piyy  +  t4 * pinn * pinn  -  2. * (pitx * pitx  +  pity * pity  -  pixy * pixy  +  t2 * (pitn * pitn  -  pixn * pixn  -  piyn * piyn))));

				precision piu0 = fabs(pitt * ut  -  pitx * ux  -  pity * uy  -  t2 * pitn * un);
				precision piu1 = fabs(pitx * ut  -  pixx * ux  -  pixy * uy  -  t2 * pixn * un);
				precision piu2 = fabs(pity * ut  -  pixy * ux  -  piyy * uy  -  t2 * piyn * un);
				precision piu3 = fabs(pitn * ut  -  pixn * ux  -  piyn * uy  -  t2 * pinn * un) * t;

				precision piz0 = fabs(zt * pitt  -  t2 * zn * pitn);
				precision piz1 = fabs(zt * pitx  -  t2 * zn * pixn);
				precision piz2 = fabs(zt * pity  -  t2 * zn * piyn);
				precision piz3 = fabs(zt * pitn  -  t2 * zn * pinn) * t;

				precision trpi = fabs(pitt  -  pixx  -  piyy  -  t2 * pinn);

				fprintf(RpiInverse, "%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, pi_mag / pt);
				fprintf(piu_ortho0, "%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, piu0 / (ut * pi_mag));
				fprintf(piu_ortho1, "%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, piu1 / (ut * pi_mag));
				fprintf(piu_ortho2, "%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, piu2 / (ut * pi_mag));
				fprintf(piu_ortho3, "%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, piu3 / (ut * pi_mag));
				fprintf(piz_ortho0, "%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, piz0 / (t * zn * pi_mag));
				fprintf(piz_ortho1, "%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, piz1 / (t * zn * pi_mag));
				fprintf(piz_ortho2, "%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, piz2 / (t * zn * pi_mag));
				fprintf(piz_ortho3, "%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, piz3 / (t * zn * pi_mag));
				fprintf(pi_trace,   "%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, trpi / pi_mag);
			}
		}
	}
	fclose(RpiInverse);
	fclose(piu_ortho0);
	fclose(piu_ortho1);
	fclose(piu_ortho2);
	fclose(piu_ortho3);
	fclose(piz_ortho0);
	fclose(piz_ortho1);
	fclose(piz_ortho2);
	fclose(piz_ortho3);
#endif
#endif
}


void output_aniso_bjorken(const hydro_variables * const __restrict__ q, const precision * const e, double e0, double t, lattice_parameters lattice)
{
#ifdef ANISO_HYDRO
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	FILE *energy, *plptratio, *temperature, *aLplot, *Lambdaplot;

	energy      = fopen("output/eratio.dat", "a");
	plptratio   = fopen("output/plptratio.dat", "a");
	temperature = fopen("output/T.dat", "a");
	aLplot      = fopen("output/aL.dat", "a");
	Lambdaplot	= fopen("output/Lambda.dat", "a");

	int s = central_index(lattice);

	precision e_s = e[s];
	precision pl = q[s].pl;
	precision pt = (e_s - pl) / 2.;

	precision T = effectiveTemperature(e_s) * hbarc;
	precision aL = compute_conformal_aL(pl, e_s);
	precision z = 1. / (aL * aL)  -  1.;
	precision t_200 = 1.;

	if(z > delta)
	{
		precision sqrtz = sqrt(z);
		precision t = atan(sqrtz) / sqrtz;
		t_200 = 1.  +  (1. + z) * t;
	}
	else if(z < -delta && z > -1.)
	{
		precision sqrtmz = sqrt(-z);
		precision t = atanh(sqrtmz) / sqrtmz;
		t_200 = 1.  +  (1. + z) * t;
	}
	else if(fabs(z) <= delta)
	{
		precision z2 = z  * z;
		precision z3 = z2 * z;
		precision z4 = z3 * z;
		precision z5 = z4 * z;
		precision z6 = z5 * z;
		t_200 = 2. + 0.6666666666666667*z - 0.1333333333333333*z2 + 0.05714285714285716*z3 - 0.031746031746031744*z4 + 0.020202020202020193*z5 - 0.013986013986013984*z6;
	}

	precision Lambda = pow(2. * e_s / (aL * aL * EOS_FACTOR * t_200), 0.25) * hbarc;

	fprintf(energy, 	"%.3f\t%.8f\n", t, e_s / e0);
	fprintf(plptratio, 	"%.3f\t%.8f\n", t, pl / pt);
	fprintf(temperature,"%.3f\t%.8f\n", t, T);
	fprintf(aLplot,		"%.3f\t%.8f\n", t, aL);
	fprintf(Lambdaplot,	"%.3f\t%.8f\n", t, Lambda);

	fclose(energy);
	fclose(plptratio);
	fclose(temperature);
	fclose(aLplot);
	fclose(Lambdaplot);
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

	FILE *energy, *plptratio;
	FILE *uxplot, *urplot;
	char fname1[255], fname2[255];
	char fname3[255], fname4[255];

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
				precision un = u[s].un;
				precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t2 * un * un);
				precision utperp = sqrt(1.  +  ux * ux  +  uy * uy);
				precision zt = t * un / utperp;
				precision zn = ut / t / utperp;
			#ifdef PIMUNU
				precision pitt = q[s].pitt;
				precision pitn = q[s].pitn;
				precision pinn = q[s].pinn;
			#else
				precision pitt = 0, pitn = 0, pinn = 0;
			#endif
				// pl = zmu.z.nu.Tmunu
				precision pl = e_s/3.  +  zt * zt * pitt  +  t4 * zn * zn * pinn  +  2. * t2 * zt * zn * pitn;
			#endif

				precision pt = (e_s - pl) / 2.;

				fprintf(energy, 	"%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, e_s);
				fprintf(plptratio,	"%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, pl / pt);
				fprintf(uxplot, 	"%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, ux);
				fprintf(urplot, 	"%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, ur);
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

	FILE *energy, *plptratio;
	FILE *uxplot, *uyplot, *urplot;
	char fname1[255], fname2[255];
	char fname3[255], fname4[255], fname5[255];

	sprintf(fname1, "output/e_%.3f.dat", t);
	sprintf(fname2, "output/plpt_%.3f.dat", t);
	sprintf(fname3, "output/ux_%.3f.dat", t);
	sprintf(fname4, "output/uy_%.3f.dat", t);
	sprintf(fname5, "output/ur_%.3f.dat", t);

	energy      = fopen(fname1, "w");
	plptratio 	= fopen(fname2, "w");
	uxplot    	= fopen(fname3, "w");
	uyplot    	= fopen(fname4, "w");
	urplot		= fopen(fname5, "w");

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
			#if (PT_MATCHING == 1)
				precision pt = q[s].pt;
			#else
				precision pt = (e_s - pl) / 2.;
			#endif
			#else
				precision un = u[s].un;
				precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t2 * un * un);
				precision utperp = sqrt(1.  +  ux * ux  +  uy * uy);
				precision zt = t * un / utperp;
				precision zn = ut / t / utperp;

				precision p = equilibriumPressure(e_s);
			#ifdef PIMUNU
				precision pitt = q[s].pitt;
				precision pitn = q[s].pitn;
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
				// pt = -Ximunu.Tmunu / 2 (assumes pimunu is traceless and orthogonal to u)
				precision pt = (3. * (p + Pi)  -  pl) / 2.;
			#endif

				fprintf(energy, 	"%.3f\t%.3f\t%.3f\t%.8e\n", x, y, z, e_s);
				fprintf(plptratio,	"%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, pl / pt);
				fprintf(uxplot, 	"%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, ux);
				fprintf(uyplot, 	"%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, uy);
				fprintf(urplot, 	"%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, ur);
			}
		}
	}
	fclose(energy);
	fclose(plptratio);
	fclose(uxplot);
	fclose(uyplot);
	fclose(urplot);
}


void output_dynamical_variables(double t, double dt, lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro)
{
	int initial_type = initial.initialConditionType;

	if(initial_type == 1)
	{
		precision T0 = initial.initialCentralTemperatureGeV;
		precision e0 = equilibriumEnergyDensity(T0 / hbarc);

		output_aniso_bjorken(q, e, e0, t, lattice);
	}
	else if(initial_type == 2 || initial_type == 3)
	{
		output_gubser(q, u, e, t, lattice);
		output_residual_gradients(q, u, up, e, t, dt, lattice, hydro);
		output_residual_shear_validity(q, u, e, t, dt, lattice);
	}
	else if(initial_type == 4)
	{
		output_optical_glauber(q, u, e, t, lattice);
		output_residual_gradients(q, u, up, e, t, dt, lattice, hydro);
		output_residual_shear_validity(q, u, e, t, dt, lattice);
	}
}


void output_semi_analytic_solution_if_any(lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro)
{
	switch(initial.initialConditionType)
	{
		case 1:	// anisotropic Bjorken
		{
			printf("\nRunning semi-analytic anisotropic Bjorken solution...\n");
			run_semi_analytic_aniso_bjorken(lattice, initial, hydro);
			break;
		}
		default:
		{
			break;
		}
	}
}













