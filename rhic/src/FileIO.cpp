#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../include/FileIO.h"
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"
#include "../include/EquationofState.h"
#include "../include/TransportCoefficients.h"

const double hbarc = 0.197326938;


inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}

int center_index(int nx, int ny, int nz, int ncx, int ncy, int ncz)
{
	// (not sure what this means)
	int ictr = (nx % 2 == 0) ? ncx/2 : (ncx-1)/2;
	int jctr = (ny % 2 == 0) ? ncy/2 : (ncy-1)/2;
	int kctr = (nz % 2 == 0) ? ncz/2 : (ncz-1)/2;

	return ictr  +  ncx * (jctr  +  ncy * kctr);
}



double compute_conformal_aL(double pl, double e)
{
	precision x   = pl / e;				// x = pl / e
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


// void output_ur(const precision * const ux, const precision * const uy, double t, int nx, int ny, int nz, double dx, double dy, double dz)
// {
// 	FILE * output;
// 	char fname[255];
// 	sprintf(fname, "output/ur_%.3f.dat", t);

// 	output = fopen(fname, "w");

// 	for(int k = 2; k < nz + 2; k++)
// 	{
// 		double z = (k - 2.0 - (nz - 1.0)/2.0) * dz;

// 		for(int j = 2; j < ny + 2; j++)
// 		{
// 			double y = (j - 2.0 - (ny - 1.0)/2.0) * dy;

// 			for(int i = 2; i < nx + 2; i++)
// 			{
// 				double x = (i - 2.0 - (nx - 1.0)/2.0) * dx;

// 				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

// 				double ux_s = ux[s];
// 				double uy_s = uy[s];

// 				fprintf(output, "%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, sqrt(ux_s*ux_s + uy_s*uy_s));
// 			}
// 		}
// 	}
// 	fclose(output);
// }


void output_shear_validity(const conserved_variables * const __restrict__ q, const fluid_velocity * const __restrict__ u, const precision * const e, double t, int nx, int ny, int nz, double dx, double dy, double dz)
{
	FILE *RpiInverse;
	FILE *piu_ortho0, *piu_ortho1, *piu_ortho2, *piu_ortho3;	
	FILE *piz_ortho0, *piz_ortho1, *piz_ortho2, *piz_ortho3;	// I forgot tracelessness
	FILE *pi_trace;

	char fname1[255];
	char fname2[255], fname3[255], fname4[255], fname5[255];
	char fname6[255], fname7[255], fname8[255], fname9[255];
	char fname10[255];

	sprintf(fname1, "output/RpiInv1_%.3f.dat", t);
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
		double z = (k - 2.0 - (nz - 1.0)/2.0) * dz;

		for(int j = 2; j < ny + 2; j++)
		{
			double y = (j - 2.0 - (ny - 1.0)/2.0) * dy;

			for(int i = 2; i < nx + 2; i++)
			{
				double x = (i - 2.0 - (nx - 1.0)/2.0) * dx;

				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

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

			#if (PT_MATCHING == 1)
				precision pt = q[s].pt;
			#else
				precision e_s = e[s];
				precision pl = q[s].pl;
				precision pt = (e_s - pl) / 2.;
			#endif

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
}


void output_gubser_test(const conserved_variables * const __restrict__ q, const fluid_velocity  * const __restrict__ u, const precision * const e, double t, int nx, int ny, int nz, double dx, double dy, double dz)
{
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
				precision pl = q[s].pl;
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


void output_bjorken_test(const conserved_variables * const __restrict__ q, const precision * const e, double e0, double t, int nx, int ny, int nz)
{
	FILE *energy, *plptratio, *temperature, *aLplot, *Lambdaplot;

	energy      = fopen("output/eratio.dat", "a");
	plptratio   = fopen("output/plptratio.dat", "a");
	temperature = fopen("output/T.dat", "a");
	aLplot      = fopen("output/aL.dat", "a");
	Lambdaplot	= fopen("output/Lambda.dat", "a");

	int s = center_index(nx, ny, nz, nx + 4, ny + 4, nz + 4);

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
}


void output_dynamical_variables(double t, int nx, int ny, int nz, double dx, double dy, double dz, int initialConditionType, double e0)
{
	if(initialConditionType == 1)
	{
		output_bjorken_test(q, e, e0, t, nx, ny, nz);
		return;
	}
	if(initialConditionType == 2 || initialConditionType == 3)
	{
		output_gubser_test(q, u, e, t, nx, ny, nz, dx, dy, dz);

	#ifdef PIMUNU
		output_shear_validity(q, u, e, t, nx, ny, nz, dx, dy, dz);
	#endif

		return;
	}
}



