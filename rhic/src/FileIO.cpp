
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../include/FileIO.h"
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"
#include "../include/EquationofState.h"
#include "../include/Hypergeometric.h"

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



void output(const precision * const var, double t, const char * var_name, int nx, int ny, int nz, double dx, double dy, double dz)
{
	FILE * output;
	char fname[255];
	sprintf(fname, "output/%s_%.3f.dat", var_name, t);

	output = fopen(fname, "w");

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

				fprintf(output, "%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, var[s]);
			}
		}
	}
	fclose(output);
}

void output_ur(const precision * const ux, const precision * const uy, double t, int nx, int ny, int nz, double dx, double dy, double dz)
{
	FILE * output;
	char fname[255];
	sprintf(fname, "output/ur_%.3f.dat", t);

	output = fopen(fname, "w");

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

				double ux_s = ux[s];
				double uy_s = uy[s];

				fprintf(output, "%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, sqrt(ux_s*ux_s + uy_s*uy_s));
			}
		}
	}
	fclose(output);
}


void output_aL(const precision * const e, const precision * const pl, double t, int nx, int ny, int nz, double dx, double dy, double dz)
{
	FILE * output;
	char fname[255];
	sprintf(fname, "output/aL_%.3f.dat", t);

	output = fopen(fname, "w");

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

				precision aL = compute_conformal_aL(pl[s], e[s]);

				fprintf(output, "%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, aL);
			}
		}
	}
	fclose(output);
}


void output_gubser_test(const precision * const e, const precision * const pl, const precision * const ux, const precision * const uy, double t, int nx, int ny, int nz, double dx, double dy, double dz)
{
	FILE * energy;
	FILE * plptratio;
	FILE * temperature;
	FILE * aLplot;
	FILE * Lambdaplot;
	FILE * uxplot;
	FILE * uyplot;
	FILE * urplot;

	char fname1[255];
	char fname2[255];
	char fname3[255];
	char fname4[255];
	char fname5[255];
	char fname6[255];
	char fname7[255];
	char fname8[255];

	sprintf(fname1, "output/e_%.3f.dat", t);
	sprintf(fname2, "output/plpt_%.3f.dat", t);
	sprintf(fname3, "output/T_%.3f.dat", t);
	sprintf(fname4, "output/aL_%.3f.dat", t);
	sprintf(fname5, "output/Lambda_%.3f.dat", t);
	sprintf(fname6, "output/ux_%.3f.dat", t);
	sprintf(fname7, "output/uy_%.3f.dat", t);
	sprintf(fname8, "output/ur_%.3f.dat", t);

	energy      = fopen(fname1, "w");
	plptratio   = fopen(fname2, "w");
	temperature = fopen(fname3, "w");
	aLplot      = fopen(fname4, "w");
	Lambdaplot  = fopen(fname5, "w");
	uxplot      = fopen(fname6, "w");
	uyplot      = fopen(fname7, "w");
	urplot      = fopen(fname8, "w");

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

				precision ux_s = ux[s];
				precision uy_s = uy[s];
				precision ur = sqrt(ux_s * ux_s  +  uy_s * uy_s);

				precision e_s = e[s];
				precision pl_s = pl[s];
				precision pt_s = 0.5 * (e_s - pl_s);

				precision T = effectiveTemperature(e_s) * hbarc;
				precision aL = compute_conformal_aL(pl_s, e_s);

				precision xi = 1.0 / (aL * aL)  -  1.0;

				precision hype = 1.0;

				if(xi > delta)
				{
					hype = atan(sqrt(xi)) /sqrt(xi);
				}
				else if(xi < - delta && xi > -1.0)
				{
					hype = atanh(sqrt(-xi)) / sqrt(-xi);
				}

				precision Lambda = pow(2.0 * e_s / (aL * aL * EOS_FACTOR * t_200(xi, hype)), 0.25) * hbarc;

				fprintf(energy, "%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, e_s);
				fprintf(plptratio, "%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, pl_s / pt_s);
				fprintf(temperature, "%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, T);
				fprintf(aLplot, "%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, aL);
				fprintf(Lambdaplot, "%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, Lambda);
				fprintf(uxplot, "%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, ux_s);
				fprintf(uyplot, "%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, uy_s);
				fprintf(urplot, "%.3f\t%.3f\t%.3f\t%.8f\n", x, y, z, ur);
			}
		}
	}
	fclose(energy);
	fclose(plptratio);
	fclose(temperature);
	fclose(aLplot);
	fclose(Lambdaplot);
	fclose(uxplot);
	fclose(uyplot);
	fclose(urplot);
}


void output_bjorken_test(const precision * const e, const precision * const pl, double e0, double t, int nx, int ny, int nz)
{
	FILE * energy;
	FILE * plptratio;
	FILE * temperature;
	FILE * aLplot;
	FILE * Lambdaplot;

	energy      = fopen("output/eratio.dat", "a");
	plptratio   = fopen("output/plptratio.dat", "a");
	temperature = fopen("output/T.dat", "a");
	aLplot      = fopen("output/aL.dat", "a");
	Lambdaplot	= fopen("output/Lambda.dat", "a");

	int s = center_index(nx, ny, nz, nx + 4, ny + 4, nz + 4);

	precision e_s = e[s];
	precision pl_s = pl[s];
	precision pt_s = 0.5 * (e_s - pl_s);

	precision T = effectiveTemperature(e_s) * hbarc;
	precision aL = compute_conformal_aL(pl_s, e_s);

	precision z = 1.0 / (aL * aL)  -  1.0;

	precision hype = 1.0;

	if(z > delta) hype = atan(sqrt(z)) /sqrt(z);
	else if(z < - delta && z > -1.0) hype = atanh(sqrt(-z)) / sqrt(-z);

	precision Lambda = pow(2.0 * e_s / (aL * aL * EOS_FACTOR * t_200(z, hype)), 0.25) * hbarc;

	fprintf(energy, "%.3f\t%.8f\n", t, e_s);
	fprintf(plptratio, "%.3f\t%.8f\n", t, pl_s / pt_s);
	fprintf(temperature, "%.3f\t%.8f\n", t, T);
	fprintf(aLplot, "%.3f\t%.8f\n", t, aL);
	fprintf(Lambdaplot, "%.3f\t%.8f\n", t, Lambda);

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
		output_bjorken_test(e, q->pl, e0, t, nx, ny, nz);
		return;
	}
	if(initialConditionType == 2 || initialConditionType == 3)
	{
		output_gubser_test(e, q->pl, u->ux, u->uy, t, nx, ny, nz, dx, dy, dz);
		return;
	}

	output(e, t, "e", nx, ny, nz, dx, dy, dz);

	output(u->ut, t, "ut", nx, ny, nz, dx, dy, dz);
	output(u->ux, t, "ux", nx, ny, nz, dx, dy, dz);
	output(u->uy, t, "uy", nx, ny, nz, dx, dy, dz);
	output(u->un, t, "un", nx, ny, nz, dx, dy, dz);

	output_ur(u->ux, u->uy, t, nx, ny, nz, dx, dy, dz);
	output_aL(e, q->pl, t, nx, ny, nz, dx, dy, dz);

	output(q->pl, t, "pl", nx, ny, nz, dx, dy, dz);

#if (PT_MATCHING == 1)
	output(q->pt, t, "pt", nx, ny, nz, dx, dy, dz);
#endif

#ifdef PIMUNU
	output(q->pitt, t, "pitt", nx, ny, nz, dx, dy, dz);
	output(q->pitx, t, "pitx", nx, ny, nz, dx, dy, dz);
	output(q->pity, t, "pity", nx, ny, nz, dx, dy, dz);
	output(q->pitn, t, "pitn", nx, ny, nz, dx, dy, dz);
	output(q->pixx, t, "pixx", nx, ny, nz, dx, dy, dz);
	output(q->pixy, t, "pixy", nx, ny, nz, dx, dy, dz);
	output(q->pixn, t, "pixn", nx, ny, nz, dx, dy, dz);
	output(q->piyy, t, "piyy", nx, ny, nz, dx, dy, dz);
	output(q->piyn, t, "piyn", nx, ny, nz, dx, dy, dz);
	output(q->pinn, t, "pinn", nx, ny, nz, dx, dy, dz);
#endif
#ifdef WTZMU
	output(q->WtTz, t, "WtTz", nx, ny, nz, dx, dy, dz);
	output(q->WxTz, t, "WxTz", nx, ny, nz, dx, dy, dz);
	output(q->WyTz, t, "WyTz", nx, ny, nz, dx, dy, dz);
	output(q->WnTz, t, "WnTz", nx, ny, nz, dx, dy, dz);
#endif
}


