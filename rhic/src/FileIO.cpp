/*
 * FileIO.c
 *
 *  Created on: Oct 24, 2015
 *      Author: bazow
 */
#include <stdlib.h>
#include <stdio.h>
#include "../include/FileIO.h"
#include "../include/Parameters.h"
#include "../include/DynamicalVariables.h"


inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


void output(const PRECISION * const var, double t, const char * var_name, int nx, int ny, int nz, double dx, double dy, double dz)
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


void output_dynamical_variables(double t, int nx, int ny, int nz, double dx, double dy, double dz)
{
	output(e, t, "e", nx, ny, nz, dx, dy, dz);

	output(u->ut, t, "ut", nx, ny, nz, dx, dy, dz);
	output(u->ux, t, "ux", nx, ny, nz, dx, dy, dz);
	output(u->uy, t, "uy", nx, ny, nz, dx, dy, dz);
	output(u->un, t, "un", nx, ny, nz, dx, dy, dz);

	output(q->pl, t, "pl", nx, ny, nz, dx, dy, dz);
	output(q->pt, t, "pt", nx, ny, nz, dx, dy, dz);

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


