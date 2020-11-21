#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string>
#include <string.h>
#include <iostream>
#include <iomanip>
using namespace std;
#include <sstream>
#include <fstream>



int main()
{
	int nx;
	int ny;
	int nz;

    char fname1[255];
	char fname2[255];
	char fname3[255];
	char fname4[255];
	char fname5[255];
	char fname6[255];

    double times[8] = {0.01, 0.51, 1.01, 2.01, 3.01, 4.01, 5.01, 6.01};	// hard coded (assumes output_interval = 0.5)

    double x_shift;
    bool got_shift = false;

    // conformal viscous hydro
    for(int n = 0; n < 8; n++)
    {
    	double t = times[n];

    	printf("Scanning t = %.3f files\n", t);

    	// open 3d files
    	FILE * energy_3d;
    	FILE * plpt_3d;
    	FILE * piperp_3d;
    	FILE * Wperp_3d;
    	FILE * ux_3d;
    	FILE * un_3d;

    	sprintf(fname1, "../../../../output/e_%.3f.dat", t);
		sprintf(fname2, "../../../../output/plpt_%.3f.dat", t);
		sprintf(fname3, "../../../../output/piperp_pt_%.3f.dat", t);
		sprintf(fname4, "../../../../output/WTz_pl_2pt_%.3f.dat", t);
		sprintf(fname5, "../../../../output/ux_%.3f.dat", t);
		sprintf(fname6, "../../../../output/un_%.3f.dat", t);

		energy_3d = fopen(fname1, "r");
		plpt_3d   = fopen(fname2, "r");
		piperp_3d = fopen(fname3, "r");
		Wperp_3d  = fopen(fname4, "r");
		ux_3d     = fopen(fname5, "r");
		un_3d     = fopen(fname6, "r");

    	fscanf(energy_3d, "%d\n%d\n%d", &nx, &ny, &nz);
    	fscanf(plpt_3d,   "%d\n%d\n%d", &nx, &ny, &nz);
    	fscanf(piperp_3d, "%d\n%d\n%d", &nx, &ny, &nz);
    	fscanf(Wperp_3d,  "%d\n%d\n%d", &nx, &ny, &nz);
    	fscanf(ux_3d,     "%d\n%d\n%d", &nx, &ny, &nz);
    	fscanf(un_3d,     "%d\n%d\n%d", &nx, &ny, &nz);

    	// x-axis
    	FILE * energy_x;
    	FILE * plpt_x;
    	FILE * piperp_x;
    	FILE * Wperp_x;				// may need to do something special for WTz (expect to be zero along both x and eta axis)
    	FILE * ux_x;
    	FILE * un_x;

    	sprintf(fname1, "x_axis/e_%.3f.dat", t);
		sprintf(fname2, "x_axis/plpt_%.3f.dat", t);
		sprintf(fname3, "x_axis/piperp_pt_%.3f.dat", t);
		sprintf(fname4, "x_axis/WTz_pl_2pt_%.3f.dat", t);
		sprintf(fname5, "x_axis/ux_%.3f.dat", t);
		sprintf(fname6, "x_axis/un_%.3f.dat", t);

		energy_x = fopen(fname1, "w");
		plpt_x   = fopen(fname2, "w");
		piperp_x = fopen(fname3, "w");
		Wperp_x  = fopen(fname4, "w");
		ux_x     = fopen(fname5, "w");
		un_x     = fopen(fname6, "w");

		// z-azis
    	FILE * energy_z;
    	FILE * plpt_z;
    	FILE * piperp_z;
    	FILE * Wperp_z;				// x shifted from 0 for Wperp plot along eta-axis;
    	FILE * ux_z;
    	FILE * un_z;

    	sprintf(fname1, "z_axis/e_%.3f.dat", t);
		sprintf(fname2, "z_axis/plpt_%.3f.dat", t);
		sprintf(fname3, "z_axis/piperp_pt_%.3f.dat", t);
		sprintf(fname4, "z_axis/WTz_pl_2pt_%.3f.dat", t);
		sprintf(fname5, "z_axis/ux_%.3f.dat", t);
		sprintf(fname6, "z_axis/un_%.3f.dat", t);

		energy_z = fopen(fname1, "w");
		plpt_z   = fopen(fname2, "w");
		piperp_z = fopen(fname3, "w");
		Wperp_z  = fopen(fname4, "w");
		ux_z     = fopen(fname5, "w");
		un_z     = fopen(fname6, "w");

    	for(int k = 0; k < nz; k++)
    	{
    		for(int j = 0; j < ny; j++)
    		{
    			for(int i = 0; i < nx; i++)
    			{
    				double x, y, z;

    				double energy, plpt, piperp, Wperp, ux, un;

    				fscanf(energy_3d, "%lf\t%lf\t%lf\t%lf", &x, &y, &z, &energy);
    				fscanf(plpt_3d,   "%lf\t%lf\t%lf\t%lf", &x, &y, &z, &plpt);
    				fscanf(piperp_3d, "%lf\t%lf\t%lf\t%lf", &x, &y, &z, &piperp);
    				fscanf(Wperp_3d,  "%lf\t%lf\t%lf\t%lf", &x, &y, &z, &Wperp);
    				fscanf(ux_3d,     "%lf\t%lf\t%lf\t%lf", &x, &y, &z, &ux);
    				fscanf(un_3d,     "%lf\t%lf\t%lf\t%lf", &x, &y, &z, &un);

    				if(j == (ny - 1)/2 && k == (nz - 1)/2)
    				{
    					fprintf(energy_x, "%.2f\t%.4e\n", x, energy);
    					fprintf(plpt_x,   "%.2f\t%.4e\n", x, plpt);
    					fprintf(piperp_x, "%.2f\t%.4e\n", x, piperp);
    					fprintf(Wperp_x,  "%.2f\t%.4e\n", x, Wperp);
    					fprintf(ux_x,     "%.2f\t%.4e\n", x, ux);
    					fprintf(un_x,     "%.2f\t%.4e\n", x, un);
       				}

       				if(i == (nx - 1)/2 && j == (ny - 1)/2)
    				{
    					fprintf(energy_z, "%.2f\t%.4e\n", z, energy);
    					fprintf(plpt_z,   "%.2f\t%.4e\n", z, plpt);
    					fprintf(piperp_z, "%.2f\t%.4e\n", z, piperp);
    					fprintf(ux_z,     "%.2f\t%.4e\n", z, ux);
    					fprintf(un_z,     "%.2f\t%.4e\n", z, un);
       				}

                    // offset x to get nonzero Wperp
                    if(i == (int)ceil(0.75*(nx - 1)) && j == (ny - 1)/2)
                    {
                        if(!got_shift)
                        {
                            x_shift = x;
                            got_shift = true;
                        }

                        fprintf(Wperp_z,  "%.2f\t%.4e\n", z, Wperp);
                    }
    			}
    		}
    	}

    	fclose(energy_x);
    	fclose(plpt_x);
    	fclose(piperp_x);
    	fclose(Wperp_x);
    	fclose(ux_x);
    	fclose(un_x);

    	fclose(energy_z);
    	fclose(plpt_z);
    	fclose(piperp_z);
    	fclose(Wperp_z);
    	fclose(ux_z);
    	fclose(un_z);

    	fclose(energy_3d);
    	fclose(plpt_3d);
    	fclose(piperp_3d);
    	fclose(Wperp_3d);
    	fclose(ux_3d);
    	fclose(un_3d);

    }

    printf("\nShifted x = %.2f for Wperp plot\n\n", x_shift);

}

