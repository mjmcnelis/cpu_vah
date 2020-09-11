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
    double times[4] = {0.05, 2.05, 5.05, 10.05};		// selected times to plot

    for(int n = 0; n < 4; n++)
    {
    	double t = times[n];

    	printf("Viscous anisotropic hydro: scanning t = %.3f files\n", t);

    	// open 3d energy density file
    	FILE * energy_3d;
    	sprintf(fname1, "../../../../output/e_%.3f.dat", t);
		energy_3d = fopen(fname1, "r");

        // get grid points
    	fscanf(energy_3d, "%d\n%d\n%d", &nx, &ny, &nz);

    	// extract energy density in transverse plane (eta_s = 0)
    	FILE * energy_xy;
    	sprintf(fname1, "xy_plane/e_%.3f.dat", t);
		energy_xy = fopen(fname1, "w");

    	for(int k = 0; k < nz; k++)
    	{
    		for(int j = 0; j < ny; j++)
    		{
    			for(int i = 0; i < nx; i++)
    			{
    				double x;
                    double y;
                    double z;
                    double energy;

    				fscanf(energy_3d, "%lf\t%lf\t%lf\t%lf", &x, &y, &z, &energy);

    				if(k == (nz - 1)/2)
    				{
    					fprintf(energy_xy, "%.2f\t%.2f\t%.4e\n", x, y, energy);
       				}
    			}
    		}
    	}

    	fclose(energy_xy);
    }
}


