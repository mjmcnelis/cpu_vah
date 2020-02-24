#include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h> 
#include "../include/Macros.h"
#include "../include/Precision.h"
#include "../include/Parameters.h"
#include "../include/FreezeoutFinder.h"


freezeout_finder::freezeout_finder(lattice_parameters lattice)
{
	nx = lattice.lattice_points_x;
	ny = lattice.lattice_points_y;
	nz = lattice.lattice_points_eta;

	dt = lattice.fixed_time_step;
	dx = lattice.lattice_spacing_x;
	dy = lattice.lattice_spacing_y;
	dz = lattice.lattice_spacing_eta;
	tau_coarse_factor = lattice.tau_coarse_factor;

	// initialize cornelius variables for freezeout surface finding (see example_4d() in example_cornelius)
#ifdef BOOST_INVARIANT
	dimension = 3;
	lattice_spacing = new double[dimension];
	lattice_spacing[0] = dt * tau_coarse_factor;  
    lattice_spacing[1] = dx;                              
    lattice_spacing[2] = dy;

	if(!(nx > 1 && ny > 1))
	{
		printf("freezeout_finder::freezeout_finder error: 2d spatial grid needs to have finite size\n");
		exit(-1);
 	}

 	//printf("%lf\t%lf\t%lf\n", lattice_spacing[0], lattice_spacing[1], lattice_spacing[2]);

#else
	dimension = 4;
	lattice_spacing = new double[dimension];
	lattice_spacing[0] = dt * tau_coarse_factor;  
    lattice_spacing[1] = dx;                              
    lattice_spacing[2] = dy;
    lattice_spacing[3] = dz;

    if(!(nx > 1 && ny > 1 && nz > 1))
	{
		printf("freezeout_finder::freezeout_finder error: 3d spatial grid needs to have finite size\n");
		exit(-1);
 	}
#endif
}

freezeout_finder::~freezeout_finder()
{

}


