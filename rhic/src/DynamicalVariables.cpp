#include <stdlib.h>
#include <stdio.h>
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"
#include "../include/Parameters.h"


hydro_variables *q, *Q, *qI;
fluid_velocity *u, *up, *uI;
precision *e, *lambda, *aT, *aL;

// need anisotropic variables


void allocate_memory(lattice_parameters lattice)
{
	// allocate memory for the dynamical variables

	int nx = lattice.lattice_points_x;				// physical grid points
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	int ncx = nx + 4;								// computational grid points = physical + ghost + white
	int ncy = ny + 4;
	int ncz = nz + 4;

	int length = ncx * ncy * ncz;					// size of each array

	size_t bytes = sizeof(precision);

	e = (precision *)calloc(length, bytes);

#ifdef ANISO_HYDRO 									// anisotropic variables 
#ifdef LATTICE_QCD
	lambda = (precision *)calloc(length, bytes);
	aT = (precision *)calloc(length, bytes);
	aL = (precision *)calloc(length, bytes);
#endif
#endif

	u  = (fluid_velocity *)calloc(length, sizeof(fluid_velocity));
	up = (fluid_velocity *)calloc(length, sizeof(fluid_velocity));
	uI = (fluid_velocity *)calloc(length, sizeof(fluid_velocity));

	q  = (hydro_variables *)calloc(length, sizeof(hydro_variables));
	Q  = (hydro_variables *)calloc(length, sizeof(hydro_variables));
	qI = (hydro_variables *)calloc(length, sizeof(hydro_variables));
}


void swap_hydro_variables(hydro_variables ** hydro_1, hydro_variables ** hydro_2)
{
	hydro_variables * hydro_temp = *hydro_1;		// swap hydro_1 <-> hydro_2
	*hydro_1 = *hydro_2;
	*hydro_2 = hydro_temp;
}


void swap_fluid_velocity(fluid_velocity ** velocity_1, fluid_velocity ** velocity_2)
{
	fluid_velocity * velocity_temp = *velocity_1;	// swap velocity_1 <-> velocity_2
	*velocity_1 = *velocity_2;
	*velocity_2 = velocity_temp;
}


void free_memory()
{
	free(e);
	free(u);
	free(up);
	free(uI);
	free(q);
	free(qI);
	free(Q);
}



