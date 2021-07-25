#include <stdlib.h>
#include <stdio.h>
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"
#include "../include/Parameters.h"
#include "../include/RungeKutta.h"

hydro_variables *q, *q1, *q2;			// dynamical variables
fluid_velocity *u, *up;					// fluid velocity
precision *e;                           // energy density
precision *lambda, *aT, *aL;            // anisotropic variables

hydro_variables *q3;					// additional variables for third-order Runge-Kutta
hydro_variables *q4;
fluid_velocity *uI;


#ifdef ANISO_HYDRO
#ifdef LATTICE_QCD
	int *aniso_regulation;				// for monitoring regulation of (lambda, aT, aL) in aniso hydro
#endif
#endif

#ifdef MONITOR_REGULATIONS
	int *viscous_regulation;			// for monitoring regulation of (piperp, Wperp) or (pimunu, Pi)
#endif

#ifdef MONITOR_PLPT
	int *plpt_regulation;				// for monitoring regulation of pl, pt in aniso hydro
#endif

#ifdef MONITOR_B
	int *b_regulation;					// for monitoring regulation of b in aniso hydro
#endif

#ifdef MONITOR_TTAUMU
	float *Tmunu_violations;			// for monitoring violations in T^{\tau\mu}
#endif

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

	aniso_regulation = (int *)calloc(length, sizeof(int));
#endif
#endif

#ifdef MONITOR_REGULATIONS
	viscous_regulation = (int *)calloc(length, sizeof(int));
#endif

#ifdef MONITOR_PLPT
	plpt_regulation = (int *)calloc(length, sizeof(int));
#endif

#ifdef MONITOR_B
	b_regulation = (int *)calloc(length, sizeof(int));
#endif

#ifdef MONITOR_TTAUMU
	Tmunu_violations = (float *)calloc(length, sizeof(float));
#endif

	u  = (fluid_velocity *)calloc(length, sizeof(fluid_velocity));
	up = (fluid_velocity *)calloc(length, sizeof(fluid_velocity));

	q  = (hydro_variables *)calloc(length, sizeof(hydro_variables));
	q1 = (hydro_variables *)calloc(length, sizeof(hydro_variables));
	q2 = (hydro_variables *)calloc(length, sizeof(hydro_variables));

	// allocate additional variables for third-order Runge-Kutta
#if (STAGES > 2)
	q3 = (hydro_variables *)calloc(length, sizeof(hydro_variables));
#endif
#if (STAGES > 3)
	q4 = (hydro_variables *)calloc(length, sizeof(hydro_variables));
#endif

#if (STAGES > 2)
	uI = (fluid_velocity *)calloc(length, sizeof(fluid_velocity));
#endif
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
	free(q);
	free(q1);
	free(q2);

	free(q3);
	free(q4);
	free(uI);
}



