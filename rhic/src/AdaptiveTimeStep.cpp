#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../include/Precision.h"
#include "../include/FluxTerms.h"
#include "../include/DynamicalVariables.h"
#include "../include/SourceTerms.h"
#include "../include/NeighborCells.h"
#include "../include/EquationOfState.h"
#include "../include/TransportCoefficients.h"

const precision fraction = 0.1;

// set the adaptive time step
precision set_time_step(hydro_time_scales dt_hydro, precision dt_min)
{
	precision dt_CFL 	= dt_hydro.dt_CFL;
	precision dt_micro 	= dt_hydro.dt_micro;
	precision dt_rate	= dt_hydro.dt_rate;

	printf("%lf\t%lf\t%lf\n", dt_CFL, dt_micro, dt_rate);

	precision dt = fmin(dt_CFL, fraction * fmin(dt_micro, dt_rate));

	if(dt < dt_min)
	{
		printf("set_time_step error: dt = %.6f < %.3f\n", dt, dt_min);
	}

	return dt_min * fmax(1., floor(dt / dt_min));	// round down dt to numerical precision
}


inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}

inline precision central_derivative(const precision * const __restrict__ f, int n, precision dx)
{
	return (f[n + 1] - f[n]) / (2. * dx);		// f[n] = fm  |	 f[n+1] = fp
}


precision compute_CFL_time_scale(const precision * const __restrict__ vxi, const precision * const __restrict__ vyj, const precision * const __restrict__ vnk, precision ut, precision ux, precision uy, precision un, precision dx, precision dy, precision dn, precision Theta)
{
	precision vx = ux / ut;
	precision vy = uy / ut;
	precision vn = un / ut;

	precision ax = compute_max_local_propagation_speed(vxi, vx, Theta);
	precision ay = compute_max_local_propagation_speed(vyj, vy, Theta);
	precision an = compute_max_local_propagation_speed(vnk, vn, Theta);

	return 0.125 / fmax(fabs(ax / dx), fmax(fabs(ay / dy), fabs(an / dn)));
	//return  0.125 / fmax(fabs(vx / dx), fmax(fabs(vy / dy), fabs(vn / dn)));
}



hydro_time_scales compute_hydro_time_scales(precision t, const hydro_variables * const __restrict__ q, const precision * const __restrict__ e, const fluid_velocity * const __restrict__ u, const fluid_velocity * const __restrict__ up, int nx, int ny, int nz, precision dt, precision dt_prev, precision dx, precision dy, precision dn, precision etabar_const, precision dt_min, hydro_parameters hydro)
{
	hydro_time_scales dt_hydro;						// default value for hydro time scales
	dt_hydro.dt_CFL = 100;
	dt_hydro.dt_micro = 100;
	dt_hydro.dt_rate = 100;

	precision Theta = hydro.flux_limiter;

	precision t2 = t * t;
	int stride_y = nx + 4;							// strides for neighbor cells along x, y, n (stride_x = 1)
	int stride_z = (nx + 4) * (ny + 4);				// stride formulas based from linear_column_index()

	precision ui1[6];								// fluid velocity of neighbor cells along x [i-1, i+1]
	precision uj1[6];								// fluid velocity of neighbor cells along y [j-1, j+1]
	precision uk1[6];								// fluid velocity of neighbor cells along n [k-1, k+1]

	precision vxi[4];								// vx of neighbor cells along x [i-2, i-1, i+1, i+2]
	precision vyj[4];								// vy of neighbor cells along y [j-2, j-1, j+1, j+2]
	precision vnk[4];								// vn of neighbor cells along n [k-2, k-1, k+1, k+2]

	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

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

				precision e_s = e[s];

				precision ux = u[s].ux;		// current fluid velocity
				precision uy = u[s].uy;
				precision un = u[s].un;
				precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t2 * un * un);

				precision ux_p = up[s].ux;	// previous fluid velocity
				precision uy_p = up[s].uy;
				precision un_p = up[s].un;

				get_fluid_velocity_neighbor_cells(u[simm], u[sim], u[sip], u[sipp], u[sjmm], u[sjm], u[sjp], u[sjpp], u[skmm], u[skm], u[skp], u[skpp], ui1, uj1, uk1, vxi, vyj, vnk, t2);


				// compute the hydrodynamic time scales for fluid cell s
				precision dt_CFL = compute_CFL_time_scale(vxi, vyj, vnk, ut, ux, uy, un, dx, dy, dn, Theta);

				// update the minimum hydro time scales of the entire grid
				dt_hydro.dt_CFL	= fmin(dt_CFL, dt_hydro.dt_CFL);
			}
		}
	}
	return dt_hydro;
}






