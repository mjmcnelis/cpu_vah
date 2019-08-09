#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"
#include "../include/SourceTerms.h"
#include "../include/NeighborCells.h"

const precision fraction = 0.125;	

// set the adaptive time step
precision set_time_step(precision t_max, precision dt_min)
{
	precision CFL_factor = 0.125;	

	precision dt = fmin(fraction, CFL_factor) * t_max;

	if(dt < dt_min)
	{
		printf("set_time_step error: dt = %lf < %lf\n", dt, dt_min);
	}

	return dt_min * fmax(1., floor(dt / dt_min));	// round down dt to numerical precision
}


inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


// max time step to stabilize time evolution
precision compute_max_time_step(precision t, const hydro_variables * const __restrict__ q, const precision * const __restrict__ e, const fluid_velocity * const __restrict__ u, const fluid_velocity * const __restrict__ up, int nx, int ny, int nz, precision dt, precision dt_prev, precision dx, precision dy, precision dn, precision etabar_const)
{
	precision dt_max = 1./0.;		// default value for max time step

	precision t2 = t * t;
	int stride_y = nx + 4;							// strides for neighbor cells along x, y, n (stride_x = 1)
	int stride_z = (nx + 4) * (ny + 4);				// stride formulas based from linear_column_index()

	precision q_s[NUMBER_CONSERVED_VARIABLES];		// current variables at cell s
	precision   S[NUMBER_CONSERVED_VARIABLES];		// external source terms in hydrodynamic equations

	precision e1[6];								// primary variables of neighbors [i-1, i+1, j-1, j+1, k-1, k+1]

	precision qi1[2 * NUMBER_CONSERVED_VARIABLES];	// conserved variables of neighbor cells along x [i-1, i+1]
	precision qj1[2 * NUMBER_CONSERVED_VARIABLES];	// conserved variables of neighbor cells along y [j-1, j+1]
	precision qk1[2 * NUMBER_CONSERVED_VARIABLES];	// conserved variables of neighbor cells along n [k-1, k+1]

	precision qi2[2 * NUMBER_CONSERVED_VARIABLES];	// conserved variables of neighbor cells along x [i-2, i+2]
	precision qj2[2 * NUMBER_CONSERVED_VARIABLES];	// conserved variables of neighbor cells along y [j-2, j+2]
	precision qk2[2 * NUMBER_CONSERVED_VARIABLES];	// conserved variables of neighbor cells along n [k-2, k+2]

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

				q_s[0] = q[s].ttt;			// hydro variables of cell s
				q_s[1] = q[s].ttx;
				q_s[2] = q[s].tty;
				q_s[3] = q[s].ttn;

				int a = 4;					// counter
			#ifdef ANISO_HYDRO
				q_s[a] = q[s].pl;	a++;
			#if (PT_MATCHING == 1)
				q_s[a] = q[s].pt;	a++;
			#endif
			#endif

			#ifdef PIMUNU
				q_s[a] = q[s].pitt; a++;
				q_s[a] = q[s].pitx; a++;
				q_s[a] = q[s].pity; a++;
				q_s[a] = q[s].pitn; a++;
				q_s[a] = q[s].pixx; a++;
				q_s[a] = q[s].pixy; a++;
				q_s[a] = q[s].pixn; a++;
				q_s[a] = q[s].piyy; a++;
				q_s[a] = q[s].piyn; a++;
				q_s[a] = q[s].pinn; a++;
			#endif

			#ifdef WTZMU
				q_s[a] = q[s].WtTz; a++;
				q_s[a] = q[s].WxTz; a++;
				q_s[a] = q[s].WyTz; a++;
				q_s[a] = q[s].WnTz;
			#endif

			#ifdef PI
				q_s[a] = q[s].Pi;
			#endif

				precision e_s = e[s];

				precision ux = u[s].ux;
				precision uy = u[s].uy;
				precision un = u[s].un;
				precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t2 * un * un);

				precision ux_p = up[s].ux;	// previous fluid velocity
				precision uy_p = up[s].uy;
				precision un_p = up[s].un;

			#if (PT_MATCHING == 0)
				get_primary_neighbor_cells(e, e1, sim, sip, sjm, sjp, skm, skp);
			#endif

				get_fluid_velocity_neighbor_cells(u[simm], u[sim], u[sip], u[sipp], u[sjmm], u[sjm], u[sjp], u[sjpp], u[skmm], u[skm], u[skp], u[skpp], ui1, uj1, uk1, vxi, vyj, vnk, t2);

				get_hydro_neighbor_cells(q[sim],  q[sip],  qi1);
				get_hydro_neighbor_cells(q[simm], q[sipp], qi2);

				get_hydro_neighbor_cells(q[sjm],  q[sjp],  qj1);
				get_hydro_neighbor_cells(q[sjmm], q[sjpp], qj2);

				get_hydro_neighbor_cells(q[skm],  q[skp],  qk1);
				get_hydro_neighbor_cells(q[skmm], q[skpp], qk2);


				// update the max time step
			#ifdef ANISO_HYDRO
				dt_max = fmin(dt_max, source_terms_aniso_hydro(S, q_s, e_s, t, qi1, qj1, qk1, e1, ui1, uj1, uk1, ux, uy, un, ux_p, uy_p, un_p, dt_prev, dx, dy, dn, etabar_const));
			#else
				dt_max = fmin(dt_max, source_terms_viscous_hydro(S, q_s, e_s, t, qi1, qj1, qk1, e1, ui1, uj1, uk1, ux, uy, un, ux_p, uy_p, un_p, dt_prev, dx, dy, dn, etabar_const));
			#endif
			}
		}
	}

	return dt_max;
}






