#include <stdlib.h>
#include <math.h>
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"
#include "../include/GhostCells.h"
#include "../include/FluxTerms.h"
#include "../include/SourceTerms.h"
#include "../include/InferredVariables.h"
#include "../include/NeighborCells.h"
#include "../include/Regulation.h"

inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


void euler_step(precision t, const hydro_variables * const __restrict__ q, hydro_variables * const __restrict__ Q,
const precision * const __restrict__ e, const fluid_velocity * const __restrict__ u, const fluid_velocity * const __restrict__ up, int nx, int ny, int nz, precision dt, precision dx, precision dy, precision dn, precision etabar)
{
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

	precision  Hx_plus[NUMBER_CONSERVED_VARIABLES];	// Hx_{i + 1/2}
	precision Hx_minus[NUMBER_CONSERVED_VARIABLES];	// Hx_{i - 1/2}

	precision  Hy_plus[NUMBER_CONSERVED_VARIABLES];	// Hy_{j + 1/2}
	precision Hy_minus[NUMBER_CONSERVED_VARIABLES];	// Hy_{j - 1/2}

	precision  Hn_plus[NUMBER_CONSERVED_VARIABLES];	// Hn_{k + 1/2}
	precision Hn_minus[NUMBER_CONSERVED_VARIABLES];	// Hn_{k - 1/2}

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

				q_s[0] = q[s].ttt;			// conserved variables of cell s
				q_s[1] = q[s].ttx;
				q_s[2] = q[s].tty;
				q_s[3] = q[s].ttn;
				q_s[4] = q[s].pl;

				int a = 5;					// counter
			#if (PT_MATCHING == 1)
				q_s[a] = q[s].pt;	a++;
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
				get_conserved_neighbor_cells(q[simm], q[sim], q[sip], q[sipp], q[sjmm], q[sjm], q[sjp], q[sjpp], q[skmm], q[skm], q[skp], q[skpp], qi1, qj1, qk1, qi2, qj2, qk2);

				// compute the external source terms (S)
				source_terms(S, q_s, e_s, t, qi1, qj1, qk1, e1, ui1, uj1, uk1, ux, uy, un, ux_p, uy_p, un_p, dt, dx, dy, dn, etabar);

				// compute the flux terms
				flux_terms(Hx_plus, Hx_minus, q_s, qi1, qi2, vxi, ux/ut);
				flux_terms(Hy_plus, Hy_minus, q_s, qj1, qj2, vyj, uy/ut);
				flux_terms(Hn_plus, Hn_minus, q_s, qk1, qk2, vnk, un/ut);

				// add euler step
				for(int n = 0; n < NUMBER_CONSERVED_VARIABLES; n++)
				{
					q_s[n] += dt * (S[n]  +  (Hx_minus[n] - Hx_plus[n]) / dx  +  (Hy_minus[n] - Hy_plus[n]) / dy  +  (Hn_minus[n] - Hn_plus[n]) / dn);
				}

				// update Q
				Q[s].ttt = q_s[0];
				Q[s].ttx = q_s[1];
				Q[s].tty = q_s[2];
				Q[s].ttn = q_s[3];
				Q[s].pl  = q_s[4];

				a = 5;	// reset counter
			#if (PT_MATCHING == 1)
				Q[s].pt = q_s[a]; a++;
			#endif
			#ifdef PIMUNU
				Q[s].pitt = q_s[a]; a++;
				Q[s].pitx = q_s[a]; a++;
				Q[s].pity = q_s[a]; a++;
				Q[s].pitn = q_s[a]; a++;
				Q[s].pixx = q_s[a]; a++;
				Q[s].pixy = q_s[a]; a++;
				Q[s].pixn = q_s[a]; a++;
				Q[s].piyy = q_s[a]; a++;
				Q[s].piyn = q_s[a]; a++;
				Q[s].pinn = q_s[a]; a++;
			#endif
			#ifdef WTZMU
				Q[s].WtTz = q_s[a]; a++;
				Q[s].WxTz = q_s[a]; a++;
				Q[s].WyTz = q_s[a]; a++;
				Q[s].WnTz = q_s[a];
			#endif
			}
		}
	}
}


void runge_kutta2(const hydro_variables * const __restrict__ q, hydro_variables * const __restrict__ Q, int nx, int ny, int nz)
 {
 	// runge kutta time evolution update: (q + Q) / 2  -> store in Q
	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				Q[s].ttt = (q[s].ttt  +  Q[s].ttt) / 2.;
				Q[s].ttx = (q[s].ttx  +  Q[s].ttx) / 2.;
				Q[s].tty = (q[s].tty  +  Q[s].tty) / 2.;
				Q[s].ttn = (q[s].ttn  +  Q[s].ttn) / 2.;
				Q[s].pl  = (q[s].pl   +  Q[s].pl)  / 2.;
			#if (PT_MATCHING == 1)
				Q[s].pt  = (q[s].pt   +  Q[s].pt)  / 2.;
			#endif
			#ifdef PIMUNU
				Q[s].pitt = (q[s].pitt  +  Q[s].pitt) / 2.;
				Q[s].pitx = (q[s].pitx  +  Q[s].pitx) / 2.;
				Q[s].pity = (q[s].pity  +  Q[s].pity) / 2.;
				Q[s].pitn = (q[s].pitn  +  Q[s].pitn) / 2.;
				Q[s].pixx = (q[s].pixx  +  Q[s].pixx) / 2.;
				Q[s].pixy = (q[s].pixy  +  Q[s].pixy) / 2.;
				Q[s].pixn = (q[s].pixn  +  Q[s].pixn) / 2.;
				Q[s].piyy = (q[s].piyy  +  Q[s].piyy) / 2.;
				Q[s].piyn = (q[s].piyn  +  Q[s].piyn) / 2.;
				Q[s].pinn = (q[s].pinn  +  Q[s].pinn) / 2.;
			#endif
			#ifdef WTZMU
				Q[s].WtTz = (q[s].WtTz  +  Q[s].WtTz) / 2.;
				Q[s].WxTz = (q[s].WxTz  +  Q[s].WxTz) / 2.;
				Q[s].WyTz = (q[s].WyTz  +  Q[s].WyTz) / 2.;
				Q[s].WnTz = (q[s].WnTz  +  Q[s].WnTz) / 2.;
			#endif
			}
		}
	}
}


void evolve_hydro_one_time_step(precision t, precision dt, int nx, int ny, int nz, precision dx, precision dy, precision dz, precision etabar)
{
	// first intermediate time step (compute qS = q + dt.(S - dHx/dx - dHy/dy - dHn/dn))
	euler_step(t, q, qI, e, u, up, nx, ny, nz, dt, dx, dy, dz, etabar);

	// next time step
	t += dt;

	// compute uS, e
	set_inferred_variables(qI, e, uI, t, nx, ny, nz);

	// regulate dissipative components of qS
	regulate_dissipative_currents(t, qI, e, uI, nx, ny, nz);	

	// set ghost cells for qS, uS, e,
	set_ghost_cells(qI, e, uI, nx, ny, nz);

	// second intermediate time step (compute Q = qS + dt.(S - dHx/dx - dHy/dy - dHn/dn))
	euler_step(t, qI, Q, e, uI, u, nx, ny, nz, dt, dx, dy, dz, etabar);

	// 2nd order runge kutta: (q + Q) / 2 -> Q
	runge_kutta2(q, Q, nx, ny, nz);

	// swap up and u
	swap_fluid_velocity(&up, &u);

	// maybe this can all be grouped together in one function somehow
	// compute u, e,
	set_inferred_variables(Q, e, u, t, nx, ny, nz);

	// regulate dissipative components of Q
	regulate_dissipative_currents(t, Q, e, u, nx, ny, nz);

	// set ghost cells for Q, u, e,
	set_ghost_cells(Q, e, u, nx, ny, nz);

	// swap q and Q
	set_current_hydro_variables();
}




