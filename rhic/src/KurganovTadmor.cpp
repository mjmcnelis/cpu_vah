/*
 * FullyDiscreteKurganovTadmorScheme.cpp
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>

#include "../include/KurganovTadmor.h"
#include "../include/Parameters.h"
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"
#include "../include/GhostCells.h"
#include "../include/FluxTerms.h"
#include "../include/SourceTerms.h"
#include "../include/InferredVariables.h"
#include "../include/EquationOfState.h"
#include "../include/NeighborCells.h"
#include "../include/Regulation.h"


using namespace std;

inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


void euler_step(precision t, const CONSERVED_VARIABLES * const __restrict__ q, CONSERVED_VARIABLES * const __restrict__ Q,
const precision * const __restrict__ e, const FLUID_VELOCITY * const __restrict__ u, const FLUID_VELOCITY * const __restrict__ up, int nx, int ny, int nz, precision dt, precision dx, precision dy, precision dn, precision etabar)
{
	// compute the euler step
	int stride_x = 1;									// strides for neighbor cells along x, y, n
	int stride_y = nx + 4;								// stride formulas based from linear_column_index()
	int stride_z = (nx + 4) * (ny + 4);

	precision ql[NUMBER_CONSERVED_VARIABLES];			// current variables at cell s
	precision S[NUMBER_CONSERVED_VARIABLES];			// external source terms in hydrodynamic equations

	precision e1[6];									// primary variables of neighbors [i-1, i+1, j-1, j+1, k-1, k+1]

	precision qi1[2 * NUMBER_CONSERVED_VARIABLES];		// conserved variables of neighbor cells along x [i-1, i+1]
	precision qj1[2 * NUMBER_CONSERVED_VARIABLES];		// conserved variables of neighbor cells along y [j-1, j+1]
	precision qk1[2 * NUMBER_CONSERVED_VARIABLES];		// conserved variables of neighbor cells along n [k-1, k+1]

	precision qi2[2 * NUMBER_CONSERVED_VARIABLES];		// conserved variables of neighbor cells along x [i-2, i+2]
	precision qj2[2 * NUMBER_CONSERVED_VARIABLES];		// conserved variables of neighbor cells along y [j-2, j+2]
	precision qk2[2 * NUMBER_CONSERVED_VARIABLES];		// conserved variables of neighbor cells along n [k-2, k+2]

	precision ui1[8];									// fluid velocity of neighbor cells along x [i-1, i+1]
	precision uj1[8];									// fluid velocity of neighbor cells along y [j-1, j+1]
	precision uk1[8];									// fluid velocity of neighbor cells along n [k-1, k+1]

	// for flux terms
	precision vxi[4];									// vx of neighbor cells along x [i-2, i-1, i+1, i+2]
	precision vyj[4];									// vy of neighbor cells along y [j-2, j-1, j+1, j+2]
	precision vnk[4];									// vn of neighbor cells along n [k-2, k-1, k+1, k+2]

	precision Hx_plus[NUMBER_CONSERVED_VARIABLES];		// Hx_{i + 1/2}
	precision Hx_minus[NUMBER_CONSERVED_VARIABLES];		// Hx_{i - 1/2}

	precision Hy_plus[NUMBER_CONSERVED_VARIABLES];		// Hy_{j + 1/2}
	precision Hy_minus[NUMBER_CONSERVED_VARIABLES];		// Hy_{j - 1/2}

	precision Hn_plus[NUMBER_CONSERVED_VARIABLES];		// Hn_{k + 1/2}
	precision Hn_minus[NUMBER_CONSERVED_VARIABLES];		// Hn_{k - 1/2}


	// loop over physical grid points
	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				int sim  = s - stride_x;		// neighbor cell indices (x)
				int simm = sim - stride_x;
				int sip  = s + stride_x;
				int sipp = sip + stride_x;

				int sjm  = s - stride_y;		// neighbor cell indices (y)
				int sjmm = sjm - stride_y;
				int sjp  = s + stride_y;
				int sjpp = sjp + stride_y;

				int skm  = s - stride_z;		// neighbor cell indices (n)
				int skmm = skm - stride_z;
				int skp  = s + stride_z;
				int skpp = skp + stride_z;

				int r = 0;

				ql[0] = q->ttt[s];				// conserved variables of cell s
				ql[1] = q->ttx[s];
				ql[2] = q->tty[s];
				ql[3] = q->ttn[s];

				ql[4] = q->pl[s];

			#if (PT_MATCHING == 1)
				ql[5] = q->pt[s];		
			#endif

				precision e_s = e[s];

				precision ut = u->ut[s];		// current fluid velocity
				precision ux = u->ux[s];
				precision uy = u->uy[s];
				precision un = u->un[s];

				precision ut_p = up->ut[s];		// previous fluid velocity
				precision ux_p = up->ux[s];
				precision uy_p = up->uy[s];
				precision un_p = up->un[s];


				// get primary variables of neighbor cells
			#if (PT_MATCHING == 0)
				get_primary_neighbor_cells(e, e1, sim, sip, sjm, sjp, skm, skp);
			#endif

				get_u_neighbor_cells(u->ut, u->ux, u->uy, u->un, ui1, uj1, uk1, sim, sip, sjm, sjp, skm, skp);


				// get spatial fluid velocity of neighbor cells
				get_v_neighbor_cells(u->ut, u->ux, u->uy, u->un, vxi, vyj, vnk, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);

				// get conserved variables of neighbor cells
				get_q_neighbor_cells(q->ttt, qi1, qj1, qk1, qi2, qj2, qk2, &r, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);
				get_q_neighbor_cells(q->ttx, qi1, qj1, qk1, qi2, qj2, qk2, &r, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);
				get_q_neighbor_cells(q->tty, qi1, qj1, qk1, qi2, qj2, qk2, &r, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);
				get_q_neighbor_cells(q->ttn, qi1, qj1, qk1, qi2, qj2, qk2, &r, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);

				get_q_neighbor_cells(q->pl,  qi1, qj1, qk1, qi2, qj2, qk2, &r, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);

			#if (PT_MATCHING == 1)
				get_q_neighbor_cells(q->pt,  qi1, qj1, qk1, qi2, qj2, qk2, &r, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);
			#endif


				// compute the external source terms (S)
				source_terms(S, ql, e_s, t, qi1, qj1, qk1, e1, ui1, uj1, uk1, ut, ux, uy, un, ut_p, ux_p, uy_p, un_p, dt, dx, dy, dn, etabar);


				// compute the flux terms Hx_plus, Hx_minus
				flux_terms(Hx_plus, ql, qi1, qi2, vxi, ux, ut, &right_half_cell_extrapolation_forward, &left_half_cell_extrapolation_forward);
				flux_terms(Hx_minus, ql, qi1, qi2, vxi, ux, ut, &right_half_cell_extrapolation_backward, &left_half_cell_extrapolation_backward);


				// compute the flux terms Hy_plus, Hy_minus
				flux_terms(Hy_plus, ql, qj1, qj2, vyj, uy, ut, &right_half_cell_extrapolation_forward, &left_half_cell_extrapolation_forward);
				flux_terms(Hy_minus, ql, qj1, qj2, vyj, uy, ut, &right_half_cell_extrapolation_backward, &left_half_cell_extrapolation_backward);


				// compute the flux terms Hn_plus, Hn_minus
				flux_terms(Hn_plus, ql, qk1, qk2, vnk, un, ut, &right_half_cell_extrapolation_forward, &left_half_cell_extrapolation_forward);
				flux_terms(Hn_minus, ql, qk1, qk2, vnk, un, ut, &right_half_cell_extrapolation_backward, &left_half_cell_extrapolation_backward);


				// add euler step
				for(short n = 0; n < NUMBER_CONSERVED_VARIABLES; n++)
				{
					// question: is last term 1/(t.dn) or 1/dn?
					ql[n] += dt * (S[n]  +  (Hx_minus[n] - Hx_plus[n]) / dx  +  (Hy_minus[n] - Hy_plus[n]) / dy  +  (Hn_minus[n] - Hn_plus[n]) / dn);

					//ql[n] += dt * (S[n]  +  (Hx_minus[n] - Hx_plus[n]) / dx  +  (Hy_minus[n] - Hy_plus[n]) / dy)  +  (Hn_minus[n] - Hn_plus[n]) / dn;
				}


				// update Q
				Q->ttt[s] = ql[0];
				Q->ttx[s] = ql[1];
				Q->tty[s] = ql[2];
				Q->ttn[s] = ql[3];

				Q->pl[s] = ql[4];

			#if (PT_MATCHING == 1)
				Q->pt[s] = ql[5];
			#endif

			#ifdef PIMUNU
				Q->pitt[s] = ql[5];
				Q->pitx[s] = ql[6];
				Q->pity[s] = ql[7];
				Q->pitn[s] = ql[8];
				Q->pixx[s] = ql[9];
				Q->pixy[s] = ql[10];
				Q->pixn[s] = ql[11];
				Q->piyy[s] = ql[12];
				Q->piyn[s] = ql[13];
				Q->pinn[s] = ql[14];
			#endif
			#ifdef WTZMU
				Q->WtTz[s] = ql[15];
				Q->WxTz[s] = ql[16];
				Q->WyTz[s] = ql[17];
				Q->WnTz[s] = ql[18];
			#endif
			}
		}
	}
}


void convex_combination(const CONSERVED_VARIABLES * const __restrict__ q, CONSERVED_VARIABLES * const __restrict__ Q, int nx, int ny, int nz)
 {
 	// rungeKutta2 update: (q + Q) / 2  -> store in Q

 	// loop over physical grid points
	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				Q->ttt[s] = (q->ttt[s]  +  Q->ttt[s]) / 2.0;
				Q->ttx[s] = (q->ttx[s]  +  Q->ttx[s]) / 2.0;
				Q->tty[s] = (q->tty[s]  +  Q->tty[s]) / 2.0;
				Q->ttn[s] = (q->ttn[s]  +  Q->ttn[s]) / 2.0;

				Q->pl[s] = (q->pl[s]  +  Q->pl[s]) / 2.0;

			#if (PT_MATCHING == 1)
				Q->pt[s] = (q->pt[s]  +  Q->pt[s]) / 2.0;
			#endif

			#ifdef PIMUNU
				Q->pitt[s] = (q->pitt[s]  +  Q->pitt[s]) / 2.0;
				Q->pitx[s] = (q->pitx[s]  +  Q->pitx[s]) / 2.0;
				Q->pity[s] = (q->pity[s]  +  Q->pity[s]) / 2.0;
				Q->pitn[s] = (q->pitn[s]  +  Q->pitn[s]) / 2.0;
				Q->pixx[s] = (q->pixx[s]  +  Q->pixx[s]) / 2.0;
				Q->pixy[s] = (q->pixy[s]  +  Q->pixy[s]) / 2.0;
				Q->pixn[s] = (q->pixn[s]  +  Q->pixn[s]) / 2.0;
				Q->piyy[s] = (q->piyy[s]  +  Q->piyy[s]) / 2.0;
				Q->piyn[s] = (q->piyn[s]  +  Q->piyn[s]) / 2.0;
				Q->pinn[s] = (q->pinn[s]  +  Q->pinn[s]) / 2.0;
			#endif
			#ifdef WTZMU
				Q->WtTz[s] = (q->WtTz[s]  +  Q->WtTz[s]) / 2.0;
				Q->WxTz[s] = (q->WxTz[s]  +  Q->WxTz[s]) / 2.0;
				Q->WyTz[s] = (q->WyTz[s]  +  Q->WyTz[s]) / 2.0;
				Q->WnTz[s] = (q->WnTz[s]  +  Q->WnTz[s]) / 2.0;
			#endif
			}
		}
	}
}


// main algorithm
void runge_kutta2(precision t, precision dt, int nx, int ny, int nz, precision dx, precision dy, precision dz, precision etabar)
{
	// first intermediate time step (compute qS = q + dt.(S - dHx/dx - dHy/dy - dHz/dz))
	euler_step(t, q, qS, e, u, up, nx, ny, nz, dt, dx, dy, dz, etabar);

	// next time step
	t += dt;

	// compute uS, e,
	set_inferred_variables(qS, e, uS, t, nx, ny, nz);

	// regulate dissipative components of qS
	regulate_dissipative_currents(t, qS, e, uS, nx, ny, nz);

	// set ghost cells for qS, uS, e,
	set_ghost_cells(qS, e, uS, nx, ny, nz);

	// second intermediate time step (compute Q = qS + dt.(S - dHx/dx - dHy/dy - dHz/dz))
	euler_step(t, qS, Q, e, uS, u, nx, ny, nz, dt, dx, dy, dz, etabar);

	// Runge-Kutta: (q + Q) / 2 -> Q
	convex_combination(q, Q, nx, ny, nz);

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
	set_current_conserved_variables();
}




