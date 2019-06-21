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


void euler_step(precision t, const CONSERVED_VARIABLES * const __restrict__ q_current, CONSERVED_VARIABLES * const __restrict__ q_updated,
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

				ql[0] = q_current->ttt[s];		// conserved variables of cell s
				ql[1] = q_current->ttx[s];
				ql[2] = q_current->tty[s];
				ql[3] = q_current->ttn[s];

				ql[4] = q_current->pl[s];
				ql[5] = q_current->pt[s];		// added pt

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
				get_primary_neighbor_cells(e, e1, sim, sip, sjm, sjp, skm, skp);

				get_u_neighbor_cells(u->ut, u->ux, u->uy, u->un, ui1, uj1, uk1, sim, sip, sjm, sjp, skm, skp);

				// get spatial fluid velocity of neighbor cells
				get_v_neighbor_cells(u->ut, u->ux, u->uy, u->un, vxi, vyj, vnk, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);

				// get conserved variables of neighbor cells
				get_q_neighbor_cells(q_current->ttt, qi1, qj1, qk1, qi2, qj2, qk2, &r, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);
				get_q_neighbor_cells(q_current->ttx, qi1, qj1, qk1, qi2, qj2, qk2, &r, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);
				get_q_neighbor_cells(q_current->tty, qi1, qj1, qk1, qi2, qj2, qk2, &r, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);
				get_q_neighbor_cells(q_current->ttn, qi1, qj1, qk1, qi2, qj2, qk2, &r, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);
				get_q_neighbor_cells(q_current->pl,  qi1, qj1, qk1, qi2, qj2, qk2, &r, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);


				get_q_neighbor_cells(q_current->pt, qi1, qj1, qk1, qi2, qj2, qk2, &r, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);





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
				q_updated->ttt[s] = ql[0];
				q_updated->ttx[s] = ql[1];
				q_updated->tty[s] = ql[2];
				q_updated->ttn[s] = ql[3];

				q_updated->pl[s] = ql[4];
				q_updated->pt[s] = ql[5];

			#ifdef PIMUNU
				q_updated->pitt[s] = ql[5];
				q_updated->pitx[s] = ql[6];
				q_updated->pity[s] = ql[7];
				q_updated->pitn[s] = ql[8];
				q_updated->pixx[s] = ql[9];
				q_updated->pixy[s] = ql[10];
				q_updated->pixn[s] = ql[11];
				q_updated->piyy[s] = ql[12];
				q_updated->piyn[s] = ql[13];
				q_updated->pinn[s] = ql[14];
			#endif
			#ifdef WTZMU
				q_updated->WtTz[s] = ql[15];
				q_updated->WxTz[s] = ql[16];
				q_updated->WyTz[s] = ql[17];
				q_updated->WnTz[s] = ql[18];
			#endif
			}
		}
	}
}


void convex_combination(const CONSERVED_VARIABLES * const __restrict__ ql, CONSERVED_VARIABLES * const __restrict__ Ql, int nx, int ny, int nz)
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

				Ql->ttt[s] = (ql->ttt[s]  +  Ql->ttt[s]) / 2.0;
				Ql->ttx[s] = (ql->ttx[s]  +  Ql->ttx[s]) / 2.0;
				Ql->tty[s] = (ql->tty[s]  +  Ql->tty[s]) / 2.0;
				Ql->ttn[s] = (ql->ttn[s]  +  Ql->ttn[s]) / 2.0;

				Ql->pl[s] = (ql->pl[s]  +  Ql->pl[s]) / 2.0;
				Ql->pt[s] = (ql->pt[s]  +  Ql->pt[s]) / 2.0;

			#ifdef PIMUNU
				Ql->pitt[s] = (ql->pitt[s]  +  Ql->pitt[s]) / 2.0;
				Ql->pitx[s] = (ql->pitx[s]  +  Ql->pitx[s]) / 2.0;
				Ql->pity[s] = (ql->pity[s]  +  Ql->pity[s]) / 2.0;
				Ql->pitn[s] = (ql->pitn[s]  +  Ql->pitn[s]) / 2.0;
				Ql->pixx[s] = (ql->pixx[s]  +  Ql->pixx[s]) / 2.0;
				Ql->pixy[s] = (ql->pixy[s]  +  Ql->pixy[s]) / 2.0;
				Ql->pixn[s] = (ql->pixn[s]  +  Ql->pixn[s]) / 2.0;
				Ql->piyy[s] = (ql->piyy[s]  +  Ql->piyy[s]) / 2.0;
				Ql->piyn[s] = (ql->piyn[s]  +  Ql->piyn[s]) / 2.0;
				Ql->pinn[s] = (ql->pinn[s]  +  Ql->pinn[s]) / 2.0;
			#endif
			#ifdef WTZMU
				Ql->WtTz[s] = (ql->WtTz[s]  +  Ql->WtTz[s]) / 2.0;
				Ql->WxTz[s] = (ql->WxTz[s]  +  Ql->WxTz[s]) / 2.0;
				Ql->WyTz[s] = (ql->WyTz[s]  +  Ql->WyTz[s]) / 2.0;
				Ql->WnTz[s] = (ql->WnTz[s]  +  Ql->WnTz[s]) / 2.0;
			#endif
			}
		}
	}
}


// main algorithm
// void rungeKutta2(precision t, precision dt, CONSERVED_VARIABLES * __restrict__ q, CONSERVED_VARIABLES * __restrict__ Q, int nx, int ny, int nz, int ncx, int ncy, precision dx, precision dy, precision dz, precision etabar)
void rungeKutta2(precision t, precision dt, int nx, int ny, int nz, precision dx, precision dy, precision dz, precision etabar)
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




