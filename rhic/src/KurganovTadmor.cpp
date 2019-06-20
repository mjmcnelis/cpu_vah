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


void euler_step(PRECISION t, const CONSERVED_VARIABLES * const __restrict__ q_current, CONSERVED_VARIABLES * const __restrict__ q_updated,
const PRECISION * const __restrict__ e, const PRECISION * const __restrict__ p, const FLUID_VELOCITY * const __restrict__ u, const FLUID_VELOCITY * const __restrict__ up, int nx, int ny, int nz, int ncx, int ncy, PRECISION dt, PRECISION dx, PRECISION dy, PRECISION dn, PRECISION etabar)
{
	// compute the euler step
	int stride_x = 1;									// strides for neighbor cells along x, y, n
	int stride_y = ncx;									// stride formulas based from linear_column_index()
	int stride_z = ncx * ncy;

	PRECISION Q_s[NUMBER_CONSERVED_VARIABLES];			// current variables at cell s
	PRECISION S[NUMBER_CONSERVED_VARIABLES];			// external source terms in hydrodynamic equations

	PRECISION e1[6];									// primary variables of neighbors [i-1, i+1, j-1, j+1, k-1, k+1]
	PRECISION p1[6];

	PRECISION qi1[2 * NUMBER_CONSERVED_VARIABLES];		// conserved variables of neighbor cells along x [i-1, i+1]
	PRECISION qj1[2 * NUMBER_CONSERVED_VARIABLES];		// conserved variables of neighbor cells along y [j-1, j+1]
	PRECISION qk1[2 * NUMBER_CONSERVED_VARIABLES];		// conserved variables of neighbor cells along n [k-1, k+1]

	PRECISION qi2[2 * NUMBER_CONSERVED_VARIABLES];		// conserved variables of neighbor cells along x [i-2, i+2]
	PRECISION qj2[2 * NUMBER_CONSERVED_VARIABLES];		// conserved variables of neighbor cells along y [j-2, j+2]
	PRECISION qk2[2 * NUMBER_CONSERVED_VARIABLES];		// conserved variables of neighbor cells along n [k-2, k+2]

	PRECISION ui1[8];									// fluid velocity of neighbor cells along x [i-1, i+1]
	PRECISION uj1[8];									// fluid velocity of neighbor cells along y [j-1, j+1]
	PRECISION uk1[8];									// fluid velocity of neighbor cells along n [k-1, k+1]

	// for flux terms
	PRECISION vxi[4];									// vx of neighbor cells along x [i-2, i-1, i+1, i+2]
	PRECISION vyj[4];									// vy of neighbor cells along y [j-2, j-1, j+1, j+2]
	PRECISION vnk[4];									// vn of neighbor cells along n [k-2, k-1, k+1, k+2]

	PRECISION Hx_plus[NUMBER_CONSERVED_VARIABLES];		// Hx_{i + 1/2}
	PRECISION Hx_minus[NUMBER_CONSERVED_VARIABLES];		// Hx_{i - 1/2}

	PRECISION Hy_plus[NUMBER_CONSERVED_VARIABLES];		// Hy_{j + 1/2}
	PRECISION Hy_minus[NUMBER_CONSERVED_VARIABLES];		// Hy_{j - 1/2}

	PRECISION Hn_plus[NUMBER_CONSERVED_VARIABLES];		// Hn_{k + 1/2}
	PRECISION Hn_minus[NUMBER_CONSERVED_VARIABLES];		// Hn_{k - 1/2}


	// loop over physical grid points
	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, ncx, ncy);

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

				Q_s[0] = q_current->ttt[s];		// conserved variables of cell s
				Q_s[1] = q_current->ttx[s];
				Q_s[2] = q_current->tty[s];
				Q_s[3] = q_current->ttn[s];

				Q_s[4] = q_current->pl[s];		// todo: add pt

				PRECISION e_s = e[s];
				PRECISION p_s = p[s];

				PRECISION ut = u->ut[s];		// current fluid velocity
				PRECISION ux = u->ux[s];
				PRECISION uy = u->uy[s];
				PRECISION un = u->un[s];

				PRECISION ut_p = up->ut[s];		// previous fluid velocity
				PRECISION ux_p = up->ux[s];
				PRECISION uy_p = up->uy[s];
				PRECISION un_p = up->un[s];

				// get primary variables of neighbor cells
				get_primary_neighbor_cells(e, p, e1, p1, sim, sip, sjm, sjp, skm, skp);

				get_u_neighbor_cells(u->ut, u->ux, u->uy, u->un, ui1, uj1, uk1, sim, sip, sjm, sjp, skm, skp);

				// get spatial fluid velocity of neighbor cells
				get_v_neighbor_cells(u->ut, u->ux, u->uy, u->un, vxi, vyj, vnk, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);

				// get conserved variables of neighbor cells
				get_q_neighbor_cells(q_current->ttt, qi1, qj1, qk1, qi2, qj2, qk2, &r, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);
				get_q_neighbor_cells(q_current->ttx, qi1, qj1, qk1, qi2, qj2, qk2, &r, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);
				get_q_neighbor_cells(q_current->tty, qi1, qj1, qk1, qi2, qj2, qk2, &r, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);
				get_q_neighbor_cells(q_current->ttn, qi1, qj1, qk1, qi2, qj2, qk2, &r, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);
				get_q_neighbor_cells(q_current->pl,  qi1, qj1, qk1, qi2, qj2, qk2, &r, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);

			#ifdef PIMUNU
				// fill in
			#endif
			#ifdef W_TZ_MU
				// fill in
			#endif


				// compute the external source terms (S)
				source_terms(Q_s, S, u, ut_p, ux_p, uy_p, un_p, t, e, p, s, ncx, ncy, nz + 4, etabar, dt, dx, dy, dn, i, j, k, 0, 0, 0, q_current, qi1, qj1, qk1, e1, p1, e_s, p_s, ui1, uj1, uk1, ut, ux, uy, un, ut_p, ux_p, uy_p, un_p);

				// compute the flux terms Hx_plus, Hx_minus
				flux_terms(Q_s, qi1, qi2, vxi, ux, ut, Hx_plus, &rightHalfCellExtrapolationForward, &leftHalfCellExtrapolationForward);
				flux_terms(Q_s, qi1, qi2, vxi, ux, ut, Hx_minus, &rightHalfCellExtrapolationBackwards, &leftHalfCellExtrapolationBackwards);

				// compute the flux terms Hy_plus, Hy_minus
				flux_terms(Q_s, qj1, qj2, vyj, uy, ut, Hy_plus, &rightHalfCellExtrapolationForward, &leftHalfCellExtrapolationForward);
				flux_terms(Q_s, qj1, qj2, vyj, uy, ut, Hy_minus, &rightHalfCellExtrapolationBackwards, &leftHalfCellExtrapolationBackwards);

				// compute the flux terms Hn_plus, Hn_minus
				flux_terms(Q_s, qk1, qk2, vnk, un, ut, Hn_plus, &rightHalfCellExtrapolationForward, &leftHalfCellExtrapolationForward);
				flux_terms(Q_s, qk1, qk2, vnk, un, ut, Hn_minus, &rightHalfCellExtrapolationBackwards, &leftHalfCellExtrapolationBackwards);


				// add euler step
				for(short n = 0; n < NUMBER_CONSERVED_VARIABLES; n++)
				{
					Q_s[n] += dt * (S[n]  +  (Hx_minus[n] - Hx_plus[n]) / dx  +  (Hy_minus[n] - Hy_plus[n]) / dy  +  (Hn_minus[n] - Hn_plus[n]) / dn);
				}

				// update Q
				q_updated->ttt[s] = Q_s[0];
				q_updated->ttx[s] = Q_s[1];
				q_updated->tty[s] = Q_s[2];
				q_updated->ttn[s] = Q_s[3];

				q_updated->pl[s] = Q_s[4];

			#ifdef PIMUNU
				q_updated->pitt[s] = Q_s[5];
				q_updated->pitx[s] = Q_s[6];
				q_updated->pity[s] = Q_s[7];
				q_updated->pitn[s] = Q_s[8];
				q_updated->pixx[s] = Q_s[9];
				q_updated->pixy[s] = Q_s[10];
				q_updated->pixn[s] = Q_s[11];
				q_updated->piyy[s] = Q_s[12];
				q_updated->piyn[s] = Q_s[13];
				q_updated->pinn[s] = Q_s[14];
			#endif
			#ifdef W_TZ_MU
				q_updated->WtTz[s] = Q_s[15];
				q_updated->WxTz[s] = Q_s[16];
				q_updated->WyTz[s] = Q_s[17];
				q_updated->WnTz[s] = Q_s[18];
			#endif
			}
		}
	}
}


void convex_combination(const CONSERVED_VARIABLES * const __restrict__ q, CONSERVED_VARIABLES * const __restrict__ Q, int nx, int ny, int nz, int ncx, int ncy)
 {
 	// rungeKutta2 update: (q + Q) / 2  -> store in Q

 	// loop over physical grid points
	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, ncx, ncy);

				Q->ttt[s] = (q->ttt[s]  +  Q->ttt[s]) / 2.0;
				Q->ttx[s] = (q->ttx[s]  +  Q->ttx[s]) / 2.0;
				Q->tty[s] = (q->tty[s]  +  Q->tty[s]) / 2.0;
				Q->ttn[s] = (q->ttn[s]  +  Q->ttn[s]) / 2.0;

				Q->pl[s] = (q->pl[s]  +  Q->pl[s]) / 2.0;
				//Q->pt[s] = (q->pt[s]  +  Q->pt[s]) / 2.0;  // add pt

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
			#ifdef W_TZ_MU
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
void rungeKutta2(PRECISION t, PRECISION dt, CONSERVED_VARIABLES * __restrict__ q, CONSERVED_VARIABLES * __restrict__ Q, int nx, int ny, int nz, int ncx, int ncy, PRECISION dx, PRECISION dy, PRECISION dz, PRECISION etabar)
{
	// first intermediate time step (compute qS = q + dt.(S - dHx/dx - dHy/dy - dHz/dz))
	euler_step(t, q, qS, e, p, u, up, nx, ny, nz, ncx, ncy, dt, dx, dy, dz, etabar);

	// next time step
	t += dt;

	// compute uS, e, p
	set_inferred_variables(qS, e, p, uS, t, nx, ny, nz, ncx, ncy);

	// regulate dissipative components of qS
	regulate_dissipative_currents(t, qS, e, p, uS, nx, ny, nz, ncx, ncy);

	// set ghost cells for qS, uS, e, p
	set_ghost_cells(qS, e, p, uS, nx, ny, nz, ncx, ncy);

	// second intermediate time step (compute Q = qS + dt.(S - dHx/dx - dHy/dy - dHz/dz))
	euler_step(t, qS, Q, e, p, uS, u, nx, ny, nz, ncx, ncy, dt, dx, dy, dz, etabar);

	// Runge-Kutta: (q + Q) / 2 -> Q
	convex_combination(q, Q, nx, ny, nz, ncx, ncy);

	// swap up and u
	swap_fluid_velocity(&up, &u);

	// maybe this can all be grouped together in one function somehow
	// compute u, e, p
	set_inferred_variables(Q, e, p, u, t, nx, ny, nz, ncx, ncy);

	// regulate dissipative components of Q
	regulate_dissipative_currents(t, Q, e, p, u, nx, ny, nz, ncx, ncy);

	// set ghost cells for Q, u, e, p
	set_ghost_cells(Q, e, p, u, nx, ny, nz, ncx, ncy);

	// swap q and Q
	set_current_conserved_variables();
}




