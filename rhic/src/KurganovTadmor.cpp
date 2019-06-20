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
#include "../include/FluxTerms.h"
#include "../include/SourceTerms.h"
#include "../include/EnergyMomentumTensor.h"
#include "../include/EquationOfState.h"
#include "../include/AnisotropicDistributionFunctions.h"
#include "../include/NeighborCells.h"

using namespace std;

inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}

void setNeighborCellsJK2(const PRECISION * const __restrict__ Q, PRECISION * const __restrict__ I, int s, int n, int smm, int sm, int sp, int spp)
{
	//PRECISION data_ns = in[s];
	// *(I + n)      = Q[smm];		// store Q[s] + neighbor cells [s - 2, s + 2] in I
	// *(I + n + 1)  = Q[sm];
	// *(I + n + 2)  = Q[s];
	// *(I + n + 3)  = Q[sp];
	// *(I + n + 4)  = Q[spp];

	// what's the difference?
	I[n]      = Q[smm];			// store Q[s] + neighbor cells [s - 2, s + 2] in I
	I[n + 1]  = Q[sm];
	I[n + 2]  = Q[s];
	I[n + 3]  = Q[sp];
	I[n + 4]  = Q[spp];
}

void euler_step(PRECISION t, const CONSERVED_VARIABLES * const __restrict__ q_current, CONSERVED_VARIABLES * const __restrict__ q_updated,
const PRECISION * const __restrict__ e, const PRECISION * const __restrict__ p, const FLUID_VELOCITY * const __restrict__ u, const FLUID_VELOCITY * const __restrict__ up, int nx, int ny, int nz, int ncx, int ncy, PRECISION dt, PRECISION dx, PRECISION dy, PRECISION dn, PRECISION etabar)
{
	// compute the euler step

	int nxR = nx + 2;									// right physical boundary indices
	int nyR = ny + 2;
	int nzR = nz + 2;

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



	PRECISION I[5 * NUMBER_CONSERVED_VARIABLES];		// old formula
	PRECISION J[5 * NUMBER_CONSERVED_VARIABLES];
	PRECISION K[5 * NUMBER_CONSERVED_VARIABLES];


	// loop over physical grid points
	for(int k = 2; k < nzR; k++)
	{
		for(int j = 2; j < nyR; j++)
		{
			for(int i = 2; i < nxR; i++)
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
			/*
			#ifdef PIMUNU
				Q[5] = Q_current->pitt[s];		// there are a lot of variables
				Q[6] = Q_current->pitx[s];		// is it possible to use independent components?
				Q[7] = Q_current->pity[s];		// for (e, p, u) reconstruction formulas
				Q[8] = Q_current->pitn[s];		// only need pitt, pitx, pity, pitn, pixx, pixy, Wt, Wx, Wy, Wn maybe?
				Q[9] = Q_current->pixx[s];
				Q[10] = Q_current->pixy[s];		// can reconstruct pixn, piyn, pinn
				Q[11] = Q_current->pixn[s];		// still leaves me with 16 variables
				Q[12] = Q_current->piyy[s];		// it may be possible to use the orthogonal formulas
				Q[13] = Q_current->piyn[s];		// but it sounds daunting
				Q[14] = Q_current->pinn[s];

				get_Q_neighbor_cells(Q_current->pitt, Q_I, Q_J, Q_K, &m, s, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);
				get_Q_neighbor_cells(Q_current->pitx, Q_I, Q_J, Q_K, &m, s, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);
				get_Q_neighbor_cells(Q_current->pity, Q_I, Q_J, Q_K, &m, s, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);
				get_Q_neighbor_cells(Q_current->pitn, Q_I, Q_J, Q_K, &m, s, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);
				get_Q_neighbor_cells(Q_current->pixx, Q_I, Q_J, Q_K, &m, s, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);
				get_Q_neighbor_cells(Q_current->pixy, Q_I, Q_J, Q_K, &m, s, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);
				get_Q_neighbor_cells(Q_current->pixn, Q_I, Q_J, Q_K, &m, s, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);
				get_Q_neighbor_cells(Q_current->piyy, Q_I, Q_J, Q_K, &m, s, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);
				get_Q_neighbor_cells(Q_current->piyn, Q_I, Q_J, Q_K, &m, s, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);
				get_Q_neighbor_cells(Q_current->pinn, Q_I, Q_J, Q_K, &m, s, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);
			#endif
			#ifdef W_TZ_MU
				Q[15] = Q_current->WtTz[s];
				Q[16] = Q_current->WxTz[s];
				Q[17] = Q_current->WyTz[s];
				Q[18] = Q_current->WnTz[s];

				get_Q_neighbor_cells(Q_current->WtTz, Q_I, Q_J, Q_K, &m, s, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);
				get_Q_neighbor_cells(Q_current->WxTz, Q_I, Q_J, Q_K, &m, s, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);
				get_Q_neighbor_cells(Q_current->WyTz, Q_I, Q_J, Q_K, &m, s, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);
				get_Q_neighbor_cells(Q_current->WnTz, Q_I, Q_J, Q_K, &m, s, simm, sim, sip, sipp, sjmm, sjm, sjp, sjpp, skmm, skm, skp, skpp);
			#endif
			*/

				int ptr = 0;

				setNeighborCellsJK2(q_current->ttt, I, s, ptr, simm, sim, sip, sipp);
				setNeighborCellsJK2(q_current->ttt, J, s, ptr, sjmm, sjm, sjp, sjpp);
				setNeighborCellsJK2(q_current->ttt, K, s, ptr, skmm, skm, skp, skpp);
				ptr += 5;
				setNeighborCellsJK2(q_current->ttx, I, s, ptr, simm, sim, sip, sipp);
				setNeighborCellsJK2(q_current->ttx, J, s, ptr, sjmm, sjm, sjp, sjpp);
				setNeighborCellsJK2(q_current->ttx, K, s, ptr, skmm, skm, skp, skpp);
				ptr += 5;
				setNeighborCellsJK2(q_current->tty, I, s, ptr, simm, sim, sip, sipp);
				setNeighborCellsJK2(q_current->tty, J, s, ptr, sjmm, sjm, sjp, sjpp);
				setNeighborCellsJK2(q_current->tty, K, s, ptr, skmm, skm, skp, skpp);
				ptr += 5;
				setNeighborCellsJK2(q_current->ttn, I, s, ptr, simm, sim, sip, sipp);
				setNeighborCellsJK2(q_current->ttn, J, s, ptr, sjmm, sjm, sjp, sjpp);
				setNeighborCellsJK2(q_current->ttn, K, s, ptr, skmm, skm, skp, skpp);
				ptr += 5;
				setNeighborCellsJK2(q_current->pl, I, s, ptr, simm, sim, sip, sipp);
				setNeighborCellsJK2(q_current->pl, J, s, ptr, sjmm, sjm, sjp, sjpp);
				setNeighborCellsJK2(q_current->pl, K, s, ptr, skmm, skm, skp, skpp);
				

				




				// compute the external source terms (S)

				source_terms(Q_s, S, u, ut_p, ux_p, uy_p, un_p, t, e, p, s, ncx, ncy, nz + 4, etabar, dt, dx, dy, dn, i, j, k, 0, 0, 0, q_current, qi1, qj1, qk1, e1, p1, e_s, p_s, ui1, uj1, uk1, ut, ux, uy, un, ut_p, ux_p, uy_p, un_p);

				// try the old version
				//loadSourceTerms(Q_s, S, u, up->ut[s], up->ux[s], up->uy[s], up->un[s], t, e, p, s, ncx, ncy, nz + 4, etabar, dt, dx, dy, dn, i, j, k, 0, 0, 0, q_current);



				// compute the flux terms (Hx, Hy, Hn)
				// implementing an updated version of this algorithm

				// I should try to reproduce the old algorithm (using this function)
				flux_terms(Q_s, qi1, qi2, vxi, ux, ut, Hx_plus, &rightHalfCellExtrapolationForward, &leftHalfCellExtrapolationForward);
				flux_terms(Q_s, qi1, qi2, vxi, ux, ut, Hx_minus, &rightHalfCellExtrapolationBackwards, &leftHalfCellExtrapolationBackwards);

				flux_terms(Q_s, qj1, qj2, vyj, uy, ut, Hy_plus, &rightHalfCellExtrapolationForward, &leftHalfCellExtrapolationForward);
				flux_terms(Q_s, qj1, qj2, vyj, uy, ut, Hy_minus, &rightHalfCellExtrapolationBackwards, &leftHalfCellExtrapolationBackwards);

				flux_terms(Q_s, qk1, qk2, vnk, un, ut, Hn_plus, &rightHalfCellExtrapolationForward, &leftHalfCellExtrapolationForward);
				flux_terms(Q_s, qk1, qk2, vnk, un, ut, Hn_minus, &rightHalfCellExtrapolationBackwards, &leftHalfCellExtrapolationBackwards);
				




				// PRECISION uT = getTransverseFluidVelocityMagnitude(u, s);

				// //compute the flux terms (Hx_plus, Hx_minus)
				// flux(I, Hx_plus, &rightHalfCellExtrapolationForward, &leftHalfCellExtrapolationForward, &spectralRadiusX, &Fx, t, e[s],uT,i,j,k);
				// flux(I, Hx_minus, &rightHalfCellExtrapolationBackwards, &leftHalfCellExtrapolationBackwards, &spectralRadiusX, &Fx, t, e[s], uT,i,j,k);

				// flux(J, Hy_plus, &rightHalfCellExtrapolationForward, &leftHalfCellExtrapolationForward, &spectralRadiusY, &Fy, t, e[s],uT,i,j,k);
				// flux(J, Hy_minus, &rightHalfCellExtrapolationBackwards, &leftHalfCellExtrapolationBackwards, &spectralRadiusY, &Fy, t, e[s], uT,i,j,k);

				// flux(K, Hn_plus, &rightHalfCellExtrapolationForward, &leftHalfCellExtrapolationForward, &spectralRadiusZ, &Fz, t, e[s],uT,i,j,k);
				// flux(K, Hn_minus, &rightHalfCellExtrapolationBackwards, &leftHalfCellExtrapolationBackwards, &spectralRadiusZ, &Fz, t, e[s], uT,i,j,k);


// IDEAL is defined if turn PIMUNU, W_TZ_MU off
// figure out how to add these other source terms in source_terms
// #ifndef IDEAL
// 				// loadSourceTermsX(I, H, u, s, dx);
				// loadSourceTermsY(J, H, u, s, dy);
				// loadSourceTermsZ(K, H, u, s, dz);
// 				// for (unsigned int n = 0; n < NUMBER_CONSERVATION_LAWS; ++n)
// 				// {
// 				// 	*(flux_X + n) += *(H + n);
// 				// 	*(flux_X + n) *= dt;
// 				// }
// #endif

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
 	// rungeKutta2 update: (q + Q) / 2  -> store in Q (looks okay)
 	int nxR = nx + 2;
 	int nyR = ny + 2;	// right physical boundary indices
 	int nzR = nz + 2;

	for(int k = 2; k < nzR; k++)
	{
		for(int j = 2; j < nyR; j++)
		{
			for(int i = 2; i < nxR; i++)
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


void regulate_residual_dissipative_currents(PRECISION t, CONSERVED_VARIABLES * const __restrict__ Q_current, PRECISION * const __restrict__ e, PRECISION * const __restrict__ p,
const FLUID_VELOCITY * const __restrict__ u, int nx, int ny, int nz, int ncx, int ncy, PRECISION dx, PRECISION dy, PRECISION dz, VALIDITY_DOMAIN * const __restrict__ validityDomain)
{
	PRECISION eps = 1.e-7;

	PRECISION xi0 = 0.1;		// regulation parameters (maybe move these in hydro parameters?)
	PRECISION rho_max = 1.0;

	PRECISION t2 = t * t;

	int nxR = nx + 2;			// right physical boundary indices
 	int nyR = ny + 2;
 	int nzR = nz + 2;

	for(int k = 2; k < nzR; k++)
	{
		for(int j = 2; j < nyR; j++)
		{
			for(int i = 2; i < nxR; i++)
			{
				int s = linear_column_index(i, j, k, ncx, ncy);

				PRECISION e_s = e[s];
				PRECISION p_s = p[s];
				PRECISION pl = Q_current->pl[s];

				if(e_s < 0.0)
				{
					e[s] = eps;
					p[s] = eps;
				}
				if(pl  < 0.0) 
				{
					Q_current->pl[s] = eps;
				}

				PRECISION pt = transversePressureHat(e_s, p_s, pl);		// temporary

				PRECISION ut = u->ut[s];
				PRECISION ux = u->ux[s];
				PRECISION uy = u->uy[s];
				PRECISION un = u->un[s];

				PRECISION utperp = sqrt(1.0  +  ux * ux  +  uy * uy);
				PRECISION zt = t * un / utperp;
				PRECISION zn = ut / t / utperp;

				PRECISION t2un = t2 * un;
				PRECISION t2zn = t2 * zn;

				// sqrt(T_aniso.T_aniso)
				PRECISION T_aniso_mag = sqrt(e_s * e_s  +  pl * pl  +  2.0 * pt * pt);

				// regulate transverse shear stress
			#ifdef PIMUNU
				PRECISION pitt = Q_current->pitt[s];
				PRECISION pitx = Q_current->pitx[s];
				PRECISION pity = Q_current->pity[s];
				PRECISION pitn = Q_current->pitn[s];
				PRECISION pixx = Q_current->pixx[s];
				PRECISION pixy = Q_current->pixy[s];
				PRECISION pixn = Q_current->pixn[s];
				PRECISION piyy = Q_current->piyy[s];
				PRECISION piyn = Q_current->piyn[s];
				PRECISION pinn = Q_current->pinn[s];


				// should I enforce tracelessness and orthogonality first
				// then regulate with a0?

				// sqrt(pi.pi), Tr(pi), u.pi and z.pi
				PRECISION pi_mag = sqrt(fabs(pitt * pitt  +  pixx * pixx  +  piyy * piyy  +  t2 * t2 * pinn * pinn  -  2.0 * (pitx * pitx  +  pity * pity  -  pixy * pixy  +  t2 * (pitn * pitn  -  pixn * pixn  -  piyn * piyn))));

				PRECISION Trpi = pitt  -  pixx - piyy -  t2 * pinn;

				PRECISION piu0 = pitt * ut  -  pitx * ux  -  pity * uy  -  pitn * t2un;
				PRECISION piu1 = pitx * ut  -  pixx * ux  -  pixy * uy  -  pixn * t2un;
				PRECISION piu2 = pity * ut  -  pixy * ux  -  piyy * uy  -  piyn * t2un;
				PRECISION piu3 = pitn * ut  -  pixn * ux  -  piyn * uy  -  pinn * t2un;

				PRECISION piz0 = zt * pitt  - t2zn * pitn;
				PRECISION piz1 = zt * pitx  - t2zn * pixn;
				PRECISION piz2 = zt * pity  - t2zn * piyn;
				PRECISION piz3 = zt * pitn  - t2zn * pinn;	// the dimensions of piu3, piz3 aren't consistent with the others...

				PRECISION denom_pi = xi0 * rho_max * pi_mag;

				PRECISION a0 = pi_mag / rho_max / T_aniso_mag;
				PRECISION a1 = fabs(Trpi / denom_pi);
				PRECISION a2 = fabs(piu0 / denom_pi);
				PRECISION a3 = fabs(piu1 / denom_pi);
				PRECISION a4 = fabs(piu2 / denom_pi);
				PRECISION a5 = fabs(piu3 / denom_pi);
				PRECISION a6 = fabs(piz0 / denom_pi);
				PRECISION a7 = fabs(piz1 / denom_pi);
				PRECISION a8 = fabs(piz2 / denom_pi);
				PRECISION a9 = fabs(piz3 / denom_pi);

				// compute the rho factor
				PRECISION rho_pi = fmax(a0, fmax(a1, fmax(a2, fmax(a3, fmax(a4, fmax(a5, fmax(a6, fmax(a7, fmax(a8, a9)))))))));

				PRECISION factor_pi;
				if(rho_pi > eps) factor_pi = tanh(rho_pi) / rho_pi;
				else factor_pi = 1.0;

				// regulate
				Q_current->pitt[s] *= factor_pi;
				Q_current->pitx[s] *= factor_pi;
				Q_current->pity[s] *= factor_pi;
				Q_current->pitn[s] *= factor_pi;
				Q_current->pixx[s] *= factor_pi;
				Q_current->pixy[s] *= factor_pi;
				Q_current->pixn[s] *= factor_pi;
				Q_current->piyy[s] *= factor_pi;
				Q_current->piyn[s] *= factor_pi;
				Q_current->pinn[s] *= factor_pi;
			#endif
			// regulate Wmu
			#ifdef W_TZ_MU
				PRECISION Wt = Q_current->WtTz[s];
				PRECISION Wx = Q_current->WxTz[s];
				PRECISION Wy = Q_current->WyTz[s];
				PRECISION Wn = Q_current->WnTz[s];

				// sqrt(2.W.W), W.u and W.z
				PRECISION W_mag = sqrt(fabs(Wt * Wt  -  Wx * Wx  -  Wy * Wy  -  t2 * Wn * Wn));

				PRECISION Wu = Wt * ut  -  Wx * ux  -  Wy * uy  -  Wn * t2un;
				PRECISION Wz = Wt * zt  -  Wn * t2zn;

				PRECISION denom_W = xi0 * rho_max * W_mag;

				PRECISION b0 = W_mag / rho_max / T_aniso_mag;
				PRECISION b1 = fabs(Wu / denom_W);
				PRECISION b2 = fabs(Wz / denom_W);

				// compute the rho factor
				PRECISION rho_W = fmax(b0, fmax(b1, b2));

				PRECISION factor_W;
				if(rho_W > eps) factor_W = tanh(rho_W) / rho_W;
				else factor_W = 1.0;

				// regulate
				Q_current->WtTz[s] *= factor_W;
				Q_current->WxTz[s] *= factor_W;
				Q_current->WyTz[s] *= factor_W;
				Q_current->WnTz[s] *= factor_W;
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

	// compute uS, e, p (move source code here?)
	set_inferred_variables(qS, e, p, u, uS, t, nx, ny, nz, ncx, ncy, fTSol_1);
	//set_inferred_variables_new(qS, e, p, uS, t, nx, ny, nz, ncx, ncy);

	// regulate dissipative components of qS
//#if defined PIMUNU || (defined PIMUNU && defined W_TZ_MU)
	regulate_residual_dissipative_currents(t, qS, e, p, uS, nx, ny, nz, ncx, ncy, dx, dy, dz, validityDomain);
//#endif

	// set ghost cells for qS, uS, e, p
	set_ghost_cells(qS, e, p, uS, nx, ny, nz, ncx, ncy);

	// second intermediate time step (compute Q = qS + dt.(S - dHx/dx - dHy/dy - dHz/dz))
	euler_step(t, qS, Q, e, p, uS, u, nx, ny, nz, ncx, ncy, dt, dx, dy, dz, etabar);

	// Runge-Kutta: (q + Q) / 2 -> Q
	convex_combination(q, Q, nx, ny, nz, ncx, ncy);

	// swap up and u
	swap_fluid_velocity(&up, &u);



	// this can all be grouped together in one function somehow
	// compute u, e, p
	set_inferred_variables(Q, e, p, up, u, t, nx, ny, nz, ncx, ncy, fTSol_2);
	//set_inferred_variables_new(Q, e, p, u, t, nx, ny, nz, ncx, ncy);

	// regulate dissipative components of Q
//#if defined PIMUNU || (defined PIMUNU && defined W_TZ_MU)
	regulate_residual_dissipative_currents(t, Q, e, p, u, nx, ny, nz, ncx, ncy, dx, dy, dz, validityDomain);
//#endif

	// set ghost cells for Q, u, e, p
	set_ghost_cells(Q, e, p, u, nx, ny, nz, ncx, ncy);



	// swap q and Q
	set_current_conserved_variables();
}




