
#include <math.h>

#include "../include/Regulation.h"
#include "../include/DynamicalVariables.h"


inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}

void regulate_dissipative_currents(PRECISION t, CONSERVED_VARIABLES * const __restrict__ Q_current, PRECISION * const __restrict__ e, const FLUID_VELOCITY * const __restrict__ u, int nx, int ny, int nz)
{
	PRECISION eps = 1.e-7;

	PRECISION xi0 = 0.1;		// regulation parameters (maybe move these in hydro parameters?)
	PRECISION rho_max = 1.0;

	PRECISION t2 = t * t;

	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				PRECISION e_s = e[s];
				PRECISION pl = Q_current->pl[s];

				// is this necessary?
				if(e_s < 0.0) e[s] = eps;
				
				if(pl < 0.0) Q_current->pl[s] = eps;

				PRECISION pt = 0.5 * (e_s - pl);	// temporary

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