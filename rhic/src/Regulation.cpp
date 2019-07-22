
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include "../include/Regulation.h"
#include "../include/Precision.h"
#include "../include/Projections.h"
#include "../include/DynamicalVariables.h"

using namespace std;

inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}

void regulate_dissipative_currents(precision t, CONSERVED_VARIABLES * const __restrict__ q, precision * const __restrict__ e, const FLUID_VELOCITY * const __restrict__ u, int nx, int ny, int nz)
{
	precision eps = 1.e-7;		// is there a more effective way to regulate pl, pt?
	precision rho_max = 1.0;

	precision t2 = t * t;
	precision t4 = t2 * t2;

	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				precision e_s = e[s];
				precision pl  = q->pl[s];

				if(pl < eps)
				{
					q->pl[s] = eps;
					pl = eps;
				}

			#if (PT_MATCHING == 1)
				precision pt  = q->pt[s];
				if(pt < eps)
				{
					q->pt[s] = eps;
					pt = eps;
				}
			#endif

			#if (NUMBER_OF_RESIDUAL_CURRENTS != 0)
				precision ut = u->ut[s];
				precision ux = u->ux[s];
				precision uy = u->uy[s];
				precision un = u->un[s];

				precision utperp = sqrt(1.0  +  ux * ux  +  uy * uy);
				precision zt = t * un / utperp;
				precision zn = ut / t / utperp;

			#ifdef CONFORMAL_EOS
				precision pt = 0.5 * (e_s - pl);
			#endif

				precision T_aniso_mag = sqrt(e_s * e_s  +  pl * pl  +  2.0 * pt * pt);

				transverse_projection Xi(ut, ux, uy, un, zt, zn, t2);
			#endif

			// reproject and regulate pi components to ensure tracelessness and orthogonality to u, z (does the order matter?)
			#ifdef PIMUNU
				precision pitt = q->pitt[s];
				precision pitx = q->pitx[s];
				precision pity = q->pity[s];
				precision pitn = q->pitn[s];
				precision pixx = q->pixx[s];
				precision pixy = q->pixy[s];
				precision pixn = q->pixn[s];
				precision piyy = q->piyy[s];
				precision piyn = q->piyn[s];
				precision pinn = q->pinn[s];

				double_transverse_projection Xi2(Xi, t2, t4);

				Xi2.double_transverse_project_tensor(pitt, pitx, pity, pitn, pixx, pixy, pixn, piyy, piyn, pinn);

				precision pi_mag = sqrt(fabs(pitt * pitt  +  pixx * pixx  +  piyy * piyy  +  t4 * pinn * pinn  -  2.0 * (pitx * pitx  +  pity * pity  -  pixy * pixy  +  t2 * (pitn * pitn  -  pixn * pixn  -  piyn * piyn))));

				precision rho_pi = pi_mag / (rho_max * T_aniso_mag);

				precision factor_pi;
				if(rho_pi > 1.e-5) factor_pi = tanh(rho_pi) / rho_pi;
				else factor_pi = 1.0  -  rho_pi * rho_pi / 3.;
				
				q->pitt[s] = factor_pi * pitt;
				q->pitx[s] = factor_pi * pitx;
				q->pity[s] = factor_pi * pity;
				q->pitn[s] = factor_pi * pitn;
				q->pixx[s] = factor_pi * pixx;
				q->pixy[s] = factor_pi * pixy;
				q->pixn[s] = factor_pi * pixn;
				q->piyy[s] = factor_pi * piyy;
				q->piyn[s] = factor_pi * piyn;
				q->pinn[s] = factor_pi * pinn;
			#endif

			// reproject and regulate WTz components to ensure orthogonality to u, z (does the order matter?)
			#ifdef WTZMU
				precision WtTz = q->WtTz[s];
				precision WxTz = q->WxTz[s];
				precision WyTz = q->WyTz[s];
				precision WnTz = q->WnTz[s];

				Xi.transverse_project_vector(WtTz, WxTz, WyTz, WnTz);

				precision W_mag = sqrt(2.0 * fabs(WtTz * WtTz  -  WxTz * WxTz  -  WyTz * WyTz  -  t2 * WnTz * WnTz));

				precision rho_W = W_mag / (rho_max * T_aniso_mag);

				precision factor_W;
				if(rho_W > 1.e-5) factor_W = tanh(rho_W) / rho_W;
				else factor_W = 1.0  -  rho_W * rho_W / 3.;	// 2nd-order expansion				
				
				q->WtTz[s] = factor_W * WtTz;
				q->WxTz[s] = factor_W * WxTz;
				q->WyTz[s] = factor_W * WyTz;
				q->WnTz[s] = factor_W * WnTz;
			#endif
			}
		}
	}
}




