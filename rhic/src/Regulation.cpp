
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include "../include/Regulation.h"
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"

using namespace std;

inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}

void regulate_dissipative_currents(precision t, CONSERVED_VARIABLES * const __restrict__ ql, precision * const __restrict__ el, const FLUID_VELOCITY * const __restrict__ ul, int nx, int ny, int nz)
{
	precision epsilon = 1.e-7;

	precision xi0 = 0.1;		// regulation parameters (maybe move these in hydro parameters?)
	precision rho_max = 1.0;

	precision t2 = t * t;

	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				precision e_s = el[s];
				precision pl  = ql->pl[s];
				precision pt  = ql->pt[s];

				// if(e_s < eps) // already regulated in inferred variables
				// {
				// 	el[s] = eps;
				// 	e_s = eps;
				// }
				if(pl < epsilon || std::isnan(pl))
				{
					ql->pl[s] = epsilon;
					pl = epsilon;
				}
			// #ifdef CONFORMAL_EOS
			// 	ql->pt[s] = 0.5 * (e_s - pl);
			// 	pt = 0.5 * (e_s - pl);
			// #endif
				if(pt < epsilon || std::isnan(pt))
				{
					ql->pt[s] = epsilon;
					pt = epsilon;
				}


				precision ut = ul->ut[s];
				precision ux = ul->ux[s];
				precision uy = ul->uy[s];
				precision un = ul->un[s];

				precision utperp = sqrt(1.0  +  ux * ux  +  uy * uy);
				precision zt = t * un / utperp;
				precision zn = ut / t / utperp;

				precision t2un = t2 * un;
				precision t2zn = t2 * zn;

				// sqrt(T_aniso.T_aniso)
				precision T_aniso_mag = sqrt(e_s * e_s  +  pl * pl  +  2.0 * pt * pt);

				// regulate transverse shear stress
			#ifdef PIMUNU
				precision pitt = ql->pitt[s];
				precision pitx = ql->pitx[s];
				precision pity = ql->pity[s];
				precision pitn = ql->pitn[s];
				precision pixx = ql->pixx[s];
				precision pixy = ql->pixy[s];
				precision pixn = ql->pixn[s];
				precision piyy = ql->piyy[s];
				precision piyn = ql->piyn[s];
				precision pinn = ql->pinn[s];


				// should I enforce tracelessness and orthogonality first
				// then regulate with a0?

				// sqrt(pi.pi), Tr(pi), u.pi and z.pi
				precision pi_mag = sqrt(fabs(pitt * pitt  +  pixx * pixx  +  piyy * piyy  +  t2 * t2 * pinn * pinn  -  2.0 * (pitx * pitx  +  pity * pity  -  pixy * pixy  +  t2 * (pitn * pitn  -  pixn * pixn  -  piyn * piyn))));

				precision Trpi = pitt  -  pixx - piyy -  t2 * pinn;

				precision piu0 = pitt * ut  -  pitx * ux  -  pity * uy  -  pitn * t2un;
				precision piu1 = pitx * ut  -  pixx * ux  -  pixy * uy  -  pixn * t2un;
				precision piu2 = pity * ut  -  pixy * ux  -  piyy * uy  -  piyn * t2un;
				precision piu3 = pitn * ut  -  pixn * ux  -  piyn * uy  -  pinn * t2un;

				precision piz0 = zt * pitt  - t2zn * pitn;
				precision piz1 = zt * pitx  - t2zn * pixn;
				precision piz2 = zt * pity  - t2zn * piyn;
				precision piz3 = zt * pitn  - t2zn * pinn;	// the dimensions of piu3, piz3 aren't consistent with the others...

				precision denom_pi = xi0 * rho_max * pi_mag;

				precision a0 = pi_mag / rho_max / T_aniso_mag;
				precision a1 = fabs(Trpi / denom_pi);
				precision a2 = fabs(piu0 / denom_pi);
				precision a3 = fabs(piu1 / denom_pi);
				precision a4 = fabs(piu2 / denom_pi);
				precision a5 = fabs(piu3 / denom_pi);
				precision a6 = fabs(piz0 / denom_pi);
				precision a7 = fabs(piz1 / denom_pi);
				precision a8 = fabs(piz2 / denom_pi);
				precision a9 = fabs(piz3 / denom_pi);

				// compute the rho factor
				precision rho_pi = fmax(a0, fmax(a1, fmax(a2, fmax(a3, fmax(a4, fmax(a5, fmax(a6, fmax(a7, fmax(a8, a9)))))))));

				precision factor_pi;
				if(rho_pi > eps) factor_pi = tanh(rho_pi) / rho_pi;
				else factor_pi = 1.0;

				// regulate
				ql->pitt[s] *= factor_pi;
				ql->pitx[s] *= factor_pi;
				ql->pity[s] *= factor_pi;
				ql->pitn[s] *= factor_pi;
				ql->pixx[s] *= factor_pi;
				ql->pixy[s] *= factor_pi;
				ql->pixn[s] *= factor_pi;
				ql->piyy[s] *= factor_pi;
				ql->piyn[s] *= factor_pi;
				ql->pinn[s] *= factor_pi;
			#endif
			// regulate Wmu
			#ifdef WTZMU
				precision WtTz = ql->WtTz[s];
				precision WxTz = ql->WxTz[s];
				precision WyTz = ql->WyTz[s];
				precision WnTz = ql->WnTz[s];

				// sqrt(2.W.W), W.u and W.z
				precision W_mag = sqrt(2.0 * fabs(WtTz * WtTz  -  WxTz * WxTz  -  WyTz * WyTz  -  t2 * WnTz * WnTz));

				precision Wu = WtTz * ut  -  WxTz * ux  -  WyTz * uy  -  WnTz * t2un;
				precision Wz = WtTz * zt  -  WnTz * t2zn;

				precision denom_W = xi0 * rho_max * W_mag;

				precision b0 = W_mag / rho_max / T_aniso_mag;
				precision b1 = fabs(Wu / denom_W);
				precision b2 = fabs(Wz / denom_W);

				// compute the rho factor
				precision rho_W = fmax(b0, fmax(b1, b2));

				precision factor_W;
				if(rho_W > eps) factor_W = tanh(rho_W) / rho_W;
				else factor_W = 1.0;

				// regulate
				ql->WtTz[s] *= factor_W;
				ql->WxTz[s] *= factor_W;
				ql->WyTz[s] *= factor_W;
				ql->WnTz[s] *= factor_W;
			#endif
			}
		}
	}
}