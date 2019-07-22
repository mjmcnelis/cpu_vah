
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include "../include/Regulation.h"
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"

using namespace std;

int conformal_error = 0;
double bulkPi_error = 1.e-14;

inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}

void regulate_dissipative_currents(precision t, CONSERVED_VARIABLES * const __restrict__ q, precision * const __restrict__ e, const FLUID_VELOCITY * const __restrict__ u, int nx, int ny, int nz)
{
	precision eps = 1.e-8;	// is there a more effective way to regulate pl, pt?

	//precision xi0 = 0.1;		// regulation parameters (maybe move these in hydro parameters?)
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

			// regulate bulk pressure if conformal
			#ifdef CONFORMAL_EOS
				#if (PT_MATCHING == 1)
					precision peq = equilibriumPressure(e_s);
					precision pt_s  = q->pt[s];

					q->pl[s] = peq  +  2./3. * (pl - pt_s);
					q->pt[s] = peq  -  1./3. * (pl - pt_s);
				#endif
			#endif

				if(pl < eps)
				{
					printf("Regulation: pl = %lf is regulated to %lf", pl, eps);
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
			#else
				precision pt = 0.5 * (e_s - pl);
			#endif

				precision ut = u->ut[s];
				precision ux = u->ux[s];
				precision uy = u->uy[s];
				precision un = u->un[s];

				precision utperp = sqrt(1.0  +  ux * ux  +  uy * uy);
				precision zt = t * un / utperp;
				precision zn = ut / t / utperp;

				precision t2un = t2 * un;
				precision t2zn = t2 * zn;

				// sqrt(T_aniso.T_aniso)
				precision T_aniso_mag = sqrt(e_s * e_s  +  pl * pl  +  2.0 * pt * pt);

				// transverse projector (for new regulation scheme)
			#if (NUMBER_OF_RESIDUAL_CURRENTS != 0)
				precision Xitt = 1.0  -  ut * ut  +  zt * zt;
				precision Xitx = - ut * ux;
				precision Xity = - ut * uy;
				precision Xitn = - ut * un  +  zt * zn;
				precision Xixx = - 1.0  -  ux * ux;
				precision Xixy = - ux * uy;
				precision Xixn = - ux * un;
				precision Xiyy = - 1.0  -  uy * uy;
				precision Xiyn = - uy * un;
				precision Xinn = - 1.0 / t2  -  un * un  +  zn * zn;
			#endif


				// regulate transverse shear stress
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

				// sqrt(pi.pi), Tr(pi), u.pi and z.pi
				/*
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
				if(rho_pi > 1.e-5) factor_pi = tanh(rho_pi) / rho_pi;
				else factor_pi = 1.0  -  rho_pi * rho_pi / 3.;

				// regulate
				q->pitt[s] *= factor_pi;
				q->pitx[s] *= factor_pi;
				q->pity[s] *= factor_pi;
				q->pitn[s] *= factor_pi;
				q->pixx[s] *= factor_pi;
				q->pixy[s] *= factor_pi;
				q->pixn[s] *= factor_pi;
				q->piyy[s] *= factor_pi;
				q->piyn[s] *= factor_pi;
				q->pinn[s] *= factor_pi;
				*/

				// new regulation scheme
				precision pi_mag = sqrt(fabs(pitt * pitt  +  pixx * pixx  +  piyy * piyy  +  t2 * t2 * pinn * pinn  -  2.0 * (pitx * pitx  +  pity * pity  -  pixy * pixy  +  t2 * (pitn * pitn  -  pixn * pixn  -  piyn * piyn))));

				precision rho_pi = pi_mag / rho_max / T_aniso_mag;

				precision factor_pi;
				if(rho_pi > 1.e-5)
				{
					factor_pi = tanh(rho_pi) / rho_pi;
				}
				else
				{
					factor_pi = 1.0  -  rho_pi * rho_pi / 3.;
				}

				// Xi^{\mu\nu}_{\alpha\beta} (can invoke some symmetries)
				// precision Xitt_tt =  0.5 * Xitt * Xitt;
				// precision Xitt_tx = -0.5 * Xitt * Xitx;
				// precision Xitt_ty = -0.5 * Xitt * Xity;
				// precision Xitt_tn = -0.5 * Xitt * Xitn;

				// precision Xitt_xx =  Xitx * Xitx  -  0.5 * Xitt * Xixx;
				// precision Xitt_xy =  Xitx * Xity  -  0.5 * Xitt * Xixy;
				// precision Xitt_xn =  t2 * (Xitx * Xitn  -  0.5 * Xitt * Xixn);


				// precision Xixx_xx = 0.5 * Xixx * Xixx;
				// precision Xiyy_yy = 0.5 * Xiyy * Xiyy;

				// precision Xinn_nn = 0.5 * t4 * Xinn * Xinn;


				precision pitt_pro = factor_pi * (Xitt_tt * pitt  +  Xitt_xx * pixx  +  Xitt_yy * piyy  +  Xitt_nn * pinn  +  2.0 * (Xitt_tx * pitx  +  Xitt_ty * pity  +  Xitt_tn * pitn  +  Xitt_xy * pixy  +  Xitt_xn * pixn  +  Xitt_yn * piyn));
				precision pitx_pro = factor_pi * (Xitx_tt * pitt  +  Xitx_xx * pixx  +  Xitx_yy * piyy  +  Xitx_nn * pinn  +  2.0 * (Xitx_tx * pitx  +  Xitx_ty * pity  +  Xitx_tn * pitn  +  Xitx_xy * pixy  +  Xitx_xn * pixn  +  Xitx_yn * piyn));
				precision pity_pro = factor_pi * (Xity_tt * pitt  +  Xity_xx * pixx  +  Xity_yy * piyy  +  Xity_nn * pinn  +  2.0 * (Xity_tx * pitx  +  Xity_ty * pity  +  Xity_tn * pitn  +  Xity_xy * pixy  +  Xity_xn * pixn  +  Xity_yn * piyn));
				precision pitn_pro = factor_pi * (Xitn_tt * pitt  +  Xitn_xx * pixx  +  Xitn_yy * piyy  +  Xitn_nn * pinn  +  2.0 * (Xitn_tx * pitx  +  Xitn_ty * pity  +  Xitn_tn * pitn  +  Xitn_xy * pixy  +  Xitn_xn * pixn  +  Xitn_yn * piyn));
				precision pixx_pro = factor_pi * (Xixx_tt * pitt  +  Xixx_xx * pixx  +  Xixx_yy * piyy  +  Xixx_nn * pinn  +  2.0 * (Xixx_tx * pitx  +  Xixx_ty * pity  +  Xixx_tn * pitn  +  Xixx_xy * pixy  +  Xixx_xn * pixn  +  Xixx_yn * piyn));
				precision pixy_pro = factor_pi * (Xixy_tt * pitt  +  Xixy_xx * pixx  +  Xixy_yy * piyy  +  Xixy_nn * pinn  +  2.0 * (Xixy_tx * pitx  +  Xixy_ty * pity  +  Xixy_tn * pitn  +  Xixy_xy * pixy  +  Xixy_xn * pixn  +  Xixy_yn * piyn));
				precision pixn_pro = factor_pi * (Xixn_tt * pitt  +  Xixn_xx * pixx  +  Xixn_yy * piyy  +  Xixn_nn * pinn  +  2.0 * (Xixn_tx * pitx  +  Xixn_ty * pity  +  Xixn_tn * pitn  +  Xixn_xy * pixy  +  Xixn_xn * pixn  +  Xixn_yn * piyn));
				precision piyy_pro = factor_pi * (Xiyy_tt * pitt  +  Xiyy_xx * pixx  +  Xiyy_yy * piyy  +  Xiyy_nn * pinn  +  2.0 * (Xiyy_tx * pitx  +  Xiyy_ty * pity  +  Xiyy_tn * pitn  +  Xiyy_xy * pixy  +  Xiyy_xn * pixn  +  Xiyy_yn * piyn));
				precision piyn_pro = factor_pi * (Xiyn_tt * pitt  +  Xiyn_xx * pixx  +  Xiyn_yy * piyy  +  Xiyn_nn * pinn  +  2.0 * (Xiyn_tx * pitx  +  Xiyn_ty * pity  +  Xiyn_tn * pitn  +  Xiyn_xy * pixy  +  Xiyn_xn * pixn  +  Xiyn_yn * piyn));
				precision pinn_pro = factor_pi * (Xinn_tt * pitt  +  Xinn_xx * pixx  +  Xinn_yy * piyy  +  Xinn_nn * pinn  +  2.0 * (Xinn_tx * pitx  +  Xinn_ty * pity  +  Xinn_tn * pitn  +  Xinn_xy * pixy  +  Xinn_xn * pixn  +  Xinn_yn * piyn));

				q->pitt[s] = pitt_pro;
				q->pitx[s] = pitx_pro;
				q->pity[s] = pity_pro;
				q->pitn[s] = pitn_pro;
				q->pixx[s] = pixx_pro;
				q->pixy[s] = pixy_pro;
				q->pixn[s] = pixn_pro;
				q->piyy[s] = piyy_pro;
				q->piyn[s] = piyn_pro;
				q->pinn[s] = pinn_pro;

			#endif

			// regulate Wmu
			#ifdef WTZMU
				precision WtTz = q->WtTz[s];
				precision WxTz = q->WxTz[s];
				precision WyTz = q->WyTz[s];
				precision WnTz = q->WnTz[s];

				precision W_mag = sqrt(2.0 * fabs(WtTz * WtTz  -  WxTz * WxTz  -  WyTz * WyTz  -  t2 * WnTz * WnTz));

				precision rho_W = W_mag / rho_max / T_aniso_mag;

				precision factor_W;
				if(rho_W > 1.e-5)
				{
					factor_W = tanh(rho_W) / rho_W;
				}
				else
				{
					factor_W = 1.0  -  rho_W * rho_W / 3.;	// 2nd-order expansion
				}

				// regulate and reproject WTz components to ensure orthogonality to u, z (kill 2 birds with 1 stone)
				precision WtTz_pro = factor_W * (Xitt * WtTz  -  Xitx * WxTz  -  Xity * WyTz  -  t2 * Xitn * WnTz);
				precision WxTz_pro = factor_W * (Xitx * WtTz  -  Xixx * WxTz  -  Xixy * WyTz  -  t2 * Xixn * WnTz);
				precision WyTz_pro = factor_W * (Xity * WtTz  -  Xixy * WxTz  -  Xiyy * WyTz  -  t2 * Xiyn * WnTz);
				precision WnTz_pro = factor_W * (Xitn * WtTz  -  Xixn * WxTz  -  Xiyn * WyTz  -  t2 * Xinn * WnTz);

				q->WtTz[s] = WtTz_pro;
				q->WxTz[s] = WxTz_pro;
				q->WyTz[s] = WyTz_pro;
				q->WnTz[s] = WnTz_pro;

			#endif

			}
		}
	}
}




