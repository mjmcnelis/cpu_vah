#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include "../include/Precision.h"
#include "../include/Projections.h"
#include "../include/DynamicalVariables.h"
using namespace std;

#define REGULATION_SCHEME 1		// 0 = no residual shear stress regulation 
								// 1 = old regulation scheme 
								// 2 = new regulation scheme

#define TEST_PIMUNU 0			// 1 = test piT orthogonality and tracelessness
#define TEST_WTZMU 0			// 1 = test WTz orthogonality 

precision piu_error = 1.e-13;
precision piz_error = 1.e-13;
precision trpi_error = 1.e-13;

precision WTzu_error = 1.e-13;
precision WTzz_error = 1.e-13;

inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}

void test_pimunu_properties(precision pitt, precision pitx, precision pity, precision pitn, precision pixx, precision pixy, precision pixn, precision piyy, precision piyn, precision pinn, precision ut, precision ux, precision uy, precision un, precision zt, precision zn, precision t2)
{
	precision trpi = fabs(pitt  -  pixx  -  piyy  -  t2 * pinn);

	precision piu0 = fabs(pitt * ut  -  pitx * ux  -  pity * uy  -  t2 * pitn * un);
	precision piu1 = fabs(pitx * ut  -  pixx * ux  -  pixy * uy  -  t2 * pixn * un);
	precision piu2 = fabs(pity * ut  -  pixy * ux  -  piyy * uy  -  t2 * piyn * un);
	precision piu3 = fabs(pitn * ut  -  pixn * ux  -  piyn * uy  -  t2 * pinn * un);

	precision piz0 = fabs(zt * pitt  -  t2 * zn * pitn);
	precision piz1 = fabs(zt * pitx  -  t2 * zn * pixn);
	precision piz2 = fabs(zt * pity  -  t2 * zn * piyn);
	precision piz3 = fabs(zt * pitn  -  t2 * zn * pinn);

	precision piu = fmax(piu0, fmax(piu1, fmax(piu2, piu3)));
	precision piz = fmax(piz0, fmax(piz1, fmax(piz2, piz3)));

    if(piu > piu_error) 
    {
    	piu_error = piu;
    	printf("test_pimunu_properties error: piT is not orthogonal to u (%.6g)\n", piu);
    }
    if(piz > piz_error)  
    {
    	piz_error = piz;
    	printf("test_pimunu_properties error: piT is not orthogonal to z (%.6g)\n", piz);
    }
    if(trpi > trpi_error) 
    {
    	trpi_error = trpi;
    	printf("test_pimunu_properties error: piT is not traceless (%.6g)\n", trpi);
    }
}


void test_WTzmu_properties(precision WtTz, precision WxTz, precision WyTz, precision WnTz, precision ut, precision ux, precision uy, precision un, precision zt, precision zn, precision t2)
{
	precision WTzu = fabs(WtTz * ut  -  WxTz * ux  -  WyTz * uy  -  t2 * WnTz * un);
	precision WTzz = fabs(WtTz * zt  -  t2 * WnTz * zn);

    if(WTzu > WTzu_error) 
    {
    	WTzu_error = WTzu;
    	printf("test_WTzmu_properties error: WTz is not orthogonal to u (%.6g)\n", WTzu);
    }
    if(WTzz > WTzz_error)  
    {
    	WTzz_error = WTzz;
    	printf("test_pimunu_properties error: WTz is not orthogonal to z (%.6g)\n", WTzz);
    }
}


void regulate_dissipative_currents(precision t, CONSERVED_VARIABLES * const __restrict__ q, precision * const __restrict__ e, const FLUID_VELOCITY * const __restrict__ u, int nx, int ny, int nz)
{
	precision eps = 1.e-7;		// is there a more effective way to regulate pl, pt?
	precision xi0 = 0.1;
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

			// #if (REGULATION_SCHEME == 0)
			// 	continue;
			// #endif

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

			#if (TEST_PIMUNU == 1)
				test_pimunu_properties(pitt, pitx, pity, pitn, pixx, pixy, pixn, piyy, piyn, pinn, ut, ux, uy, un, zt, zn, t2);
			#endif

				precision rho_pi;

			// old regulation scheme: regulate all components the same factor (either if too viscous, not orthogonal or not traceless)
			#if (REGULATION_SCHEME == 1)
				precision pi_mag = sqrt(fabs(pitt * pitt  +  pixx * pixx  +  piyy * piyy  +  t4 * pinn * pinn  -  2.0 * (pitx * pitx  +  pity * pity  -  pixy * pixy  +  t2 * (pitn * pitn  -  pixn * pixn  -  piyn * piyn))));

				precision trpi = fabs(pitt  -  pixx  -  piyy -  t2 * pinn);

				precision piu0 = fabs(pitt * ut  -  pitx * ux  -  pity * uy  -  t2 * pitn * un);
				precision piu1 = fabs(pitx * ut  -  pixx * ux  -  pixy * uy  -  t2 * pixn * un);
				precision piu2 = fabs(pity * ut  -  pixy * ux  -  piyy * uy  -  t2 * piyn * un);
				precision piu3 = fabs(pitn * ut  -  pixn * ux  -  piyn * uy  -  t2 * pinn * un);

				precision piz0 = fabs(zt * pitt  -  t2 * zn * pitn);
				precision piz1 = fabs(zt * pitx  -  t2 * zn * pixn);
				precision piz2 = fabs(zt * pity  -  t2 * zn * piyn);
				precision piz3 = fabs(zt * pitn  -  t2 * zn * pinn);	// dimensions of piu3, piz3 aren't consistent with the others...

				precision denom_pi = xi0 * rho_max * pi_mag;

				precision a0 = pi_mag / (rho_max * T_aniso_mag);
				precision a1 = trpi / denom_pi;
				precision a2 = piu0 / denom_pi;
				precision a3 = piu1 / denom_pi;
				precision a4 = piu2 / denom_pi;
				precision a5 = piu3 / denom_pi;
				precision a6 = piz0 / denom_pi;
				precision a7 = piz1 / denom_pi;
				precision a8 = piz2 / denom_pi;
				precision a9 = piz3 / denom_pi;

				rho_pi = fmax(a0, fmax(a1, fmax(a2, fmax(a3, fmax(a4, fmax(a5, fmax(a6, fmax(a7, fmax(a8, a9)))))))));

			// new regulation scheme: reproject and regulate pi components to enforce tracelessness / orthogonality
			#elif (REGULATION_SCHEME == 2)

				double_transverse_projection Xi2(Xi, t2, t4);

				Xi2.double_transverse_project_tensor(pitt, pitx, pity, pitn, pixx, pixy, pixn, piyy, piyn, pinn);

				precision pi_mag = sqrt(fabs(pitt * pitt  +  pixx * pixx  +  piyy * piyy  +  t4 * pinn * pinn  -  2.0 * (pitx * pitx  +  pity * pity  -  pixy * pixy  +  t2 * (pitn * pitn  -  pixn * pixn  -  piyn * piyn))));

				rho_pi = pi_mag / (rho_max * T_aniso_mag);
			#endif

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

			#if (TEST_WTZMU == 1)
				test_WTzmu_properties(WtTz, WxTz, WyTz, WnTz, ut, ux, uy, un, zt, zn, t2);
			#endif

				precision rho_W;

			// old regulation scheme: regulate all components the same factor (either if too viscous or not orthogonal)
			#if (REGULATION_SCHEME == 1)

				precision WTz_mag = sqrt(fabs(WtTz * WtTz  -  WxTz * WxTz  -  WyTz * WyTz  -  t2 * WnTz * WnTz));

				precision WTzu = fabs(WtTz * ut  -  WxTz * ux  -  WyTz * uy  -  t2 * WnTz * un);
				precision WTzz = fabs(WtTz * zt  -  t2 * WnTz * zn);

				precision denom_W = xi0 * rho_max * W_mag;

				precision b0 = WTz_mag / (rho_max * T_aniso_mag);
				precision b1 = WTzu / denom_W;
				precision b2 = WTzz / denom_W;

				rho_W = fmax(b0, fmax(b1, b2));

			// new regulation scheme: reproject and regulate WTz components to enforce orthogonality
			#elif (REGULATION_SCHEME == 2)

				Xi.transverse_project_vector(WtTz, WxTz, WyTz, WnTz);

				precision WTz_mag = sqrt(2.0 * fabs(WtTz * WtTz  -  WxTz * WxTz  -  WyTz * WyTz  -  t2 * WnTz * WnTz));

				rho_W = WTz_mag / (rho_max * T_aniso_mag);
			#endif

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




