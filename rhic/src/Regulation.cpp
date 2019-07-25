#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"
using namespace std;

#define TEST_PIMUNU 0					// 1 to test piT orthogonality and tracelessness
#define TEST_WTZMU 0					// 1 to test WTz orthogonality

#define XI0 	0.1						// regulation parameters
#define RHO_MAX 5.0

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

	precision piz0 = fabs(pitt * zt  -  t2 * pitn * zn);
	precision piz1 = fabs(pitx * zt  -  t2 * pixn * zn);
	precision piz2 = fabs(pity * zt  -  t2 * piyn * zn);
	precision piz3 = fabs(pitn * zt  -  t2 * pinn * zn);

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

	precision xi0 = XI0;
	precision rho_max = RHO_MAX;

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

				precision utperp = sqrt(1.  +  ux * ux  +  uy * uy);
				precision zt = t * un / utperp;
				precision zn = ut / t / utperp;

			#ifdef CONFORMAL_EOS
				precision pt = (e_s - pl) / 2.;
			#endif

				precision T_aniso_mag = sqrt(e_s * e_s  +  pl * pl  +  2. * pt * pt);
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

				precision pi_mag = sqrt(fabs(pitt * pitt  +  pixx * pixx  +  piyy * piyy  +  t4 * pinn * pinn  -  2. * (pitx * pitx  +  pity * pity  -  pixy * pixy  +  t2 * (pitn * pitn  -  pixn * pixn  -  piyn * piyn))));

				precision trpi = fabs(pitt  -  pixx  -  piyy  -  t2 * pinn);

				precision piu0 = fabs(pitt * ut  -  pitx * ux  -  pity * uy  -  t2 * pitn * un) / ut;
				precision piu1 = fabs(pitx * ut  -  pixx * ux  -  pixy * uy  -  t2 * pixn * un) / ut;
				precision piu2 = fabs(pity * ut  -  pixy * ux  -  piyy * uy  -  t2 * piyn * un) / ut;
				precision piu3 = fabs(pitn * ut  -  pixn * ux  -  piyn * uy  -  t2 * pinn * un) * t / ut;

				precision piz0 = fabs(zt * pitt  -  t2 * zn * pitn) / (t * zn);
				precision piz1 = fabs(zt * pitx  -  t2 * zn * pixn) / (t * zn);
				precision piz2 = fabs(zt * pity  -  t2 * zn * piyn) / (t * zn);
				precision piz3 = fabs(zt * pitn  -  t2 * zn * pinn) / zn;	

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

				precision rho_pi = fmax(a0, fmax(a1, fmax(a2, fmax(a3, fmax(a4, fmax(a5, fmax(a6, fmax(a7, fmax(a8, a9)))))))));

				precision factor_pi;
				if(rho_pi > 1.e-5) factor_pi = tanh(rho_pi) / rho_pi;
				else factor_pi = 1.  -  rho_pi * rho_pi / 3.;

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

				precision WTz_mag = sqrt(fabs(WtTz * WtTz  -  WxTz * WxTz  -  WyTz * WyTz  -  t2 * WnTz * WnTz));

				precision WTzu = fabs(WtTz * ut  -  WxTz * ux  -  WyTz * uy  -  t2 * WnTz * un);
				precision WTzz = fabs(WtTz * zt  -  t2 * WnTz * zn);

				precision denom_W = xi0 * rho_max * WTz_mag;

				precision b0 = WTz_mag / (rho_max * T_aniso_mag);
				precision b1 = WTzu / denom_W;
				precision b2 = WTzz / denom_W;

				precision rho_W = fmax(b0, fmax(b1, b2));

				precision factor_W;
				if(rho_W > 1.e-5) factor_W = tanh(rho_W) / rho_W;
				else factor_W = 1.  -  rho_W * rho_W / 3.;	// 2nd-order expansion

				q->WtTz[s] *= factor_W;
				q->WxTz[s] *= factor_W;
				q->WyTz[s] *= factor_W;
				q->WnTz[s] *= factor_W;
			#endif
			}
		}
	}
}




