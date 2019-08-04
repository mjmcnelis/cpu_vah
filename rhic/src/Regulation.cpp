#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"
#include "../include/Projections.h"
using namespace std;

#define XI0 	0.1						// regulation parameters
#define RHO_MAX 5.0
#define REGULATION_SCHEME 1				// 1 = new regulation scheme (old regulation otherwise)


inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


void regulate_residual_currents(precision t, hydro_variables * const __restrict__ q, precision * const __restrict__ e, const fluid_velocity * const __restrict__ u, int nx, int ny, int nz)
{
#ifdef ANISO_HYDRO
	precision eps = 1.e-6;	// enforce pl, pt to be positive

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
				precision pl  = q[s].pl;

				if(pl < eps)
				{
					q[s].pl = eps;
					pl = eps;
				}
			#if (PT_MATCHING == 1)
				precision pt  = q[s].pt;
				if(pt < eps)
				{
					q[s].pt = eps;
					pt = eps;
				}
			#endif
			#if (NUMBER_OF_RESIDUAL_CURRENTS != 0)
				precision ux = u[s].ux;
				precision uy = u[s].uy;
				precision un = u[s].un;
				precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t2 * un * un);

				precision utperp = sqrt(1.  +  ux * ux  +  uy * uy);
				precision zt = t * un / utperp;
				precision zn = ut / t / utperp;
			#ifdef CONFORMAL_EOS
				precision pt = (e_s - pl) / 2.;
			#endif
				precision T_aniso_mag = sqrt(e_s * e_s  +  pl * pl  +  2. * pt * pt);
			#endif

			#ifdef PIMUNU
				precision pitt = q[s].pitt;
				precision pitx = q[s].pitx;
				precision pity = q[s].pity;
				precision pitn = q[s].pitn;
				precision pixx = q[s].pixx;
				precision pixy = q[s].pixy;
				precision pixn = q[s].pixn;
				precision piyy = q[s].piyy;
				precision piyn = q[s].piyn;
				precision pinn = q[s].pinn;

				precision pi_mag = sqrt(fabs(pitt * pitt  +  pixx * pixx  +  piyy * piyy  +  t4 * pinn * pinn  -  2. * (pitx * pitx  +  pity * pity  -  pixy * pixy  +  t2 * (pitn * pitn  -  pixn * pixn  -  piyn * piyn))));

				precision trpi = fabs(pitt  -  pixx  -  piyy  -  t2 * pinn);

			#if (REGULATION_SCHEME == 1)
				precision piu0 = fabs(pitt * ut  -  pitx * ux  -  pity * uy  -  t2 * pitn * un) / ut;
				precision piu1 = fabs(pitx * ut  -  pixx * ux  -  pixy * uy  -  t2 * pixn * un) / ut;
				precision piu2 = fabs(pity * ut  -  pixy * ux  -  piyy * uy  -  t2 * piyn * un) / ut;
				precision piu3 = fabs(pitn * ut  -  pixn * ux  -  piyn * uy  -  t2 * pinn * un) * t / ut;

				precision piz0 = fabs(zt * pitt  -  t2 * zn * pitn) / (t * zn);
				precision piz1 = fabs(zt * pitx  -  t2 * zn * pixn) / (t * zn);
				precision piz2 = fabs(zt * pity  -  t2 * zn * piyn) / (t * zn);
				precision piz3 = fabs(zt * pitn  -  t2 * zn * pinn) / zn;
			#else
				precision piu0 = fabs(pitt * ut  -  pitx * ux  -  pity * uy  -  t2 * pitn * un);
				precision piu1 = fabs(pitx * ut  -  pixx * ux  -  pixy * uy  -  t2 * pixn * un);
				precision piu2 = fabs(pity * ut  -  pixy * ux  -  piyy * uy  -  t2 * piyn * un);
				precision piu3 = fabs(pitn * ut  -  pixn * ux  -  piyn * uy  -  t2 * pinn * un) * t;

				precision piz0 = fabs(zt * pitt  -  t2 * zn * pitn);
				precision piz1 = fabs(zt * pitx  -  t2 * zn * pixn);
				precision piz2 = fabs(zt * pity  -  t2 * zn * piyn);
				precision piz3 = fabs(zt * pitn  -  t2 * zn * pinn) * t;
			#endif

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

				if(factor_pi > 1.)
				{
					printf("Regulation error: factor_pi = %.6g\n", factor_pi);
				}

				q[s].pitt = factor_pi * pitt;
				q[s].pitx = factor_pi * pitx;
				q[s].pity = factor_pi * pity;
				q[s].pitn = factor_pi * pitn;
				q[s].pixx = factor_pi * pixx;
				q[s].pixy = factor_pi * pixy;
				q[s].pixn = factor_pi * pixn;
				q[s].piyy = factor_pi * piyy;
				q[s].piyn = factor_pi * piyn;
				q[s].pinn = factor_pi * pinn;
			#endif

			#ifdef WTZMU
				precision WtTz = q[s].WtTz;
				precision WxTz = q[s].WxTz;
				precision WyTz = q[s].WyTz;
				precision WnTz = q[s].WnTz;

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

				q[s].WtTz = factor_W * WtTz;
				q[s].WxTz = factor_W * WxTz;
				q[s].WyTz = factor_W * WyTz;
				q[s].WnTz = factor_W * WnTz;
			#endif
			}
		}
	}
#endif
}




void regulate_viscous_currents(precision t, hydro_variables * const __restrict__ q, precision * const __restrict__ e, const fluid_velocity * const __restrict__ u, int nx, int ny, int nz)
{
	// may change to adopt regulation of ttt, ttx, etc
#if (NUMBER_OF_VISCOUS_CURRENTS != 0)
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
				precision p = equilibriumPressure(e_s);
				precision ux = u[s].ux;
				precision uy = u[s].uy;
				precision un = u[s].un;
				precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t2 * un * un);

				precision Teq_mag = sqrt(e_s * e_s  +  3. * p * p);

			#ifdef PIMUNU
				precision pitt = q[s].pitt;
				precision pitx = q[s].pitx;
				precision pity = q[s].pity;
				precision pitn = q[s].pitn;
				precision pixx = q[s].pixx;
				precision pixy = q[s].pixy;
				precision pixn = q[s].pixn;
				precision piyy = q[s].piyy;
				precision piyn = q[s].piyn;
				precision pinn = q[s].pinn;

				precision pi_mag = sqrt(fabs(pitt * pitt  +  pixx * pixx  +  piyy * piyy  +  t4 * pinn * pinn  -  2. * (pitx * pitx  +  pity * pity  -  pixy * pixy  +  t2 * (pitn * pitn  -  pixn * pixn  -  piyn * piyn))));

				precision trpi = fabs(pitt  -  pixx  -  piyy  -  t2 * pinn);

			#if (REGULATION_SCHEME == 1)
				precision piu0 = fabs(pitt * ut  -  pitx * ux  -  pity * uy  -  t2 * pitn * un) / ut;
				precision piu1 = fabs(pitx * ut  -  pixx * ux  -  pixy * uy  -  t2 * pixn * un) / ut;
				precision piu2 = fabs(pity * ut  -  pixy * ux  -  piyy * uy  -  t2 * piyn * un) / ut;
				precision piu3 = fabs(pitn * ut  -  pixn * ux  -  piyn * uy  -  t2 * pinn * un) * t / ut;
			#else
				precision piu0 = fabs(pitt * ut  -  pitx * ux  -  pity * uy  -  t2 * pitn * un);
				precision piu1 = fabs(pitx * ut  -  pixx * ux  -  pixy * uy  -  t2 * pixn * un);
				precision piu2 = fabs(pity * ut  -  pixy * ux  -  piyy * uy  -  t2 * piyn * un);
				precision piu3 = fabs(pitn * ut  -  pixn * ux  -  piyn * uy  -  t2 * pinn * un) * t;
			#endif

				precision denom_pi = xi0 * rho_max * pi_mag;

				precision a0 = pi_mag / (rho_max * Teq_mag);
				precision a1 = trpi / denom_pi;
				precision a2 = piu0 / denom_pi;
				precision a3 = piu1 / denom_pi;
				precision a4 = piu2 / denom_pi;
				precision a5 = piu3 / denom_pi;

				precision rho_pi = fmax(a0, fmax(a1, fmax(a2, fmax(a3, fmax(a4, a5)))));

				precision factor_pi;
				if(rho_pi > 1.e-5) factor_pi = tanh(rho_pi) / rho_pi;
				else factor_pi = 1.  -  rho_pi * rho_pi / 3.;

				q[s].pitt = factor_pi * pitt;
				q[s].pitx = factor_pi * pitx;
				q[s].pity = factor_pi * pity;
				q[s].pitn = factor_pi * pitn;
				q[s].pixx = factor_pi * pixx;
				q[s].pixy = factor_pi * pixy;
				q[s].pixn = factor_pi * pixn;
				q[s].piyy = factor_pi * piyy;
				q[s].piyn = factor_pi * piyn;
				q[s].pinn = factor_pi * pinn;
			#endif

			#ifdef PI
				precision Pi = q[s].Pi;
				precision rho_bulk = fabs(Pi / (rho_max * Teq_mag));

				precision factor_bulk;
				if(rho_bulk > 1.e-5) factor_bulk = tanh(rho_bulk) / rho_bulk;
				else factor_bulk = 1.  -  rho_bulk * rho_bulk / 3.;	// 2nd-order expansion

				q[s].Pi = factor_bulk * Pi;
			#endif
			}
		}
	}
#endif
}




