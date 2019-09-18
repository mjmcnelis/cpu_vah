#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../include/Macros.h"
#include "../include/DynamicalVariables.h"
#include "../include/InferredVariables.h"
#include "../include/Projections.h"


inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


inline precision tanh_function(precision p)
{
	return fmin(1., tanh(p) / p);
}


void regulate_residual_currents(precision t, hydro_variables * const __restrict__ q, precision * const __restrict__ e, const fluid_velocity * const __restrict__ u, lattice_parameters lattice, hydro_parameters hydro)
{
#ifdef ANISO_HYDRO
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	int regulation_scheme = hydro.regulation_scheme;
	int reprojection = hydro.reprojection;
	precision pressure_min = hydro.pressure_min;
	precision xi0 = hydro.xi0;
	precision rho_max = hydro.rho_max;

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

				if(pl < pressure_min)
				{
					q[s].pl = pressure_min;		// cutoff like energy density?
					pl = pressure_min;
				}
			#if (PT_MATCHING == 1)
				precision pt  = q[s].pt;
				if(pt < pressure_min)
				{
					q[s].pt = pressure_min;
					pt = pressure_min;
				}
			#else
				precision pt = (e_s - pl) / 2.;
			#endif

				precision ux = u[s].ux;
				precision uy = u[s].uy;
			#ifndef BOOST_INVARIANT
				precision un = u[s].un;
			#else
				precision un = 0;
			#endif
				precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t2 * un * un);

			#ifndef BOOST_INVARIANT
				precision utperp = sqrt(1.  +  ux * ux  +  uy * uy);
				precision zt = t * un / utperp;
				precision zn = ut / t / utperp;
			#else
				precision zt = 0;
				precision zn = 1. / t;
			#endif

				precision T_aniso_mag = sqrt(e_s * e_s  +  pl * pl  +  2. * pt * pt);

			#ifdef PIMUNU
				precision pitt = q[s].pitt;
				precision pitx = q[s].pitx;
				precision pity = q[s].pity;
			#ifndef BOOST_INVARIANT
				precision pitn = q[s].pitn;
			#else
				precision pitn = 0;
			#endif
				precision pixx = q[s].pixx;
				precision pixy = q[s].pixy;
			#ifndef BOOST_INVARIANT
				precision pixn = q[s].pixn;
			#else
				precision pixn = 0;
			#endif
				precision piyy = q[s].piyy;
			#ifndef BOOST_INVARIANT
				precision piyn = q[s].piyn;
			#else
				precision piyn = 0;
			#endif
				precision pinn = q[s].pinn;

				precision pi_mag = sqrt(fabs(pitt * pitt  +  pixx * pixx  +  piyy * piyy  +  t4 * pinn * pinn  -  2. * (pitx * pitx  +  pity * pity  -  pixy * pixy  +  t2 * (pitn * pitn  -  pixn * pixn  -  piyn * piyn))));

				precision factor_pi;

				// need to reorganize this
				switch(regulation_scheme)
				{
					case 0:
					{
						precision rho_pi = pi_mag / (rho_max * T_aniso_mag);

						if(!reprojection)	// include violations of orthogonality / traceless to rho_pi
						{
							precision trpi = fabs(pitt  -  pixx  -  piyy  -  t2 * pinn);

							precision piu0 = fabs(pitt * ut  -  pitx * ux  -  pity * uy  -  t2 * pitn * un) / ut;
							precision piu1 = fabs(pitx * ut  -  pixx * ux  -  pixy * uy  -  t2 * pixn * un) / ut;
							precision piu2 = fabs(pity * ut  -  pixy * ux  -  piyy * uy  -  t2 * piyn * un) / ut;
							precision piu3 = fabs(pitn * ut  -  pixn * ux  -  piyn * uy  -  t2 * pinn * un) * t / ut;

							precision piz0 = fabs(zt * pitt  -  t2 * zn * pitn) / (t * zn);
							precision piz1 = fabs(zt * pitx  -  t2 * zn * pixn) / (t * zn);
							precision piz2 = fabs(zt * pity  -  t2 * zn * piyn) / (t * zn);
							precision piz3 = fabs(zt * pitn  -  t2 * zn * pinn) / zn;

							precision violation = fmax(trpi, fmax(piu0, fmax(piu1, fmax(piu2, fmax(piu3, fmax(piz0, fmax(piz1, fmax(piz2, piz3))))))));

							rho_pi = fmax(rho_pi, violation / (xi0 * rho_max * pi_mag));
						}
						else
						{
							printf("regulate_residual_currents error: no reprojection yet\n");
							exit(-1);
						}

						factor_pi = tanh_function(rho_pi);

						break;
					}
					case 1:
					{
						printf("regulate_residual_currents error: no regulate_scheme = 1 yet\n");
						exit(-1);
						break;
					}
					default:
					{
						printf("regulate_residual_currents error: set regulate_scheme = (0,1)\n");
						exit(-1);
					}
				}
				q[s].pitt = factor_pi * pitt;
				q[s].pitx = factor_pi * pitx;
				q[s].pity = factor_pi * pity;
			#ifndef BOOST_INVARIANT
				q[s].pitn = factor_pi * pitn;
			#endif
				q[s].pixx = factor_pi * pixx;
				q[s].pixy = factor_pi * pixy;
			#ifndef BOOST_INVARIANT
				q[s].pixn = factor_pi * pixn;
			#endif
				q[s].piyy = factor_pi * piyy;
			#ifndef BOOST_INVARIANT
				q[s].piyn = factor_pi * piyn;
			#endif
				q[s].pinn = factor_pi * pinn;
			#endif

			#ifdef WTZMU
				precision WtTz = q[s].WtTz;
				precision WxTz = q[s].WxTz;
				precision WyTz = q[s].WyTz;
				precision WnTz = q[s].WnTz;

				precision WTz_mag = sqrt(fabs(WtTz * WtTz  -  WxTz * WxTz  -  WyTz * WyTz  -  t2 * WnTz * WnTz));

				precision factor_W;

				switch(regulation_scheme)
				{
					case 0:
					{
						precision rho_W = WTz_mag / (rho_max * T_aniso_mag);

						if(!reprojection)	// include violations of orthogonality to rho_W
						{
							precision WTzu = fabs(WtTz * ut  -  WxTz * ux  -  WyTz * uy  -  t2 * WnTz * un);
							precision WTzz = fabs(WtTz * zt  -  t2 * WnTz * zn);

							precision violation = fmax(WTzu, WTzz);

							rho_W = fmax(rho_W, violation / (xi0 * rho_max * WTz_mag));
						}
						else
						{
							printf("regulate_residual_currents error: no reprojection yet\n");
							exit(-1);
						}

						factor_W = tanh_function(rho_W);
						break;
					}
					case 1:
					{
						printf("regulate_residual_currents error: no regulate_scheme = 1 yet\n");
						exit(-1);
						break;
					}
					default:
					{
						printf("regulate_residual_currents error: set regulate_scheme = (0,1)\n");
						exit(-1);
					}
				}
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


void regulate_viscous_currents(precision t, hydro_variables * const __restrict__ q, precision * const __restrict__ e, const fluid_velocity * const __restrict__ u, lattice_parameters lattice, hydro_parameters hydro)
{
#if (NUMBER_OF_VISCOUS_CURRENTS != 0)
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision rho_max = hydro.rho_max;
	precision xi0 = hydro.xi0;

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
			#ifndef BOOST_INVARIANT
				precision un = u[s].un;
			#else
				precision un = 0;
			#endif
				precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t2 * un * un);

				precision Teq_mag = sqrt(e_s * e_s  +  3. * p * p);
		/*
			#ifdef PIMUNU
				precision pitt = q[s].pitt;
				precision pitx = q[s].pitx;
				precision pity = q[s].pity;
			#ifndef BOOST_INVARIANT
				precision pitn = q[s].pitn;
			#else
				precision pitn = 0;
			#endif
				precision pixx = q[s].pixx;
				precision pixy = q[s].pixy;
			#ifndef BOOST_INVARIANT
				precision pixn = q[s].pixn;
			#else
				precision pixn = 0;
			#endif
				precision piyy = q[s].piyy;
			#ifndef BOOST_INVARIANT
				precision piyn = q[s].piyn;
			#else
				precision piyn = 0;
			#endif
				precision pinn = q[s].pinn;

				precision pi_mag = sqrt(fabs(pitt * pitt  +  pixx * pixx  +  piyy * piyy  +  t4 * pinn * pinn  -  2. * (pitx * pitx  +  pity * pity  -  pixy * pixy  +  t2 * (pitn * pitn  -  pixn * pixn  -  piyn * piyn))));

				precision trpi = fabs(pitt  -  pixx  -  piyy  -  t2 * pinn);

				precision piu0 = fabs(pitt * ut  -  pitx * ux  -  pity * uy  -  t2 * pitn * un) / ut;
				precision piu1 = fabs(pitx * ut  -  pixx * ux  -  pixy * uy  -  t2 * pixn * un) / ut;
				precision piu2 = fabs(pity * ut  -  pixy * ux  -  piyy * uy  -  t2 * piyn * un) / ut;
				precision piu3 = fabs(pitn * ut  -  pixn * ux  -  piyn * uy  -  t2 * pinn * un) * t / ut;

				precision denom_pi = xi0 * rho_max * pi_mag;

				precision a0 = pi_mag / (rho_max * Teq_mag);
				precision a1 = trpi / denom_pi;
				precision a2 = piu0 / denom_pi;
				precision a3 = piu1 / denom_pi;
				precision a4 = piu2 / denom_pi;
				precision a5 = piu3 / denom_pi;

				precision rho_pi = fmax(a0, fmax(a1, fmax(a2, fmax(a3, fmax(a4, a5)))));

				precision factor_pi = tanh_function(rho_pi);

				q[s].pitt = factor_pi * pitt;
				q[s].pitx = factor_pi * pitx;
				q[s].pity = factor_pi * pity;
			#ifndef BOOST_INVARIANT
				q[s].pitn = factor_pi * pitn;
			#endif
				q[s].pixx = factor_pi * pixx;
				q[s].pixy = factor_pi * pixy;
			#ifndef BOOST_INVARIANT
				q[s].pixn = factor_pi * pixn;
			#endif
				q[s].piyy = factor_pi * piyy;
			#ifndef BOOST_INVARIANT
				q[s].piyn = factor_pi * piyn;
			#endif
				q[s].pinn = factor_pi * pinn;
			#endif



			#ifdef PI
				precision Pi = q[s].Pi;

				precision rho_bulk = fabs(Pi / (rho_max * Teq_mag));

				precision factor_bulk = tanh_function(rho_bulk);

				q[s].Pi = factor_bulk * Pi;
			#endif
	*/
			}
		}
	}
#endif
}




