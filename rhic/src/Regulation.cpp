#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include "../include/Macros.h"
#include "../include/DynamicalVariables.h"
#include "../include/InferredVariables.h"
#include "../include/Projections.h"

const precision sqrt_two   = sqrt(2.);
const precision sqrt_three = sqrt(3.);

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
	precision pressure_min = hydro.pressure_min;
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


				if(pl < pressure_min)			// regulate longitudinal and transverse pressures
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

			#if (NUMBER_OF_RESIDUAL_CURRENTS != 0)		// regulate residual currents

				precision ux = u[s].ux;
				precision uy = u[s].uy;
			#ifndef BOOST_INVARIANT
				precision un = u[s].un;
			#else
				precision un = 0;
			#endif
				precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t2 * un * un);
				precision utperp = sqrt(1.  +  ux * ux  +  uy * uy);

			#ifndef BOOST_INVARIANT
				precision zt = t * un / utperp;
				precision zn = ut / t / utperp;
			#else
				precision zt = 0;
				precision zn = 1. / t;
			#endif

				precision Taniso = sqrt_two * pt;
			#ifdef WTZMU
				Taniso = sqrt(pl * pl  +  2. * pt * pt);
			#endif

				// 2+1d: pitt, pitx, pity, pixx, pixy (piyy = f(pixx, piyy, ux, uy); pixn, piyn, pitn, pinn = 0)
				// 3+1d: pitt, pitx, pity, pitn, pixx, pixy (pixn, piyy, piyn, pinn can be reconstructed algebraically, need chain rule for derivatives)

				// independent components are pixx and pixy

			#ifdef PIMUNU
				precision pitt = q[s].pitt;
				precision pitx = q[s].pitx;
				precision pity = q[s].pity;
				precision pixx = q[s].pixx;
				precision pixy = q[s].pixy;
				precision piyy = q[s].piyy;

			#ifndef BOOST_INVARIANT
				precision pitn = q[s].pitn;
				precision pixn = q[s].pixn;
				precision piyn = q[s].piyn;
				precision pinn = q[s].pinn;
			#else
				precision pitn = 0;
				precision pixn = 0;
				precision piyn = 0;
				precision pinn = 0;
			#endif
				// enforce orthogonality and tracelessness
				piyy = (- pixx * (1.  +  uy * uy)  +  2. * pixy * ux * uy) / (1.  +  ux * ux);
				pitx = (pixx * ux  +  pixy * uy) * ut / (utperp * utperp);
				pity = (pixy * ux  +  piyy * uy) * ut / (utperp * utperp);
				pixn = pitx * un / ut;
				piyn = pity * un / ut;
				pitn = (pixn * ux  +  piyn * uy) * ut / (utperp * utperp);
				pinn = pitn * un / ut;
				pitt = (pitx * ux  +  pity * uy  +  t2 * pitn * un) / ut;
			#else
				precision pitt = 0;
				precision pitx = 0;
				precision pity = 0;
				precision pitn = 0;
				precision pixx = 0;
				precision pixy = 0;
				precision pixn = 0;
				precision piyy = 0;
				precision piyn = 0;
				precision pinn = 0;
			#endif

			#ifdef WTZMU
				precision WtTz = q[s].WtTz;
				precision WxTz = q[s].WxTz;
				precision WyTz = q[s].WyTz;
				precision WnTz = q[s].WnTz;

				// enforce orthogonality
				WtTz = (WxTz * ux  +  WyTz * uy) * ut / (utperp * utperp);
				WnTz = WtTz * un / ut;
			#else
				precision WtTz = 0;
				precision WxTz = 0;
				precision WyTz = 0;
				precision WnTz = 0;
			#endif

				precision factor_pi = 1;
				precision factor_W = 1;

				switch(regulation_scheme)
				{
					case 0:
					{
						precision pi_mag = sqrt(fabs(pitt * pitt  +  pixx * pixx  +  piyy * piyy  +  t4 * pinn * pinn  -  2. * (pitx * pitx  +  pity * pity  -  pixy * pixy  +  t2 * (pitn * pitn  -  pixn * pixn  -  piyn * piyn))));
						precision WTz_mag = sqrt(2. * fabs(WtTz * WtTz  -  WxTz * WxTz  -  WyTz * WyTz  -  t2 * WnTz * WnTz));

						factor_pi = tanh_function(pi_mag / (rho_max * Taniso));
						factor_W = tanh_function(WTz_mag / (rho_max * Taniso));

						break;
					}
					case 1:
					{
						precision Tres = sqrt(fabs(2. * (WtTz * WtTz  -  WxTz * WxTz  -  WyTz * WyTz  -  t2 * WnTz * WnTz)  +  pitt * pitt  +  pixx * pixx  +  piyy * piyy  +  t4 * pinn * pinn  -  2. * (pitx * pitx  +  pity * pity  -  pixy * pixy  +  t2 * (pitn * pitn  -  pixn * pixn  -  piyn * piyn))));

 						precision factor = Taniso / (1.e-8 + Tres);

 						if(factor < 1.)
 						{
 							factor_pi = factor;
 							factor_W = factor;
 						}
						break;
					}
					default:
					{
						printf("regulate_residual_currents error: set regulate_scheme = (0,1)\n");
						exit(-1);
					}
				}

			#ifdef PIMUNU
				q[s].pitt = factor_pi * pitt;
				q[s].pitx = factor_pi * pitx;
				q[s].pity = factor_pi * pity;
				q[s].pixx = factor_pi * pixx;
				q[s].pixy = factor_pi * pixy;
				q[s].piyy = factor_pi * piyy;

			#ifndef BOOST_INVARIANT
				q[s].pitn = factor_pi * pitn;
				q[s].pixn = factor_pi * pixn;
				q[s].piyn = factor_pi * piyn;
				q[s].pinn = factor_pi * pinn;
			#endif
			#endif

			#ifdef WTZMU
				q[s].WtTz = factor_W * WtTz;
				q[s].WxTz = factor_W * WxTz;
				q[s].WyTz = factor_W * WyTz;
				q[s].WnTz = factor_W * WnTz;
			#endif

			#endif
			}
		}
	}
#endif
}


void regulate_viscous_currents(precision t, hydro_variables * const __restrict__ q, precision * const __restrict__ e, const fluid_velocity * const __restrict__ u, lattice_parameters lattice, hydro_parameters hydro)
{
#ifndef ANISO_HYDRO
#if (NUMBER_OF_VISCOUS_CURRENTS != 0)
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	int regulation_scheme = hydro.regulation_scheme;
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

				equation_of_state eos(e[s]);
				precision p = eos.equilibrium_pressure();

				precision ux = u[s].ux;
				precision uy = u[s].uy;

			#ifndef BOOST_INVARIANT
				precision un = u[s].un;
			#else
				precision un = 0;
			#endif

				precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t2 * un * un);
				precision utperp = sqrt(1.  +  ux * ux  +  uy * uy);

				precision Teq = sqrt_three * p;

			#ifdef PIMUNU
				precision pitt = q[s].pitt;
				precision pitx = q[s].pitx;
				precision pity = q[s].pity;
				precision pixx = q[s].pixx;
				precision pixy = q[s].pixy;
				precision piyy = q[s].piyy;
				precision pinn = q[s].pinn;

			#ifndef BOOST_INVARIANT
				precision pitn = q[s].pitn;
				precision pixn = q[s].pixn;
				precision piyn = q[s].piyn;
			#else
				precision pitn = 0;
				precision pixn = 0;
				precision piyn = 0;
			#endif
				// enforce orthogonality and tracelessness
				pinn = (pixx * (ux * ux  -  ut * ut)  +  piyy * (uy * uy  -  ut * ut)  +  2. * (pixy * ux * uy  +  t2 * un * (pixn * ux  +  piyn * uy))) / (t2 * utperp * utperp);
         		pitn = (pixn * ux  +  piyn * uy  +  t2 * un * pinn) / ut;
				pity = (pixy * ux  +  piyy * uy  +  t2 * un * piyn) / ut;
				pitx = (pixx * ux  +  pixy * uy  +  t2 * un * pixn) / ut;
				pitt = (pitx * ux  +  pity * uy  +  t2 * un * pitn) / ut;
			#else
				precision pitt = 0;
				precision pitx = 0;
				precision pity = 0;
				precision pitn = 0;
				precision pixx = 0;
				precision pixy = 0;
				precision pixn = 0;
				precision piyy = 0;
				precision piyn = 0;
				precision pinn = 0;
			#endif

			#ifdef PI
				precision Pi = q[s].Pi;
			#else
				precision Pi = 0;
			#endif

 				precision factor_pi = 1;
 				precision factor_bulk = 1;

				// regulation of viscous pressures
				switch(regulation_scheme)
				{
					case 0:		// gradually regulate the viscous pressures at all grid points
					{
						precision pi_mag = sqrt(fabs(pitt * pitt  +  pixx * pixx  +  piyy * piyy  +  t4 * pinn * pinn  -  2. * (pitx * pitx  +  pity * pity  -  pixy * pixy  +  t2 * (pitn * pitn  -  pixn * pixn  -  piyn * piyn))));

						factor_pi = tanh_function(pi_mag / (rho_max * Teq));
						factor_bulk = tanh_function(fabs(sqrt_three * Pi / (rho_max * Teq)));

						break;
					}
					case 1:		// only do regulation if viscous Tmunu magnitude > equilibrium part
					{
						precision Tvisc = sqrt(3. * Pi * Pi  + fabs(pitt * pitt  +  pixx * pixx  +  piyy * piyy  +  t4 * pinn * pinn  -  2. * (pitx * pitx  +  pity * pity  -  pixy * pixy  +  t2 * (pitn * pitn  -  pixn * pixn  -  piyn * piyn))));

						precision factor = Teq / (1.e-10 + Tvisc);

						if(factor < 1.)
						{
							factor_pi = factor;
							factor_bulk = factor;
						}

						break;
					}
					default:
					{
						printf("regulate_viscous_currents error: set regulate_scheme = (0,1)\n");
						exit(-1);
					}
				}

				#ifdef PIMUNU
					q[s].pitt = factor_pi * pitt;
					q[s].pitx = factor_pi * pitx;
					q[s].pity = factor_pi * pity;
					q[s].pixx = factor_pi * pixx;
					q[s].pixy = factor_pi * pixy;
					q[s].piyy = factor_pi * piyy;
					q[s].pinn = factor_pi * pinn;

				#ifndef BOOST_INVARIANT
					q[s].pitn = factor_pi * pitn;
					q[s].pixn = factor_pi * pixn;
					q[s].piyn = factor_pi * piyn;
				#endif
				#endif

				#ifdef PI
					q[s].Pi = factor_bulk * Pi;
				#endif
			}
		}
	}
#endif
#endif
}




