#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include "../include/Precision.h"
#include "../include/Macros.h"
#include "../include/Hydrodynamics.h"
#include "../include/DynamicalVariables.h"
#include "../include/Parameters.h"
#include "../include/EquationOfState.h"
using namespace std;

double ttt_error = 1.e-13;
double ttx_error = 1.e-13;
double tty_error = 1.e-13;
double ttn_error = 1.e-13;

inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


void set_inferred_variables_aniso_hydro(const hydro_variables * const __restrict__ q, precision * const __restrict__ e, fluid_velocity * const __restrict__ u, precision t, lattice_parameters lattice, hydro_parameters hydro)
{
#ifdef ANISO_HYDRO
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision T_switch = hydro.freezeout_temperature_GeV;
	//precision e_switch = equilibrium_energy_density(T_switch / hbarc, hydro.conformal_eos_prefactor);
	precision e_switch = equilibrium_energy_density_new(T_switch / hbarc, hydro.conformal_eos_prefactor);

	//printf("%lf\n", e_switch * hbarc);

	precision e_min = hydro.energy_min;

	precision t2 = t * t;
	precision t4 = t2 * t2;

	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				precision ttt = q[s].ttt;
				precision ttx = q[s].ttx;
				precision tty = q[s].tty;
			#ifndef BOOST_INVARIANT
				precision ttn = q[s].ttn;
			#else
				precision ttn = 0;
			#endif

				precision pl  = q[s].pl;
				precision pt  = q[s].pt;

			#ifdef PIMUNU
				precision pitt = q[s].pitt;
				precision pitx = q[s].pitx;
				precision pity = q[s].pity;
			#ifndef BOOST_INVARIANT
				precision pitn = q[s].pitn;
			#else
				precision pitn = 0;
			#endif
			#else
				precision pitt = 0;
				precision pitx = 0;
				precision pity = 0;
				precision pitn = 0;
			#endif

			#ifdef WTZMU
				precision WtTz = q[s].WtTz;
				precision WxTz = q[s].WxTz;
				precision WyTz = q[s].WyTz;
				precision WnTz = q[s].WnTz;
			#else
				precision WtTz = 0;
				precision WxTz = 0;
				precision WyTz = 0;
				precision WnTz = 0;
			#endif

				precision kt = ttt  -  pitt;			// [fm^-4]
				precision kx = ttx  -  pitx;			// [fm^-4]
				precision ky = tty  -  pity;			// [fm^-4]
				precision kn = ttn  -  pitn;			// [fm^-5]

			#ifndef BOOST_INVARIANT
				precision A = kn / (kt + pl);		  	// [fm^-1] (check these formulas)
				precision B = WtTz / (t * (kt + pl));	// [fm^-1]

				precision t2A2 = t2 * A * A;
				precision t2B2 = t2 * B * B;

				// figure out how to deal with this
				precision F = (A  -  B * sqrt(fabs(1. + t2B2 - t2A2))) / (1. + t2B2);	// [fm^-1]
				precision Ft = t * F;													// [1]
				precision x = sqrt(fabs(1.  -  Ft * Ft));								// [1]

				precision zt = Ft / x;					// [1]
				precision zn = 1. / (x * t);			// [fm^-1]
			#else
				precision zt = 0;
				precision zn = 1. / t;
			#endif

				precision Mt = kt  -  2. * WtTz * zt;
				precision Mx = kx  -  WxTz * zt;
				precision My = ky  -  WyTz * zt;
				precision Mn = kn  -  WtTz * zn  -  WnTz * zt;

				// this needs a switching

				precision eprev = e[s];

				// solution for e
			#ifdef LATTICE_QCD
			#ifndef CONFORMAL_EOS
				precision Ltt = (pl - pt) * zt * zt;
				precision ut_numerator = Mt  +  pt  -  Ltt;

				precision e_s = energy_density_cutoff(e_min, Mt  -  Ltt  -  (Mx * Mx  +  My * My) / ut_numerator  -  t2 * Mn * Mn * ut_numerator / (Mt + pl) / (Mt + pl));
			#endif
			#endif

			#ifdef CONFORMAL_EOS
			#ifndef LATTICE_QCD
				precision zt2  = zt * zt;
				precision ztzn = zt * zn;

				// what a pain in the ass (so far testing shows it works though...)
				precision a = 0.25 * t2 * ztzn * ztzn  +  0.5 * (1. + zt2) * (1.  -  0.5 * zt2);
				precision b = 0.5 * t2 * Mn * ztzn  -  1.5 * t2 * pl * ztzn * ztzn  -  (0.5 * (1. + zt2) * (Mt - 1.5 * pl * zt2)  -  (1. - 0.5 * zt2) * (Mt - 0.5 * pl * (1. + 3. * zt2)));
				precision c = Mx * Mx  +  My * My  +  t2 * Mn * Mn  -  1.5 * t2 * Mn * pl * ztzn  +  2.25 * t2 * pl * pl * ztzn * ztzn  -  (Mt - 0.5 * pl * (1. + 3. * zt2)) * (Mt - 1.5 * pl * zt2);

				precision e_s = energy_density_cutoff(e_min, (sqrt(fabs(b * b  -  4. * a * c))  -  b) / (2. * a));

				pt = (e_s - pl) / 2.;

				precision ut_numerator = Mt  +  pt  -  (pl - pt) * zt2;
			#endif
			#endif

			#ifdef CONFORMAL_EOS
			#ifdef LATTICE_QCD
				printf("\nget_inferred_variables_aniso_hydro error: not eos switch here yet\n");
				exit(-1);
			#endif
			#endif

				// solution for u^mu
				precision ut_s = sqrt(fabs(ut_numerator / (e_s + pt)));
				precision ux_s = Mx / ut_s / (e_s + pt);
				precision uy_s = My / ut_s / (e_s + pt);
			#ifndef BOOST_INVARIANT
				precision un_s = F * ut_s;
			#else
				precision un_s = 0;
			#endif

				if(std::isnan(ut_s))
				{
					printf("\nget_inferred_variables_aniso_hydro error: u^mu = (%lf, %lf, %lf, %lf) is nan\n", ut_s, ux_s, uy_s, un_s);
					exit(-1);
				}
/*
				precision de_over_e = (e_s - eprev) / eprev;

				if(de_over_e > 0.1)
				{
					if(eprev > e_switch)
					{
						printf("(eprev, e_s, |de/e|) = (%.6g, %.6g, %.6g) at (%d, %d, %lf) \n", eprev, e_s, de_over_e, i, j, t);
					}

					e_s = 1.0 * eprev;

					e_regulation[s] += 1;
				}

				precision ut_max = 10.0;

				ut_s = sqrt(1.  +  ux_s * ux_s  +  uy_s * uy_s  +  t2 * un_s * un_s);

				if(ut_s > ut_max)
				{
					//printf("(ut, ut_max) = (%.5g, %.5g) at %d, %d \n", ut, ut_max, i, j);

					precision norm = fmin(1., sqrt(fabs((ut_max * ut_max - 1.) / (ut_s * ut_s - 1.))));

					ux_s *= norm;
					uy_s *= norm;
					un_s *= norm;

					ut_regulation[s] += 1;
				}
*/

			#ifdef TEST_TTAUMU
				ut_s = sqrt(1.  +  ux_s * ux_s  +  uy_s * uy_s  +  t2 * un_s * un_s);

			#ifndef BOOST_INVARIANT
				precision utperp_s = sqrt(1.  +  ux_s * ux_s  +  uy_s * uy_s);
				precision zt_s = t * un_s / utperp_s;
				precision zn_s = ut_s / (t * utperp_s);
			#else
				precision zt_s = 0;
				precision zn_s = 1. / t;
			#endif

				precision dttt = fabs((e_s + pt) * ut_s * ut_s  -  pt  +  (pl - pt) * zt_s * zt_s  +  2. * WtTz * zt_s  +  pitt  -  ttt);
				precision dttx = fabs((e_s + pt) * ut_s * ux_s  +  WxTz * zt_s  +  pitx  -  ttx);
				precision dtty = fabs((e_s + pt) * ut_s * uy_s  +  WyTz * zt_s  +  pity  -  tty);
				precision dttn = fabs((e_s + pt) * ut_s * un_s  +  (pl - pt) * zt_s * zn_s  +  WtTz * zn_s  +  WnTz * zt_s  +  pitn  -  ttn);

				if(dttt > ttt_error || dttx > ttx_error || dtty > tty_error || dttn > ttn_error)
				{
					ttt_error = fmax(dttt, ttt_error);
					ttx_error = fmax(dttx, ttx_error);
					tty_error = fmax(dtty, tty_error);
					ttn_error = fmax(dttn, ttn_error);
					printf("get_inferred_variables_aniso_hydro: |dt^{tau/mu}| = (%.6g, %.6g, %.6g, %.6g)\n", ttt_error, ttx_error, tty_error, ttn_error);
				}
			#endif

				e[s]    = e_s;		// cutoff the solution further so that e > E_MIN
				u[s].ux = ux_s;
				u[s].uy = uy_s;
			#ifndef BOOST_INVARIANT
				u[s].un = un_s;
			#endif
			}
		}
	}
#endif
}



void set_inferred_variables_viscous_hydro(const hydro_variables * const __restrict__ q, precision * const __restrict__ e, fluid_velocity * const __restrict__ u, precision t, lattice_parameters lattice, hydro_parameters hydro)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision T_switch = hydro.freezeout_temperature_GeV;
	//precision e_switch = equilibrium_energy_density(T_switch / hbarc, hydro.conformal_eos_prefactor);
	precision e_switch = equilibrium_energy_density_new(T_switch / hbarc, hydro.conformal_eos_prefactor);

	//printf("%lf\n", e_switch * hbarc);

	precision e_min = hydro.energy_min;

	precision t2 = t * t;

	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				precision ttt = q[s].ttt;
				precision ttx = q[s].ttx;
				precision tty = q[s].tty;
			#ifndef BOOST_INVARIANT
				precision ttn = q[s].ttn;
			#else
				precision ttn = 0;
			#endif

			#ifdef PIMUNU
				precision pitt = q[s].pitt;
				precision pitx = q[s].pitx;
				precision pity = q[s].pity;
			#ifndef BOOST_INVARIANT
				precision pitn = q[s].pitn;
			#else
				precision pitn = 0;
			#endif
			#else
				precision pitt = 0;
				precision pitx = 0;
				precision pity = 0;
				precision pitn = 0;
			#endif

			#ifdef PI
				precision Pi = q[s].Pi;
			#else
				precision Pi = 0;
			#endif

				precision Mt = ttt  -  pitt;
				precision Mx = ttx  -  pitx;
				precision My = tty  -  pity;
				precision Mn = ttn  -  pitn;

				precision M_squared = Mx * Mx  +  My * My  +  t2 * Mn * Mn;

				precision eprev = e[s];

			#ifdef CONFORMAL_EOS
				//precision e_s = energy_density_cutoff(e_min, - Mt  +  sqrt(fabs(4. * Mt * Mt  -  3. * M_squared)));
				precision e_s = - Mt  +  sqrt(fabs(4. * Mt * Mt  -  3. * M_squared));
			#else

				precision e_s = eprev;					// initial guess is previous energy density
				precision de;

				int n;

				const int max_iterations = 20;
				const double energy_tolerance = 1.e-5;

				for(n = 1; n <= max_iterations; n++)	// root solving algorithm (update e)
				{
					equation_of_state EoS(e_s);
					precision p = EoS.equilibrium_pressure();

					if(p + Pi <= 0.) 					// solution when have to regulate bulk pressure
					{
						e_s = Mt  -  M_squared / Mt;
						break;
					}

					precision cs2 = EoS.speed_of_sound_squared();

					precision f = (Mt - e_s) * (Mt + p + Pi)  -  M_squared;

					precision fprime = cs2 * (Mt - e_s)  -  (Mt + p + Pi);

					de = - f / fprime;

					e_s += de;

					if(e_s < e_min)						// stop iterating if e < e_min
					{
						break;
					}
					else if(fabs(de / e_s) <= energy_tolerance)
					{
						break;
					}
				}

				// if(n > max_iterations)
				// {
				// 	printf("newton method (eprev, e_s, |de/e_s|) = (%.6g, %.6g, %.6g) failed to converge within desired percentage tolerance %lf at (i, j, k) = (%d, %d, %d)\n", eprev, e_s, fabs(de / e_s), energy_tolerance, i, j, k);
				// }


			#endif

				e_s = energy_density_cutoff(e_min, e_s);

				equation_of_state eos(e_s);
				precision p = eos.equilibrium_pressure();

				precision P = fmax(0., p + Pi);								// should I smooth regulate it?

				precision ut = sqrt(fabs((Mt + P) / (e_s + P)));
				precision ux = Mx / ut / (e_s + P);
				precision uy = My / ut / (e_s + P);
			#ifndef BOOST_INVARIANT
				precision un = Mn / ut / (e_s + P);
			#else
				precision un = 0;
			#endif

				if(std::isnan(e_s) || std::isnan(ut))
				{
					printf("\nget_inferred_variables_viscous_hydro error: (e, ut, P, Mt) = (%lf, %lf, %lf, %lf) \n", e_s, ut, P, Mt);
					exit(-1);
				}

/*
				precision de_over_e = (e_s - eprev) / eprev;

				if(de_over_e > 0.1)
				{
					if(eprev > e_switch)
					{
						printf("(eprev, e_s, |de/e|) = (%.6g, %.6g, %.6g) at (%d, %d, %lf) \n", eprev, e_s, de_over_e, i, j, t);
					}

					e_s = 1.0 * eprev;

					e_regulation[s] += 1;
				}


				precision ut_max = 10.0;

				ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t2 * un * un);

				if(ut > ut_max)
				{
					//printf("(ut, ut_max) = (%.5g, %.5g) at %d, %d \n", ut, ut_max, i, j);

					precision norm = fmin(1., sqrt(fabs((ut_max * ut_max - 1.) / (ut * ut - 1.))));

					ux *= norm;
					uy *= norm;
					un *= norm;

					ut_regulation[s] += 1;
				}
*/


			#ifdef TEST_TTAUMU
				ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t2 * un * un);

				precision dttt = fabs((e_s + p + Pi) * ut * ut  -  p  -  Pi +  pitt  -  ttt);
				precision dttx = fabs((e_s + p + Pi) * ut * ux  +  pitx  -  ttx);
				precision dtty = fabs((e_s + p + Pi) * ut * uy  +  pity  -  tty);
				precision dttn = fabs((e_s + p + Pi) * ut * un  +  pitn  -  ttn);

				if(dttt > ttt_error || dttx > ttx_error || dtty > tty_error || dttn > ttn_error)
				{
					ttt_error = fmax(dttt, ttt_error);
					ttx_error = fmax(dttx, ttx_error);
					tty_error = fmax(dtty, tty_error);
					ttn_error = fmax(dttn, ttn_error);
					printf("get_inferred_variables_viscous_hydro: |dt^{tau/mu}| = (%.6g, %.6g, %.6g, %.6g) at %d, %d, %d \n", ttt_error, ttx_error, tty_error, ttn_error, i, j, k);
				}
			#endif

				e[s]    = e_s;		// set solution for primary variables
				u[s].ux = ux;
				u[s].uy = uy;
			#ifndef BOOST_INVARIANT
				u[s].un = un;
			#endif
			}
		}
	}
}

