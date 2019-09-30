#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include "../include/Precision.h"
#include "../include/Macros.h"
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

inline int sign(precision x)
{
	if(x > 0.) return 1;
	else if(x < 0.) return -1;
	else return 0;
}

// maybe need more quantitative way of checking convergence of iterations
// maybe regulate individual rows of pimunu
// here maybe: I don't need to reproject pimunu b/c already "fixed" velocity solution (just use "old method")

void set_inferred_variables_aniso_hydro(const hydro_variables * const __restrict__ q, precision * const __restrict__ e, fluid_velocity * const __restrict__ u, precision t, lattice_parameters lattice, hydro_parameters hydro)
{
#ifdef ANISO_HYDRO
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

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

			#if (PT_MATCHING == 1)
				precision pt  = q[s].pt;
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

				// solution for e
			#if (PT_MATCHING == 1)
				precision Ltt = (pl - pt) * zt * zt;
				precision ut_numerator = Mt  +  pt  -  Ltt;

				precision e_s = energy_density_cutoff(e_min, Mt  -  Ltt  -  (Mx * Mx  +  My * My) / ut_numerator  -  t2 * Mn * Mn * ut_numerator / (Mt + pl) / (Mt + pl));

				precision ut_numerator = Mt  +  pt  -  (pl - pt) * zt2;
			#else
				precision zt2  = zt * zt;
				precision ztzn = zt * zn;

				// what a pain in the ass (so far testing shows it works though...)
				precision a = 0.25 * t2 * ztzn * ztzn  +  0.5 * (1. + zt2) * (1.  -  0.5 * zt2);
				precision b = 0.5 * t2 * Mn * ztzn  -  1.5 * t2 * pl * ztzn * ztzn  -  (0.5 * (1. + zt2) * (Mt - 1.5 * pl * zt2)  -  (1. - 0.5 * zt2) * (Mt - 0.5 * pl * (1. + 3. * zt2)));
				precision c = Mx * Mx  +  My * My  +  t2 * Mn * Mn  -  1.5 * t2 * Mn * pl * ztzn  +  2.25 * t2 * pl * pl * ztzn * ztzn  -  (Mt - 0.5 * pl * (1. + 3. * zt2)) * (Mt - 1.5 * pl * zt2);

				precision e_s = energy_density_cutoff(e_min, (sqrt(fabs(b * b  -  4. * a * c))  -  b) / (2. * a));

				precision pt = (e_s - pl) / 2.;

				precision ut_numerator = Mt  +  pt  -  (pl - pt) * zt2;
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

	precision e_min = hydro.energy_min;

	equation_of_state eos_min(e_min);
	precision p_min = eos_min.equilibrium_pressure();
	precision cs2_min = eos_min.speed_of_sound_squared();

	//printf("p_min = %.3g\n", p_min);
	//printf("cs2_min = %.3g\n", cs2_min);

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

			#ifdef CONFORMAL_EOS
				precision e_s = energy_density_cutoff(e_min, - Mt  +  sqrt(fabs(4. * Mt * Mt  -  3. * M_squared)));
			#else
				precision e_s;



/*
				bool non_emin_solution = true;

			#ifdef PI
				if((Mt - e_min) * (Mt + p_min + Pi)  -  M_squared < 0.)
				{
					Pi = (M_squared  -  (Mt + p_min) * (Mt - e_min)) / (Mt - e_min);

					if(((cs2_min - 1.) * Mt  -  Pi) <= 0.)
					{
						e_s = e_min;
						non_emin_solution = false;
					}
				}
			#endif

				if(non_emin_solution)
				{
					e_s = e[s];				// initial guess is previous energy density
					precision eprev = e_s;	// for tracking convergence

					// debug the bulk pressure here
					int n;
					const int n_max = 100;
					for(n = 1; n <= n_max; n++)			// root solving algorithm (update e)
					{
						equation_of_state EoS(e_s);
						precision p = EoS.equilibrium_pressure();
						precision cs2 = EoS.speed_of_sound_squared();
						precision cst2 = p / e_s;

						// should I regulate Pi here?

						// why even use this function manipulation at all?
						precision A = Mt * (1. - cst2)  +  Pi;
						precision B = Mt * (Mt + Pi)  -  M_squared;
						precision H = sqrt(fabs(A * A  +  4. * cst2 * B));
						precision D = (A - H) / (2. * cst2);

						precision f = e_s + D;
						precision fprime = 1. - ((cs2 - cst2) * (B  +  D * H  -  ((cs2 - cst2) * cst2 * D * Mt) / e_s)) / (cst2 * e_s * H);

						precision de = - f / fprime;

						eprev = e_s;

						e_s += de;

						if(e_s < e_min) break;
						if(fabs(de) <= 1.e-3 * fmax(1., e_s)) break;
					}

					if(n > n_max) printf("newton method (eprev, e_s) = (%lf, %lf) failed to converge at (i, j, k) = (%d, %d, %d)\n", eprev, e_s, i, j, k);
					//printf("%d\n", n);

				}
*/


			/*
				if(Mt <= e_min)  			// since e < Mt
				{
					e_s = e_min;			// is this necessary?
				}
				else 						// check that the is bracketed
				{

					precision eL = e_min;
					precision eR = 1.001 * Mt;

					equation_of_state eosL(eL);
					equation_of_state eosR(eR);

					precision pL = eosL.equilibrium_pressure();
					precision pR = eosR.equilibrium_pressure();

					// let's try to regulate bulk pressure

					precision fL = Mt  -  eL  -  M_squared / (Mt + pL + Pi);
					precision fR = Mt  -  eR  -  M_squared / (Mt + pR + Pi);

					int sign1 = sign(fL);
					int sign2 = sign(fR);

					if(sign2 != -1) 	// the third term in fR might be negative if Pi is negative enough
					{
						printf("Root bracketing error: fR = %.4g is not negative (pR, Mt, Pi) = (%.3g, %.3g, %.3g)\n", fR, pR, Mt, Pi);
						exit(-1);
					}

					if(sign1 == sign2)		// the root is not bracketed in the region [e_min, 1.01 * Mt]
					//if(false)
					{
						e_s = e_min;
					}
					else
					{
						e_s = e[s];				// initial guess is previous energy density
						precision eprev = e_s;	// for tracking convergence

						// debug the bulk pressure here
						int n;
						const int n_max = 100;
						for(n = 1; n <= n_max; n++)			// root solving algorithm (update e)
						{
							equation_of_state EoS(e_s);
							precision p = EoS.equilibrium_pressure();
							precision cs2 = EoS.speed_of_sound_squared();
							precision cst2 = p / e_s;

							// should I regulate Pi here?

							// why even use this function manipulation at all?
							precision A = Mt * (1. - cst2)  +  Pi;
							precision B = Mt * (Mt + Pi)  -  M_squared;
							precision H = sqrt(fabs(A * A  +  4. * cst2 * B));
							precision D = (A - H) / (2. * cst2);

							precision f = e_s + D;
							precision fprime = 1. - ((cs2 - cst2) * (B  +  D * H  -  ((cs2 - cst2) * cst2 * D * Mt) / e_s)) / (cst2 * e_s * H);

							precision de = - f / fprime;

							eprev = e_s;

							e_s += de;

							if(e_s < e_min) break;
							if(fabs(de) <= 1.e-3 * fmax(1., e_s)) break;
						}

						if(n > n_max) printf("newton method (eprev, e_s) = (%lf, %lf) failed to converge at (i, j, k) = (%d, %d, %d)\n", eprev, e_s, i, j, k);
						//printf("%d\n", n);
					}
				}
			*/

				if(std::isnan(e_s))
				{
					printf("\nget_inferred_variables_viscous_hydro error: e = %lf\n", e_s);
					exit(-1);
				}

			#endif
				e_s = energy_density_cutoff(e_min, e_s);

				equation_of_state eos(e_s);
				precision p = eos.equilibrium_pressure();

				precision ut = sqrt(fabs((Mt + p + Pi) / (e_s + p + Pi)));
				precision ux = Mx / ut / (e_s + p + Pi);
				precision uy = My / ut / (e_s + p + Pi);
			#ifndef BOOST_INVARIANT
				precision un = Mn / ut / (e_s + p + Pi);
			#else
				precision un = 0;
			#endif

				if(std::isnan(ut))
				{
					printf("\nget_inferred_variables_viscous_hydro error: u^mu = (%lf, %lf, %lf, %lf)\n", ut, ux, uy, un);
					exit(-1);
				}

			#ifdef TEST_TTAUMU
				ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t2 * un * un);

				precision dttt = fabs((e_s + p + Pi) * ut * ut  -  (p + Pi)  +  pitt  -  ttt);
				precision dttx = fabs((e_s + p + Pi) * ut * ux  +  pitx  -  ttx);
				precision dtty = fabs((e_s + p + Pi) * ut * uy  +  pity  -  tty);
				precision dttn = fabs((e_s + p + Pi) * ut * un  +  pitn  -  ttn);

				if(dttt > ttt_error || dttx > ttx_error || dtty > tty_error || dttn > ttn_error)
				{
					ttt_error = fmax(dttt, ttt_error);
					ttx_error = fmax(dttx, ttx_error);
					tty_error = fmax(dtty, tty_error);
					ttn_error = fmax(dttn, ttn_error);
					printf("get_inferred_variables_viscous_hydro: |dt^{tau/mu}| = (%.6g, %.6g, %.6g, %.6g)\n", ttt_error, ttx_error, tty_error, ttn_error);
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

