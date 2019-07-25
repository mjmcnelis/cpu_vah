#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include "../include/InferredVariables.h"
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"
#include "../include/Parameters.h"
#include "../include/EquationOfState.h"

using namespace std;

#define TEST_TMUNU 0		// 1 for test

double ttt_error = 1.e-13;
double ttx_error = 1.e-13;
double tty_error = 1.e-13;
double ttn_error = 1.e-13;

inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}

void set_inferred_variables(const CONSERVED_VARIABLES * const __restrict__ q, precision * const __restrict__ e, FLUID_VELOCITY * const __restrict__ u, precision t, int nx, int ny, int nz)
{
	precision t2 = t * t;

	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				precision ttt = q->ttt[s];
				precision ttx = q->ttx[s];
				precision tty = q->tty[s];
				precision ttn = q->ttn[s];
				precision pl  = q->pl[s];
			#if (PT_MATCHING == 1)
				precision pt  = q->pt[s];
			#endif
			#ifdef PIMUNU
				precision pitt = q->pitt[s];
				precision pitx = q->pitx[s];
				precision pity = q->pity[s];
				precision pitn = q->pitn[s];
			#else
				precision pitt = 0, pitx = 0, pity = 0, pitn = 0;
			#endif
			#ifdef WTZMU
				precision WtTz = q->WtTz[s];
				precision WxTz = q->WxTz[s];
				precision WyTz = q->WyTz[s];
				precision WnTz = q->WnTz[s];
			#else
				precision WtTz = 0, WxTz = 0, WyTz = 0, WnTz = 0;
			#endif
				precision kt = ttt  -  pitt;			// [fm^-4]
				precision kx = ttx  -  pitx;			// [fm^-4]
				precision ky = tty  -  pity;			// [fm^-4]
				precision kn = ttn  -  pitn;			// [fm^-5]

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

				precision Mt = kt  -  2. * WtTz * zt;
				precision Mx = kx  -  WxTz * zt;
				precision My = ky  -  WyTz * zt;
				precision Mn = kn  -  WtTz * zn  -  WnTz * zt;

				// solution for e
			#if (PT_MATCHING == 1)
				precision Ltt = (pl - pt) * zt * zt;
				precision ut_numerator = Mt  +  pt  -  Ltt;

				precision e_s = fmax(E_MIN, Mt  -  Ltt  -  (Mx * Mx  +  My * My) / ut_numerator  -  t2 * Mn * Mn * ut_numerator / (Mt + pl) / (Mt + pl));
				// old formula (residual shear, zt = 0)
				// precision e_s = fmax(E_MIN, Mt  -  (Mx * Mx  +  My * My) / (Mt + pt));

			#else
				precision zt2  = zt * zt;
				precision ztzn = zt * zn;

				// what a pain in the ass (maybe double check)
				precision a = 0.25 * t2 * ztzn * ztzn  +  0.5 * (1. + zt2) * (1.  -  0.5 * zt2);
				precision b = 0.5 * t2 * Mn * ztzn  -  1.5 * t2 * pl * ztzn * ztzn  -  (0.5 * (1. + zt2) * (Mt - 1.5 * pl * zt2)  -  (1. - 0.5 * zt2) * (Mt - 0.5 * pl * (1. + 3. * zt2)));
				precision c = Mx * Mx  +  My * My  +  t2 * Mn * Mn  -  1.5 * t2 * Mn * pl * ztzn  +  2.25 * t2 * pl * pl * ztzn * ztzn  -  (Mt - 0.5 * pl * (1. + 3. * zt2)) * (Mt - 1.5 * pl * zt2);

				// quadratic formula
				precision e_s = fmax(E_MIN, (sqrt(fabs(b * b  -  4. * a * c))  -  b) / (2. * a));

				precision pt = (e_s - pl) / 2.;
				precision ut_numerator = Mt  +  pt  -  (pl - pt) * zt2;
				// this was the old quadratic formula (residual shear, zt = 0)
				//precision e_s = fmax(E_MIN, 0.5 * (-(Mt - pl)  +  sqrt(fabs((Mt - pl) * (Mt - pl)  +  8.0 * (Mt * (Mt - 0.5 * pl) - (Mx * Mx  +  My * My))))));
				//precision ut_numerator = Mt + pt;

			#endif

				// solution for u^mu
				precision ut_s = sqrt(fabs(ut_numerator / (e_s + pt)));
				precision ux_s = Mx / ut_s / (e_s + pt);
				precision uy_s = My / ut_s / (e_s + pt);
				precision un_s = F * ut_s;

				// if(std::isnan(ut_s))
				// {
				// 	printf("\nget_inferred_variables_test error: u^mu = (%lf, %lf, %lf, %lf) is nan\n", ut_s, ux_s, uy_s, un_s);
				// 	exit(-1);
				// }

				ut_s = sqrt(1.  +  ux_s * ux_s  +  uy_s * uy_s  +  t2 * un_s * un_s);

				// get solution for primary variables
				e[s]     = e_s;
				u->ut[s] = ut_s;
				u->ux[s] = ux_s;
				u->uy[s] = uy_s;
				u->un[s] = un_s;

				// test accuracy of reconstruction formula
			#if (TEST_TMUNU == 1)
				precision utperp_s = sqrt(1.  +  ux_s * ux_s  +  uy_s * uy_s);

				precision zt_s = t * un_s / utperp_s;
				precision zn_s = ut_s / (t * utperp_s);

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
					printf("get_inferred_variables: |dt^{tau/mu}| = (%.6g, %.6g, %.6g, %.6g)\n", ttt_error, ttx_error, tty_error, ttn_error);
				}
			#endif
			}
		}
	}
}

