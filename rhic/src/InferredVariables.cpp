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

double ttt_error = 1.e-13;
double ttx_error = 1.e-13;
double tty_error = 1.e-13;
double ttn_error = 1.e-13;

inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}

void test_energy_momentum_tensor_error(const precision * const __restrict__ q, precision e_s, precision ut, precision ux, precision uy, precision un, precision t)
{
	precision ttt = q[0];
	precision ttx = q[1];
	precision tty = q[2];
	precision ttn = q[3];
	precision pl  = q[4];

#if (PT_MATCHING == 0)
	precision pt = 0.5 * (e_s - pl);	
#else
	precision pt = q[5];
#endif

#ifdef PIMUNU
	precision pitt = q[6];
	precision pitx = q[7];
	precision pity = q[8];
	precision pitn = q[9];
#else
	precision pitt = 0.0;
	precision pitx = 0.0;
	precision pity = 0.0;
	precision pitn = 0.0;
#endif

#ifdef W_TZ_MU
	precision WtTz = q[10];  
	precision WxTz = q[11];  
	precision WyTz = q[12];  
	precision WnTz = q[13];
#else
	precision WtTz = 0.0;
	precision WxTz = 0.0;
	precision WyTz = 0.0;
	precision WnTz = 0.0;
#endif

	precision utperp = sqrt(1.0  +  ux * ux  +  uy * uy);

	precision zt = t * un / utperp;
	precision zn = ut / (t * utperp);

	precision dttt = fabs((e_s + pt) * ut * ut  -  pt  +  (pl - pt) * zt * zt  +  2.0 * WtTz * zt  +  pitt  -  ttt);
	precision dttx = fabs((e_s + pt) * ut * ux  +  WxTz * zt  +  pitx  -  ttx);
	precision dtty = fabs((e_s + pt) * ut * uy  +  WyTz * zt  +  pity  -  tty);
	precision dttn = fabs((e_s + pt) * ut * un  +  (pl - pt) * zt * zn  +  WtTz * zn  +  WnTz * zt  +  pitn  -  ttn);

	if(dttt > ttt_error || dttx > ttx_error || dtty > tty_error || dttn > ttn_error)
	{
		ttt_error = fmax(dttt, ttt_error);
		ttx_error = fmax(dttx, ttx_error);
		tty_error = fmax(dtty, tty_error);
		ttn_error = fmax(dttn, ttn_error);

		printf("get_inferred_variables: |dt^{tau/mu}| = (%.6g, %.6g, %.6g, %.6g)\n", ttt_error, ttx_error, tty_error, ttn_error);
	}

}




void test_energy_momentum_tensor_error_new(precision ttt, precision ttx, precision tty, precision ttn, precision pl, precision pt, precision pitt, precision pitx, precision pity, precision pitn, precision WtTz, precision WxTz, precision WyTz, precision WnTz, precision e_s, precision ut, precision ux, precision uy, precision un, precision t)
{
	precision utperp = sqrt(1.0  +  ux * ux  +  uy * uy);

	precision zt = t * un / utperp;
	precision zn = ut / (t * utperp);

	precision dttt = fabs((e_s + pt) * ut * ut  -  pt  +  (pl - pt) * zt * zt  +  2.0 * WtTz * zt  +  pitt  -  ttt);
	precision dttx = fabs((e_s + pt) * ut * ux  +  WxTz * zt  +  pitx  -  ttx);
	precision dtty = fabs((e_s + pt) * ut * uy  +  WyTz * zt  +  pity  -  tty);
	precision dttn = fabs((e_s + pt) * ut * un  +  (pl - pt) * zt * zn  +  WtTz * zn  +  WnTz * zt  +  pitn  -  ttn);

	if(dttt > ttt_error || dttx > ttx_error || dtty > tty_error || dttn > ttn_error)
	{
		ttt_error = fmax(dttt, ttt_error);
		ttx_error = fmax(dttx, ttx_error);
		tty_error = fmax(dtty, tty_error);
		ttn_error = fmax(dttn, ttn_error);

		printf("get_inferred_variables: |dt^{tau/mu}| = (%.6g, %.6g, %.6g, %.6g)\n", ttt_error, ttx_error, tty_error, ttn_error);
	}

}


// this is a temporary test for conformal PL matching
void get_inferred_variables_test(precision t, const precision * const __restrict__ q, precision * const __restrict__ e, precision * const __restrict__ ut, precision * const __restrict__ ux, precision * const __restrict__ uy, precision * const __restrict__ un)
{
	precision ttt = q[0];
	precision ttx = q[1];
	precision tty = q[2];
	precision ttn = q[3];
	precision pl  = q[4];

#ifdef PIMUNU
	precision pitt = q[6];
	precision pitx = q[7];
	precision pity = q[8];
	precision pitn = q[9];
#else
	precision pitt = 0.0;
	precision pitx = 0.0;
	precision pity = 0.0;
	precision pitn = 0.0;
#endif

#ifdef W_TZ_MU
	WtTz = q[10];  
	WxTz = q[11];
	WyTz = q[12];
	WnTz = q[13];
#else
	precision WtTz = 0.0;
	precision WxTz = 0.0;
	precision WyTz = 0.0;
	precision WnTz = 0.0;
#endif

	precision Mt = ttt;
	precision Mx = ttx;
	precision My = tty;
	precision Mn = ttn;

// reconstruction formula for energy density
#if (PT_MATCHING == 0)

	// this needs to be further generalized to include Mn
	precision e_s = 0.5 * (-(Mt - pl)  +  sqrt(fabs((Mt - pl) * (Mt - pl)  +  8.0 * (Mt * (Mt - 0.5 * pl) - (Mx * Mx  +  My * My)))));

	e_s = fmax(E_MIN, e_s);

	precision pt = 0.5 * (e_s - pl);							// conformal formula

#else
	precision pt  = q[5];

	precision e_s = Mt  -  (Mx * Mx  +  My * My) / (Mt + pt);	// generalize

	e_s = fmax(E_MIN, e_s);

#endif

	// fluid velocity (generalize)
	precision ut_s = sqrt(fabs((Mt + pt) / (e_s + pt)));
	precision ux_s = Mx / ut_s / (e_s + pt);
	precision uy_s = My / ut_s / (e_s + pt);
	precision un_s = 0.0;

	if(std::isnan(ut_s))
	{
		printf("\nget_inferred_variables_test error: u^mu = (%lf, %lf, %lf, %lf) is nan\n", ut_s, ux_s, uy_s, un_s);
		exit(-1);
	}

	// renormalize
	ut_s = sqrt(1.0  +  ux_s * ux_s  +  uy_s * uy_s  +  t * t * un_s * un_s);

	// get solution for primary variables
	*e  = e_s;
	*ut = ut_s;
	*ux = ux_s;
	*uy = uy_s;
	*un = un_s;

	//test_energy_momentum_tensor_error(q, e_s, ut_s, ux_s, uy_s, un_s, t);
}






/*
void set_inferred_variables(const CONSERVED_VARIABLES * const __restrict__ q, precision * const __restrict__ e, FLUID_VELOCITY * const __restrict__ u, precision t, int nx, int ny, int nz)
{
	int N = 5;  // default size for q_s array (ttt, ttx, tty, ttn, pl)

#if (PT_MATCHING == 1)
	N += 1;
#endif

#ifdef PIMUNU
	N += 4;		// add time components of pimunu and Wmu
#endif

#ifdef W_TZ_MU
	N += 4;
#endif

	const int nq = N;
	precision q_s[nq];

	// loop over the physical grid points
	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				precision e_s, ut_s, ux_s, uy_s, un_s;

				q_s[0] = q->ttt[s];
				q_s[1] = q->ttx[s];
				q_s[2] = q->tty[s];
				q_s[3] = q->ttn[s];
				q_s[4] = q->pl[s];

			#if (PT_MATCHING == 1)
				q_s[5] = q->pt[s];
			#endif

			#ifdef PIMUNU
				q_s[6] = q->pitt[s];
				q_s[7] = q->pitx[s];
				q_s[8] = q->pity[s];
				q_s[9] = q->pitn[s];
			#endif

			#ifdef W_TZ_MU
				q_s[10] = q->WtTz[s];
				q_s[11] = q->WxTz[s];
				q_s[12] = q->WyTz[s];
				q_s[13] = q->WnTz[s];
			#endif

				// compute the updated values for (e, p, u^mu)
				get_inferred_variables_test(t, q_s, &e_s, &ut_s, &ux_s, &uy_s, &un_s);




				e[s]     = e_s;
				u->ut[s] = ut_s;
				u->ux[s] = ux_s;
				u->uy[s] = uy_s;
				u->un[s] = un_s;
			}
		}
	}
}
*/



void set_inferred_variables(const CONSERVED_VARIABLES * const __restrict__ q, precision * const __restrict__ e, FLUID_VELOCITY * const __restrict__ u, precision t, int nx, int ny, int nz)
{
	// loop over the physical grid points
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
				precision pitt = 0.0;
				precision pitx = 0.0;
				precision pity = 0.0;
				precision pitn = 0.0;
			#endif

			#ifdef WTZMU
				precision WtTz = q->WtTz[s];
				precision WxTz = q->WxTz[s];
				precision WyTz = q->WyTz[s];
				precision WnTz = q->WnTz[s];
			#else
				precision WtTz = 0.0;
				precision WxTz = 0.0;
				precision WyTz = 0.0;
				precision WnTz = 0.0;
			#endif		

				precision Kt = ttt - pitt;				// [fm^-4]
				precision Kx = ttx - pitx;				// [fm^-4]
				precision Ky = tty - pity;				// [fm^-4]
				precision Kn = ttn - pitn;				// [fm^-5]

				precision A = Kn / (Kt + pl);		  	// [fm^-1] (check these formulas)
				precision B = WtTz / (t * (Kt + pl));		// [fm^-1]

				precision t2 = t * t;

				precision t2A2 = t2 * A * A;
				precision t2B2 = t2 * B * B;

				// figure out how to deal with this
				precision F = (A  -  B * sqrt(fabs(1.0 + t2B2 - t2A2))) / (1.0 + t2B2);	// [fm^-1]
				precision Ft = t * F;													// [1]
				precision x = sqrt(fabs(1.0  -  Ft * Ft));								// [1]

				precision zt = Ft / x;					// [1]
				precision zn = 1.0 / (x * t);			// [fm^-1]

				precision Mt = Kt  -  2.0 * WtTz * zt;
				precision Mx = Kx  -  WxTz * zt;
				precision My = Ky  -  WyTz * zt;
				precision Mn = Kn  -  (WtTz * zn  +  WnTz * zt);

				
				// solution for e
			#if (PT_MATCHING == 1)
				precision Ltt = (pl - pt) * zt * zt;

				precision ut_numerator = Mt  +  pt  -  Ltt;

				precision e_s = fmax(E_MIN, Mt  -  Ltt  -  (Mx * Mx  +  My * My) / ut_numerator  -  t2 * Mn * Mn * ut_numerator / (Mt + pl) / (Mt + pl));
			#else
				precision zt2  = zt * zt;
				precision ztzn = zt * zn;


				precision a = 0.25 * t2 * ztzn * ztzn  +  0.5 * (1.0 + zt2) * (1.0  -  0.5 * zt2);

				precision b = 0.5 * t2 * Mn * ztzn  -  1.5 * t2 * pl * ztzn * ztzn  -  (0.5 * (1.0 + zt2) * (Mt - 1.5 * pl * zt2)  -  (1.0 - 0.5 * zt2) * (Mt - 0.5 * pl * (1.0 + 3.0 * zt2)));

				precision c = Mx * Mx  +  My * My  +  t2 * Mn * Mn  -  1.5 * t2 * Mn * pl * ztzn  +  2.25 * t2 * pl * pl * ztzn * ztzn  -  (Mt - 0.5 * pl * (1.0 + 3.0 * zt2)) * (Mt - 1.5 * pl * zt2);


				// quadratic formula
				precision e_s = fmax(E_MIN, ( - b  +  sqrt(fabs(b * b  -  4.0 * a * c))) / (2.0 * a));
				//


				precision pt = 0.5 * (e_s - pl);

				precision ut_numerator = Mt  +  pt  -  (pl - pt) * zt2;

			#endif

				// solution for u^mu
				precision ut_s = sqrt(fabs(ut_numerator / (e_s + pt)));
				precision ux_s = Mx / ut_s / (e_s + pt);
				precision uy_s = My / ut_s / (e_s + pt);
				precision un_s = F * ut_s;

				if(std::isnan(ut_s))
				{
					printf("\nget_inferred_variables_test error: u^mu = (%lf, %lf, %lf, %lf) is nan\n", ut_s, ux_s, uy_s, un_s);
					exit(-1);
				}	

				ut_s = sqrt(1.0  +  ux_s * ux_s  +  uy_s * uy_s  +  t2 * un_s * un_s);

				
				// get solution for primary variables
				e[s]     = e_s;
				u->ut[s] = ut_s;
				u->ux[s] = ux_s;
				u->uy[s] = uy_s;
				u->un[s] = un_s;

				//test_energy_momentum_tensor_error(q, e_s, ut_s, ux_s, uy_s, un_s, t);

				test_energy_momentum_tensor_error_new(ttt, ttx, tty, ttn, pl, pt, pitt, pitx, pity, pitn, WtTz, WxTz, WyTz, WnTz, e_s, ut_s, ux_s, uy_s, un_s, t);
			}
		}
	}
}

