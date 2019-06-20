/*
 * EnergyMomentumTensor.cpp
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include "../include/InferredVariables.h"
#include "../include/DynamicalVariables.h"
#include "../include/Parameters.h"
#include "../include/EquationOfState.h"
//#include "../include/TransportCoefficients.h"

inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


// this is a temporary test for conformal PL matching
void get_inferred_variables_test(PRECISION t, const PRECISION * const __restrict__ q, PRECISION * const __restrict__ e, PRECISION * const __restrict__ ut, PRECISION * const __restrict__ ux, PRECISION * const __restrict__ uy, PRECISION * const __restrict__ un)
{
	PRECISION Ttt = q[0];
	PRECISION Ttx = q[1];
	PRECISION Tty = q[2];

	PRECISION pl  = q[4];

	PRECISION Mt = Ttt;
	PRECISION Mx = Ttx;
	PRECISION My = Tty;

	// conformal solution for (e, p)
	PRECISION e_s = 0.5 * (-(Mt - pl)  +  sqrt(fabs((Mt - pl) * (Mt - pl)  +  8.0 * (Mt * (Mt - 0.5 * pl) - (Mx * Mx  +  My * My)))));

	if(e_s < 1.e-3) e_s = 1.e-3;				// this cutoff is important

	PRECISION pt = 0.5 * (e_s - pl);			// conformal formula

	// solution for u^mu
	PRECISION ut_s = sqrt(fabs((Mt + pt) / (e_s + pt)));
	PRECISION ux_s = Mx / ut_s / (e_s + pt);
	PRECISION uy_s = My / ut_s / (e_s + pt);
	PRECISION un_s = 0.0;

	if(std::isnan(ut_s))	// I'm not sure what's going nan???
	{
		printf("\ngetInferredVariables error: u^mu = (%lf, %lf, %lf, %lf) is nan\n", ut_s, ux_s, uy_s, un_s);
		exit(-1);
	}

	// get solution for primary variables
	*e = e_s;
	*ut = sqrt(1.0  +  ux_s * ux_s  +  uy_s * uy_s);	// I guess renormalize helps?
	//*ut = ut_s;
	*ux = ux_s;
	*uy = uy_s;
	*un = un_s;
}





void get_inferred_variables(PRECISION t, const PRECISION * const __restrict__ q, PRECISION * const __restrict__ e, PRECISION * const __restrict__ ut, PRECISION * const __restrict__ ux, PRECISION * const __restrict__ uy, PRECISION * const __restrict__ un)
{
	// this is what my version looks like
	PRECISION Ttt = q[0];
	PRECISION Ttx = q[1];
	PRECISION Tty = q[2];
	PRECISION Ttn = q[3];
	PRECISION Pl = q[4];

	PRECISION Pt = 0.0;		// the conformal eos depends on e which is a problem

	PRECISION pitt = 0.0;
	PRECISION pitx = 0.0;
	PRECISION pity = 0.0;
	PRECISION pitn = 0.0;

	PRECISION Wt = 0.0;
	PRECISION Wx = 0.0;
	PRECISION Wy = 0.0;
	PRECISION Wn = 0.0;

#ifdef PIMUNU
	pitt = q[6];
	pitx = q[7];
	pity = q[8];
	pitn = q[9];
#endif

#ifdef W_TZ_MU
	Wt = q[10];  // I only need to pass an array so make sure it's the right elements
	Wx = q[11];
	Wy = q[12];
	Wn = q[13];
#endif

	// I can use Wt = t2 * Wn * zn / zt
	// well it also shows up in B (so maybe I can't)

	PRECISION Kt = Ttt - pitt;				// [fm^-4]
	PRECISION Kx = Ttx - pitx;				// [fm^-4]
	PRECISION Ky = Tty - pity;				// [fm^-4]
	PRECISION Kn = Ttn - pitn;				// [fm^-5]

	PRECISION A = Kn / (Kt + Pl);		  	// [fm^-1] (check these formulas)
	PRECISION B = Wt / (t * (Kt + Pl));		// [fm^-1]

	PRECISION t2 = t * t;

	PRECISION t2A2 = t2 * A * A;
	PRECISION t2B2 = t2 * B * B;

	// figure out how to deal with this
	PRECISION F = (A  -  B * sqrt(1.0 + t2B2 - t2A2)) / (1.0 + t2B2);	// [fm^-1]
	PRECISION Ft = t * F;												// [1]
	PRECISION x = sqrt(1.0  -  Ft * Ft);								// [1]

	if(std::isnan(F) || std::isnan(x))
	{
		printf("\ngetInferredVariables error: (F, x) = (%lf, %lf) is nan\n", F, x);
		//exit(-1);
	}

	PRECISION zt = Ft / x;					// [1]
	PRECISION zn = 1.0 / (x * t);			// [fm^-1]

	PRECISION Mt = Kt  -  2.0 * Wt * zt;
	PRECISION Mx = Kx  -  Wx * zt;
	PRECISION My = Ky  -  Wy * zt;
	PRECISION Mn = Kn  -  (Wt * zn  +  Wn * zt);

	PRECISION dP_zt2 = (Pl - Pt) * zt * zt;
	PRECISION ut_numerator = Mt + Pt - dP_zt2;


	// solution for (e, p)
	PRECISION e_s = Mt  -  dP_zt2  -  (Mx * Mx  +  My * My) / ut_numerator  -  t2 * Mn * Mn * ut_numerator / (Mt + Pl) / (Mt + Pl);

	if(e_s < 0.0)
	{
		printf("\ngetInferredVariables error: e = %lf is negative\n", e_s);
	}

	PRECISION e_plus_Pt = e_s + Pt;

	// solution for u^mu
	PRECISION ut_s = sqrt(ut_numerator / e_plus_Pt);
	PRECISION ux_s = Mx / ut_s / e_plus_Pt;
	PRECISION uy_s = My / ut_s / e_plus_Pt;
	PRECISION un_s = F * ut_s;

	if(std::isnan(ut_s) || std::isnan(ux_s) || std::isnan(uy_s) || std::isnan(un_s))
	{
		printf("\ngetInferredVariables error: u^mu = (%lf, %lf, %lf, %lf) is nan\n", ut_s, ux_s, uy_s, un_s);
	}

	// test the reconstruction of T^{tau mu}
	// double dTtt = e_plus_Pt * ut_s * ut_s  -  Pt  +  dP_zt2  +  2.0 * Wt * zt  +  pitt  - Ttt;
	// double dTtx = e_plus_Pt * ut_s * ux_s  +  Wx * zt  + pitx  -  Ttx;
	// double dTty = e_plus_Pt * ut_s * uy_s  +  Wy * zt  + pity  -  Tty;
	// double dTtn = e_plus_Pt * ut_s * un_s  +  (Pl - Pt) * zt * zn  +  Wt * zn  +  Wn * zt  +  pitn  -  Ttn;

	// double eps = 1.e-14;
	// if(fabs(dTtt) > eps || fabs(dTtx) > eps || fabs(dTty) > eps || fabs(dTtn) > eps)
	// {
	// 	printf("\ngetInferredVariables error: dT^{tau/mu} = (%lf, %lf, %lf, %lf) > 1.e-14\n", dTtt, dTtx, dTty, dTtn);
	// }

	// get solution for primary variables
	*e = e_s;
	*ut = sqrt(1.0  +  ux_s * ux_s  +  uy_s * uy_s  +  t2 * un_s * un_s);
	*ux = ux_s;
	*uy = uy_s;
	*un = un_s;
}



void set_inferred_variables(const CONSERVED_VARIABLES * const __restrict__ q, PRECISION * const __restrict__ e, FLUID_VELOCITY * const __restrict__ u, PRECISION t, int nx, int ny, int nz)
{
	int N = 6;  // default size for q_s array (T^\tau\mu, Pl, Pt)

// #ifdef PIMUNU
// 	N += 4;		// add time components of pimunu and Wmu
// #endif

// #ifdef W_TZ_MU
// 	N += 4;
// #endif
	const int nq = N;
	PRECISION q_s[nq];

	// loop over the physical grid points
	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				PRECISION e_s, ut_s, ux_s, uy_s, un_s;		// do I really need p_s?

				q_s[0] = q->ttt[s];
				q_s[1] = q->ttx[s];
				q_s[2] = q->tty[s];
				q_s[3] = q->ttn[s];
				q_s[4] = q->pl[s];
				//q_s[5] = q->pt[s];  // need to add pt to the struct

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
				//
				// test version for conformal PL matching
				get_inferred_variables_test(t, q_s, &e_s, &ut_s, &ux_s, &uy_s, &un_s);

				// set the updated values to current variables
				e[s] = e_s;

				u->ut[s] = ut_s;	// should I renormalize this?
				u->ux[s] = ux_s;
				u->uy[s] = uy_s;
				u->un[s] = un_s;
			}
		}
	}
}

