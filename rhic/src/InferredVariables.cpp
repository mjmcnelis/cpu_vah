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

inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


// this is a temporary test for conformal PL matching
void get_inferred_variables_test(precision t, const precision * const __restrict__ q, precision * const __restrict__ e, precision * const __restrict__ ut, precision * const __restrict__ ux, precision * const __restrict__ uy, precision * const __restrict__ un)
{
	precision ttt = q[0];
	precision ttx = q[1];
	precision tty = q[2];
	precision ttn = q[3];
	precision pl  = q[4];
	
	precision Mt = ttt;
	precision Mx = ttx;
	precision My = tty;
	precision Mn = ttn;

// reconstruction formula for energy density
#if (PT_MATCHING == 0)

	// this needs to be further generalized to include Mn
	precision e_s = 0.5 * (-(Mt - pl)  +  sqrt(fabs((Mt - pl) * (Mt - pl)  +  8.0 * (Mt * (Mt - 0.5 * pl) - (Mx * Mx  +  My * My)))));

	if(e_s < 1.e-3 || std::isnan(e_s)) e_s = 1.e-3;		// this cutoff is important

	precision pt = 0.5 * (e_s - pl);					// conformal formula

#else
	precision pt  = q[5];

	precision e_s = Mt  -  (Mx * Mx  +  My * My) / (Mt + pt);	// generalize 

	if(e_s < 1.e-3 || std::isnan(e_s)) e_s = 1.e-3;

#endif

	// fluid velocity (generalize)
	precision ut_s = sqrt(fabs((Mt + pt) / (e_s + pt)));
	precision ux_s = Mx / ut_s / (e_s + pt);
	precision uy_s = My / ut_s / (e_s + pt);
	precision un_s = 0.0;


	if(std::isnan(ut_s))	
	{	
		// I'm not sure what's going nan???
		printf("\ngetInferredVariables error: u^mu = (%lf, %lf, %lf, %lf) is nan\n", ut_s, ux_s, uy_s, un_s);
		exit(-1);
	}


	// get solution for primary variables
	*e  = e_s;
	*ut = sqrt(1.0  +  ux_s * ux_s  +  uy_s * uy_s  +  t * t * un_s * un_s);
	*ux = ux_s;
	*uy = uy_s;
	*un = un_s;
}





void get_inferred_variables(precision t, const precision * const __restrict__ q, precision * const __restrict__ e, precision * const __restrict__ ut, precision * const __restrict__ ux, precision * const __restrict__ uy, precision * const __restrict__ un)
{
	// this is what my version looks like
	precision Ttt = q[0];
	precision Ttx = q[1];
	precision Tty = q[2];
	precision Ttn = q[3];
	precision Pl = q[4];

	precision Pt = 0.0;		// the conformal eos depends on e which is a problem

	precision pitt = 0.0;
	precision pitx = 0.0;
	precision pity = 0.0;
	precision pitn = 0.0;

	precision Wt = 0.0;
	precision Wx = 0.0;
	precision Wy = 0.0;
	precision Wn = 0.0;

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

	precision Kt = Ttt - pitt;				// [fm^-4]
	precision Kx = Ttx - pitx;				// [fm^-4]
	precision Ky = Tty - pity;				// [fm^-4]
	precision Kn = Ttn - pitn;				// [fm^-5]

	precision A = Kn / (Kt + Pl);		  	// [fm^-1] (check these formulas)
	precision B = Wt / (t * (Kt + Pl));		// [fm^-1]

	precision t2 = t * t;

	precision t2A2 = t2 * A * A;
	precision t2B2 = t2 * B * B;

	// figure out how to deal with this
	precision F = (A  -  B * sqrt(1.0 + t2B2 - t2A2)) / (1.0 + t2B2);	// [fm^-1]
	precision Ft = t * F;												// [1]
	precision x = sqrt(1.0  -  Ft * Ft);								// [1]

	if(std::isnan(F) || std::isnan(x))
	{
		printf("\ngetInferredVariables error: (F, x) = (%lf, %lf) is nan\n", F, x);
		//exit(-1);
	}

	precision zt = Ft / x;					// [1]
	precision zn = 1.0 / (x * t);			// [fm^-1]

	precision Mt = Kt  -  2.0 * Wt * zt;
	precision Mx = Kx  -  Wx * zt;
	precision My = Ky  -  Wy * zt;
	precision Mn = Kn  -  (Wt * zn  +  Wn * zt);

	precision dP_zt2 = (Pl - Pt) * zt * zt;
	precision ut_numerator = Mt + Pt - dP_zt2;


	// solution for (e, p)
	precision e_s = Mt  -  dP_zt2  -  (Mx * Mx  +  My * My) / ut_numerator  -  t2 * Mn * Mn * ut_numerator / (Mt + Pl) / (Mt + Pl);

	if(e_s < 0.0)
	{
		printf("\ngetInferredVariables error: e = %lf is negative\n", e_s);
	}

	precision e_plus_Pt = e_s + Pt;

	// solution for u^mu
	precision ut_s = sqrt(ut_numerator / e_plus_Pt);
	precision ux_s = Mx / ut_s / e_plus_Pt;
	precision uy_s = My / ut_s / e_plus_Pt;
	precision un_s = F * ut_s;

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



void set_inferred_variables(const CONSERVED_VARIABLES * const __restrict__ q, precision * const __restrict__ e, FLUID_VELOCITY * const __restrict__ u, precision t, int nx, int ny, int nz)
{
	int N = 5;  // default size for q_s array (ttt, ttx, tty, ttn, pl)

#if (PT_MATCHING == 1)
	N += 1;
#endif
// #ifdef PIMUNU
// 	N += 4;		// add time components of pimunu and Wmu
// #endif

// #ifdef W_TZ_MU
// 	N += 4;
// #endif
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

