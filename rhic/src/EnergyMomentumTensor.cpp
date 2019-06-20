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
#include "../include/EnergyMomentumTensor.h"
#include "../include/DynamicalVariables.h"
#include "../include/Parameters.h"
#include "../include/EquationOfState.h"
#include "../include/AnisotropicDistributionFunctions.h"

#define MAX_ITERS 10000
#define VBAR 0.563624     // used for what?
#define EPS 0.1


inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}

PRECISION getTransverseFluidVelocityMagnitude(const FLUID_VELOCITY * const __restrict__ u, int s)
 {
		PRECISION u1 = u->ux[s];
		PRECISION u2 = u->uy[s];
		return sqrt(fabs(u1 * u1 + u2 * u2));
}

int transverseFluidVelocityFromConservedVariables(PRECISION t, PRECISION ePrev, PRECISION uT_0,
PRECISION MB0, PRECISION MBT, PRECISION MB3, PRECISION PL, PRECISION Pi, double Ft, double x, double *uT,
int i, int jj, int k, double xi, double yj, double zk,
int fullTimeStepInversion
) {
	PRECISION uT0 = uT_0;	// initial guess for uT

	// Constants
	double Ft2 = Ft*Ft;
	double bT = x*MBT;
	double bL = x*x*MB0-Ft2*PL;
	double b = x*x+Ft2;

	double f,fp,DF;

	for(int j = 0; j < MAX_ITERS; ++j)
	{
		double e = MB0 - t*Ft*MB3 - uT0/sqrt(1 + uT0*uT0)*x*MBT;
		if(e < 0.0) return -1;
		double p = equilibriumPressure(e);
		double PtHat = 0.5*(e-PL);
		double Pt = PtHat + 1.5*Pi;

		double deduT = -x*MBT/pow(1 + uT0*uT0,1.5);
		double dPtduT = 0.5*deduT;

		f = uT0/sqrt(1+uT0*uT0)*(bL+b*Pt) - bT;
		fp = 1/pow(1 + uT0*uT0,1.5)*(bL+b*Pt)+uT0/sqrt(1+uT0*uT0)*b*dPtduT;

		if(fabs(fp)==0.0) fp = 1.e-16;

		DF = f/fp;

		*uT = uT0 - DF;

		if(isnan(*uT) || isinf(*uT) || *uT < 0 || *uT > 9.0072e+15) return -1;

		double DUT = fabs(*uT-uT0);
		double UT = fabs(*uT);
		if(DUT <=  1.e-7 * UT) return 0;
		uT0 = *uT;
	}
	return -1;
}


void get_inferred_variables_test(PRECISION t, const PRECISION * const __restrict__ q, PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, PRECISION * const __restrict__ ut, PRECISION * const __restrict__ ux, PRECISION * const __restrict__ uy, PRECISION * const __restrict__ un)
{
	// this is a test for conformal PL matching
	PRECISION Ttt = q[0];
	PRECISION Ttx = q[1];
	PRECISION Tty = q[2];

	PRECISION pl  = q[4];

	PRECISION Mt = Ttt;
	PRECISION Mx = Ttx;
	PRECISION My = Tty;

	// solution for (e, p)
	PRECISION e_s = 0.5 * (-(Mt - pl)  +  sqrt(fabs((Mt - pl) * (Mt - pl)  +  8.0 * (Mt * (Mt - 0.5 * pl) - (Mx * Mx  +  My * My)))));
	
	if(e_s < 1.e-3) e_s = 1.e-3;				// this cutoff is important 
	PRECISION p_s = equilibriumPressure(e_s);

	PRECISION pt = 0.5 * (e_s - pl);			// conformal formula
	
	// solution for u^mu
	PRECISION ut_s = sqrt(fabs((Mt + pt) / (e_s + pt)));
	PRECISION ux_s = Mx / ut_s / (e_s + pt);
	PRECISION uy_s = My / ut_s / (e_s + pt);
	PRECISION un_s = 0.0;

	if(std::isnan(ut_s))
	{
		printf("\ngetInferredVariables error: u^mu = (%lf, %lf, %lf, %lf) is nan\n", ut_s, ux_s, uy_s, un_s);
		exit(-1);
	}

	// get solution for primary variables
	*e = e_s;
	*p = p_s;
	*ut = sqrt(1.0  +  ux_s * ux_s  +  uy_s * uy_s  +  t * t * un_s * un_s);	
	*ux = ux_s;
	*uy = uy_s;
	*un = un_s;
}





void get_inferred_variables_new(PRECISION t, const PRECISION * const __restrict__ q, PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, PRECISION * const __restrict__ ut, PRECISION * const __restrict__ ux, PRECISION * const __restrict__ uy, PRECISION * const __restrict__ un)
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
	PRECISION p_s = equilibriumPressure(e_s);

	if(e_s < 0.0 || p_s < 0.0)
	{
		printf("\ngetInferredVariables error: (e, p) = (%lf, %lf) is negative\n", e_s, p_s);
	}

	PRECISION e_plus_Pt = e_s + Pt;

	// solution for u^mu
	PRECISION ut_s = sqrt(ut_numerator / e_plus_Pt);
	PRECISION ux_s = Mx / ut_s / e_plus_Pt;
	PRECISION uy_s = My / ut_s / e_plus_Pt;
	PRECISION un_s = F * ut_s;

	if(std::isnan(ut_s))
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
	*p = p_s;
	*ut = sqrt(1.0  +  ux_s * ux_s  +  uy_s * uy_s  +  t2 * un_s * un_s);	// renormalize u
	*ux = ux_s;
	*uy = uy_s;
	*un = un_s;
}



int get_inferred_variables(PRECISION t, const PRECISION * const __restrict__ q, PRECISION ePrev, PRECISION uT_0,
PRECISION * const __restrict__ e, PRECISION * const __restrict__ p,
PRECISION * const __restrict__ ut, PRECISION * const __restrict__ ux, PRECISION * const __restrict__ uy, PRECISION * const __restrict__ un, int i, int j, int k, double xi, double yj, double zk, int fullTimeStepInversion)
{
	// so I need to pass q with the right array

	PRECISION ttt = q[0];
	PRECISION ttx = q[1];
	PRECISION tty = q[2];
	PRECISION ttn = q[3];
	PRECISION pl = q[4];
	// \pi^{\mu\nu}_{\perp}
#ifdef PIMUNU
	PRECISION pitt = q[5];
	PRECISION pitx = q[6];
	PRECISION pity = q[7];
	PRECISION pitn = q[8];
#else
	PRECISION pitt = 0;
	PRECISION pitx = 0;
	PRECISION pity = 0;
	PRECISION pitn = 0;
#endif
	// W^{\mu}_{\perp z}
#ifdef W_TZ_MU
	PRECISION WtTz = q[15];
	PRECISION WxTz = q[16];
	PRECISION WyTz = q[17];
	PRECISION WnTz = q[18];
#else
	PRECISION WtTz = 0;
	PRECISION WxTz = 0;
	PRECISION WyTz = 0;
	PRECISION WnTz = 0;
#endif
	PRECISION Pi = 0;

	PRECISION M0 = ttt - pitt;
	PRECISION M1 = ttx - pitx;
	PRECISION M2 = tty - pity;
	PRECISION M3 = ttn - pitn;

	double t2 = t * t;

	double M0PL = M0 + pl;
	if(M0PL == 0.0) M0PL = 1.e-16;

	PRECISION A = M3/M0PL;
	PRECISION At = t*A;
	double B = WtTz/M0PL/t;
	double Bt = t*B;
	double At2 = At*At;
	double Bt2 = Bt*Bt;
	double F = (A-fabs(B)*sqrt(fabs(1-At2+Bt2)))/(1+Bt2);
	double Ft = t*F;
	double Ft2 = Ft*Ft;
	PRECISION x = sqrt(fabs(1.-Ft2));

	double MB0 = M0-2*WtTz*Ft/x;
	double MB1 = M1-WxTz*Ft/x;
	double MB2 = M2-WyTz*Ft/x;
	double MB3 = M3-(1+Ft2)*WtTz/t/x;

	double MBT = sqrt(MB1*MB1+MB2*MB2);
	if(MBT==0.0) MBT=1.e-16;

	double uT;
	int status = -1;

	status = transverseFluidVelocityFromConservedVariables(t, ePrev, uT_0, MB0, MBT, MB3, pl, Pi, Ft, x, &uT, i, j, k, xi, yj, zk, fullTimeStepInversion);

	

	if(status == 0)
	{
		double C2 = 1.0+pow(uT,2.);		// this fails to work everywhere
		double C = sqrt(C2);
		double U = uT/C;

		*ux=uT*MB1/MBT;
		*uy=uT*MB2/MBT;
		*un = F*C/x;
		*ut = C/x;
		*e = MB0 - t*Ft*MB3 - U*x*MBT;
		*p = equilibriumPressure(*e);
	}
	else
	{
		*e = ePrev*.999;
		*p = equilibriumPressure(*e);
		*ux=0.0;
		*uy=0.0;
		*un = 0.0;
		*ut = 1.0;
	}


	if (isnan(*e) || isnan(*ut) || isnan(*ux) || isnan(*uy) || isnan(*un)) {
		printf("=======================================================================================\n");
		printf("found NaN in getInferredVariables.\n");
		printf("Grid point = (%d, %d, %d) = (%.3f, %.3f, %.3f)\n", i, j, k, xi, yj, zk);
		if(fullTimeStepInversion==0) printf("From semiDiscreteKurganovTadmorAlgorithm.\n");
		printf("t=%.3f\n",t);
		printf("uT=%.9f\n",uT);
		printf("ePrev=%.9f\n",ePrev);
		printf("A=%.9f;B=%.9f;F=%.9f;x=%.9f;\n",A,B,F,x);
		printf("e=%.9f;p=%.9f;\n",*e,*p);
		printf("ut=%.9f;ux=%.9f;uy=%.9f;un=%.9f;\n",*ut,*ux,*uy,*un);
		printf("MB1=%.9f\n",MB1);
		printf("MB2=%.9f\n",MB2);
		printf("MB0=%.3f,\t MBT=%.3f,\t MB3=%.3f,\tPL=%.3f,ePrev=%.3f\t,uT_0=%.3f\n", MB0, MBT, MB3,pl,ePrev,uT_0);
		printf("=======================================================================================\n");
		exit(-1);
	}


	// so this end game is to update (e,p,u^mu)

	return status;
}

void set_inferred_variables_new(const CONSERVED_VARIABLES * const __restrict__ q, PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, FLUID_VELOCITY * const __restrict__ u, PRECISION t, int nx, int ny, int nz, int ncx, int ncy)
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

	// loop over the physical grid points (maybe computing x, y, z is useful for debugging)
	for(int i = 2; i < nx + 2; i++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int k = 2; k < nz + 2; k++)
			{
				int s = linear_column_index(i, j, k, ncx, ncy);

				PRECISION e_s, p_s;
				PRECISION ut_s, ux_s, uy_s, un_s;

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
				get_inferred_variables_test(t, q_s, &e_s, &p_s, &ut_s, &ux_s, &uy_s, &un_s);

				// set the updated values to current variables
				e[s] = e_s;
				p[s] = p_s;

				u->ut[s] = ut_s;	// should I renormalize this?
				u->ux[s] = ux_s;
				u->uy[s] = uy_s;
				u->un[s] = un_s;
			}
		}
	}
}


void set_inferred_variables(const CONSERVED_VARIABLES * const __restrict__ q, PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, FLUID_VELOCITY * const __restrict__ up, FLUID_VELOCITY * const __restrict__ u,
PRECISION t, int nx, int ny, int nz, int ncx, int ncy, PRECISION *fTSolution)
{
	// loop over the physical grid points (edited on 6/7)
	for(int i = 2; i < nx + 2; i++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int k = 2; k < nz + 2; k++)
			{
				int s = linear_column_index(i, j, k, ncx, ncy);

				// this is initialized with up
				PRECISION uT = getTransverseFluidVelocityMagnitude(up, s);


				PRECISION q_s[NUMBER_CONSERVED_VARIABLES];
				PRECISION e_s, p_s;
				PRECISION ut_s, ux_s, uy_s, un_s;

				q_s[0] = q->ttt[s];
				q_s[1] = q->ttx[s];
				q_s[2] = q->tty[s];
				q_s[3] = q->ttn[s];
				q_s[4] = q->pl[s];
				// \pi^{\mu\nu}_{\perp}
// #ifdef PIMUNU
// 				q_s[5] = q->pitt[s];
// 				q_s[6] = q->pitx[s];
// 				q_s[7] = q->pity[s];
// 				q_s[8] = q->pitn[s];
// #endif
// 				// w^{\mu}_{\perp z}
// #ifdef W_TZ_MU
// 				q_s[15] = q->WtTz[s];
// 				q_s[16] = q->WxTz[s];
// 				q_s[17] = q->WyTz[s];
// 				q_s[18] = q->WnTz[s];
//#endif

				// compute the updated values for (e, p, u^mu)
				int status = get_inferred_variables(t, q_s, e[s], uT, &e_s, &p_s, &ut_s, &ux_s, &uy_s, &un_s, i, j, k, 0, 0, 0, 1);


				// debugging
				if(status == 0) fTSolution[s] = 0.0;
				else fTSolution[s] = 1.0;

				// set the updated values to current variables
				e[s] = e_s;
				p[s] = p_s;
				u->ut[s] = ut_s;
				u->ux[s] = ux_s;
				u->uy[s] = uy_s;
				u->un[s] = un_s;
			}
		}
	}
}

//===================================================================
// Components of T^{\mu\nu} in (\tau,x,y,\eta_s)-coordinates
//===================================================================   (what is this for?)
// PRECISION Ttt(PRECISION e, PRECISION p, PRECISION ut, PRECISION pitt)
// {
// 	return (e + p) * ut * ut  -  p  +  pitt;
// }

// PRECISION Ttx(PRECISION e, PRECISION p, PRECISION ut, PRECISION ux, PRECISION pitx)
// {
// 	return (e + p) * ut * ux  +  pitx;
// }

// PRECISION Tty(PRECISION e, PRECISION p, PRECISION ut, PRECISION uy, PRECISION pity)
// {
// 	return (e + p) * ut * uy  +  pity;
// }

// PRECISION Ttn(PRECISION e, PRECISION p, PRECISION ut, PRECISION un, PRECISION pitn) {
// 	return (e + p) * ut * un  +  pitn;
// }
