#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include "../include/Macros.h"
#include "../include/Hydrodynamics.h"
#include "../include/TransportAnisoNonconformal.h"
#include "../include/TransportAniso.h"
#include "../include/Precision.h"

const precision g = 51.4103536012791;				// quasiparticle degeneracy factor (Eq. 46)
												// isn't this related to conformal_prefactor?
												// g = (2.(Nc.Nc - 1) + 4.Nc.Nf.7/8) . pi^4 / 90

// put the arrays here for now (can I get away with 16 points?)
const int pbar_pts = 32;

const precision pbar_root_a0[pbar_pts] = {0.044489365833267,0.234526109519619,0.576884629301886,1.07244875381782,1.72240877644465,2.52833670642579,3.49221327302199,4.61645676974977,5.90395850417424,7.35812673318624,8.9829409242126,10.78301863254,12.7636979867427,14.9311397555226,17.2924543367153,19.8558609403361,22.6308890131968,25.6286360224592,28.8621018163235,32.3466291539647,36.100494805752,40.1457197715394,44.5092079957549,49.2243949873086,54.3337213333969,59.892509162134,65.975377287935,72.6876280906627,80.1874469779135,88.7353404178924,98.829542868284,111.751398097938};

const precision pbar_weight_a0[pbar_pts] = {0.109218341952385,0.210443107938813,0.235213229669848,0.195903335972881,0.129983786286072,0.0705786238657174,0.0317609125091751,0.0119182148348386,0.00373881629461152,0.000980803306614955,0.000214864918801364,3.92034196798795e-05,5.93454161286863e-06,7.41640457866755e-07,7.60456787912078e-08,6.35060222662581e-09,4.28138297104093e-10,2.30589949189134e-11,9.79937928872709e-13,3.23780165772927e-14,8.17182344342072e-16,1.54213383339382e-17,2.11979229016362e-19,2.05442967378805e-21,1.3469825866374e-23,5.66129413039736e-26,1.41856054546304e-28,1.91337549445422e-31,1.19224876009822e-34,2.67151121924014e-38,1.33861694210626e-42,4.51053619389897e-48};

const precision pbar_root_a2[pbar_pts] = {0.196943922146667,0.529487866050161,1.01026981913845,1.640616191672,2.42200673335506,3.35625823737525,4.44557319147359,5.69257570606939,7.10035048878373,8.67248915845674,10.413146435518,12.3271087558129,14.4198784243951,16.6977773650005,19.1680758788069,21.839153763432,24.7207039368187,27.823992811746,31.1621978174102,34.7508519173206,38.6084399084037,42.7572156420076,47.2243504952188,52.0435960848824,57.257778984273,62.9227106235616,69.1136582681551,75.9368320953467,83.5517824825995,92.221284870548,102.447989923982,115.52490220024};

const precision pbar_weight_a2[pbar_pts] = {0.00825033790777967,0.0671033262747106,0.206386098255352,0.368179392999486,0.446389764546666,0.397211321904435,0.270703020914857,0.144937243765141,0.0619302157291065,0.0213227539141068,0.00594841159169929,0.00134795257769464,0.000248166548996264,3.7053223540482e-05,4.47057760459712e-06,4.33555258401213e-07,3.35571417159735e-08,2.05432200435071e-09,9.83646900727572e-11,3.63364388210833e-12,1.01834576904109e-13,2.12110313498633e-15,3.20100105319804e-17,3.39007439648141e-19,2.41904571899768e-21,1.10270714408855e-23,2.98827103874582e-26,4.34972188455989e-29,2.92108431650778e-32,7.0533942409897e-36,3.81617106981223e-40,1.39864930768275e-45};


aniso_transport_coefficients_nonconformal::aniso_transport_coefficients_nonconformal()
{

}


aniso_transport_coefficients_nonconformal::~aniso_transport_coefficients_nonconformal()
{

}


void aniso_transport_coefficients_nonconformal::compute_hypergeometric_functions_n_equals_0(precision z)
{
	if(z > delta)
	{
		precision sqrtz = sqrt(z);
		precision t = atan(sqrtz) / sqrtz;

		t_020 = 2. * (1. - t) / z;
		t_001 = 2. * (-1. / (1. + z)  +  t) / z;
	}
	else if(z < -delta && z > -1.)
	{
		precision sqrtmz = sqrt(-z);
		precision t = atanh(sqrtmz) / sqrtmz;

		t_020 = 2. * (1. - t) / z;
		t_001 = 2. * (-1. / (1. + z)  +  t) / z;
	}
	else if(fabs(z) <= delta)
	{
		precision z2 = z  * z;
		precision z3 = z2 * z;
		precision z4 = z3 * z;
		precision z5 = z4 * z;
		precision z6 = z5 * z;

		// double check these
		t_020 = 0.6666666666666666 - 0.4*z + 0.2857142857142857*z2 - 0.2222222222222222*z3 + 0.18181818181818182*z4 - 0.15384615384615385*z5 + 0.13333333333333333*z6;

		t_001 = 1.3333333333333335 - 1.6*z + 1.7142857142857144*z2 - 1.7777777777777777*z3 + 1.8181818181818181*z4 - 1.8461538461538463*z5 + 1.8666666666666667*z6;
	}
	else
	{
		printf("Error: z = %lf is out of bounds\n", z);
		exit(-1);
	}
}


void aniso_transport_coefficients_nonconformal::compute_hypergeometric_functions_n_equals_2(precision z)
{
	precision z2 = z  * z;

	// which ones do I need directly for the transport coefficients? not the root search
	// what about the quasiparticle ones? I forget

	if(z > delta)
	{
		precision sqrtz = sqrt(z);
		precision t = atan(sqrtz) / sqrtz;

		t_240 = (3.  +  2. * z  -  3. * (1. + z) * t) / z2;
		t_221 = (-3.  +  (3. + z) * t) / z2;
		t_202 = (3.  +  z  +  (z - 3.) * (1. + z) * t) / (z2 * (1. + z));
	}
	else if(z < -delta && z > -1.)
	{
		precision sqrtmz = sqrt(-z);
		precision t = atanh(sqrtmz) / sqrtmz;

		t_240 = (3.  +  2. * z  -  3. * (1. + z) * t) / z2;
		t_221 = (-3.  +  (3. + z) * t) / z2;
		t_202 = (3.  +  z  +  (z - 3.) * (1. + z) * t) / (z2 * (1. + z));
	}
	else if(fabs(z) <= delta)
	{
		precision z3 = z2 * z;
		precision z4 = z3 * z;
		precision z5 = z4 * z;
		precision z6 = z5 * z;

		t_240 = 0.4 - 0.17142857142857149*z + 0.09523809523809523*z2 - 0.06060606060606058*z3 + 0.04195804195804195*z4 - 0.030769230769230785*z5 + 0.023529411764705882*z6;

		t_221 = 0.2666666666666668 - 0.22857142857142854*z + 0.19047619047619047*z2 - 0.1616161616161616*z3 + 0.13986013986013987*z4 - 0.12307692307692308*z5 + 0.10980392156862744*z6;

		t_202 = 1.0666666666666664 - 1.3714285714285712*z + 1.5238095238095237*z2 - 1.616161616161616*z3 + 1.6783216783216781*z4 - 1.7230769230769227*z5 + 1.756862745098039*z6;
	}
	else
	{
		printf("Error: z = %lf is out of bounds\n", z);
		exit(-1);
	}
}


void aniso_transport_coefficients_nonconformal::compute_transport_coefficients(precision e, precision pl, precision pt, precision B, precision lambda, precision aT, precision aL, precision mbar, precision mass, precision mdmde)
{
	// can I cross-check with the conformal case by setting mbar = 0?
	precision lambda2 = lambda * lambda;
	precision lambda4 = lambda2 * lambda2;
	precision aT2 = aT * aT;
	precision aL2 = aL * aL;

	precision prefactor_n_equals_0 = g * aT2 * aL * lambda2 / (4. * M_PI * M_PI);
	precision prefactor_n_equals_2 = g * aT2 * aL * lambda4 / (4. * M_PI * M_PI);

	precision mbar2 = mbar * mbar;
	precision aT2_minus_aL2 = aT2 - aL2;

	// integration variables
	precision pbar, weight, pbar2, Ebar, exp_factor, total_weight;
	precision w, w2, w3, z;



	// are there any moments in anisotropic variable search that I can recycle??


	precision I_020 = 0;								// n = 0 moments
	precision I_001 = 0;

	precision I_240 = 0;								// n = 2 moments
	precision I_221 = 0;
	precision I_202 = 0;

	for(int i = 0; i < pbar_pts; i++)					// gauss integration loop by hand (can reduce time this way)
	{
		// compute n = 0 moments
		//::::::::::::::::::::::::::::::::::::::::::::
		pbar = pbar_root_a0[i];							// pbar roots / weights for a = n = 0
		weight = pbar_weight_a0[i];

		pbar2 = pbar * pbar;
		Ebar = sqrt(pbar2 + mbar2);
		exp_factor = exp(pbar - Ebar);
		total_weight = pbar * weight * exp_factor; 		// common weight for a = n = 0

		w = sqrt(aL2  +  mbar2 / pbar2);
		w2 = w * w;
		w3 = w2 * w;
		z = aT2_minus_aL2 / w2;

		compute_hypergeometric_functions_n_equals_0(z);

		I_020 += total_weight * t_020 / w3;				// R_020 = t_020 / w3, etc
		I_001 += total_weight * t_001 / w3;
		//::::::::::::::::::::::::::::::::::::::::::::


		// compute n = 2 moments
		//::::::::::::::::::::::::::::::::::::::::::::
		pbar = pbar_root_a2[i];							// pbar roots / weights for a = n = 2
		weight = pbar_weight_a2[i];

		pbar2 = pbar * pbar;
		Ebar = sqrt(pbar2 + mbar2);
		exp_factor = exp(pbar - Ebar);
		total_weight = pbar * weight * exp_factor; 		// common weight for a = n = 2

		w = sqrt(aL2  +  mbar2 / pbar2);
		w2 = w * w;
		w3 = w2 * w;
		z = aT2_minus_aL2 / w2;

		compute_hypergeometric_functions_n_equals_2(z);

		I_240 += total_weight * t_240 / w3;
		I_221 += total_weight * t_221 / w3;
		I_202 += total_weight * t_202 / w3;
		//::::::::::::::::::::::::::::::::::::::::::::

		// then repeat for n = 4
	}

	I_020 *= aL2 * prefactor_n_equals_0;
	I_001 *= aT2 * prefactor_n_equals_0 / 2.;

	I_240 *= aL2 * aL2 * prefactor_n_equals_2;
	I_221 *= aT2 * aL2 * prefactor_n_equals_2 / 2.;
	I_202 *= aT2 * aT2 * prefactor_n_equals_2 / 8.;

	// printf("I_020 = %lf\n", I_020);
	// printf("I_001 = %lf\n", I_001);
	// printf("I_240 = %lf\n", I_240);
	// printf("I_221 = %lf\n", I_221);
	//exit(-1);

// #if (NUMBER_OF_RESIDUAL_CURRENTS != 0)
// 	precision I_421 = prefactor * t_421 * Lambda2 * aL2 * 10.;
// 	precision I_402 = prefactor * t_402 * Lambda2 * 2.5;
// 	precision I_441 = prefactor * t_441 * Lambda2 * aL2 * 10.;
// 	precision I_422 = prefactor * t_422 * Lambda2 * 2.5;
// 	precision I_403 = prefactor * t_403 * Lambda2 * 5. / (12. * aL2);
// #endif

	// pl transport coefficients (there's also a quasiparticle contribution)

	//zeta_LL = I_240  -  3. * pl;
	zeta_LL = I_240  -  3. * (pl + B)  +  mdmde * (e + pl) * (I_020  +  (2.*pt + pl - e + 4.*B) / (mass * mass));
	zeta_TL = I_221 - pl;

/*
#ifdef WTZMU
	lambda_WuL = I_441 / I_421;
	lambda_WTL = 1. - lambda_WuL;
#else
	lambda_WuL = 0;
	lambda_WTL = 0;
#endif
#ifdef PIMUNU
	lambda_piL = I_422 / I_402;
#else
	lambda_piL = 0;
#endif
*/
	// pt transport coefficients
	//zeta_LT = I_221 - pt;
	zeta_LT = I_221  -  (pt + B)  +  mdmde * (e + pl) * (I_001  +  (2.*pt + pl - e + 4.*B) / (mass * mass));
	zeta_TT = 2. * (I_202 - pt);

	// printf("zeta_LL = %lf\n", zeta_LL);
	// printf("zeta_LT = %lf\n", zeta_LT);
	// exit(-1);

/*
#ifdef WTZMU
	lambda_WTT = 2. * I_422 / I_421;
	lambda_WuT = lambda_WTT  -  1.;
#else
	lambda_WTT = 0;
	lambda_WuT = 0;
#endif

#ifdef PIMUNU
	lambda_piT = 1.  -  3. * I_403 / I_402;
#else
	lambda_piT = 0;
#endif

	// WTz transport coefficients
#ifdef WTZMU
	eta_uW = 0.5 * (pl - I_221);
	eta_TW = 0.5 * (pt - I_221);
	tau_zW = pl - pt;
	lambda_WTW = 2. * I_422 / I_421  -  1.;
	lambda_WuW = 2.  -  I_441 / I_421;
	delta_WW = lambda_WTW  -  0.5;
#ifdef PIMUNU
	lambda_piuW = I_422 / I_402;
	lambda_piTW = lambda_piuW  -  1.;
#else
	lambda_piuW = 0;
	lambda_piTW = 0;
#endif
#endif
	// piT transport coefficients (validated 1st four are same as v1)
#ifdef PIMUNU
	eta_T 		= pt  -  I_202;
	tau_pipi 	= 2.  -  4. * I_403 / I_402;
	delta_pipi	= (3. * tau_pipi  + 2.) / 4.;
	lambda_pipi = I_422 / I_402  -  1.;
#ifdef WTZMU
	lambda_Wupi = lambda_WTW  -  1.;
	lambda_WTpi = lambda_WuW  +  2.;
#else
	lambda_Wupi = 0;
	lambda_WTpi = 0;
#endif
#endif

	//test_kinetic_solution(e, pl, pt, z, aL2, prefactor);

*/
}





