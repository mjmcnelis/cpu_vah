#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include "../include/Macros.h"
#include "../include/Hydrodynamics.h"
#include "../include/TransportAnisoNonconformal.h"
#include "../include/TransportAniso.h"
#include "../include/DynamicalVariables.h"
#include "../include/EquationOfState.h"
#include "../include/Precision.h"
#include "../include/Parameters.h"


const double g = 51.4103536012791;				// quasiparticle degeneracy factor (Eq. 46)
												// isn't this related to conformal_prefactor?
												// g = (2.(Nc.Nc - 1) + 4.Nc.Nf.7/8) . pi^4 / 90


aniso_transport_coefficients_nonconformal::aniso_transport_coefficients_nonconformal()
{

}


aniso_transport_coefficients_nonconformal::~aniso_transport_coefficients_nonconformal()
{

}


void aniso_transport_coefficients_nonconformal::compute_hypergeometric_functions_n_equals_2(precision z)
{
	precision z2 = z  * z;
	precision z3 = z2 * z;

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


void aniso_transport_coefficients_nonconformal::compute_transport_coefficients(precision pl, precision pt, precision lambda, precision aT, precision aL, precision mbar)
{
	// can I cross-check with the conformal case by setting mbar = 0?

	precision lambda2 = lambda * lambda;
	precision lambda4 = lambda2 * lambda2;
	precision aT2 = aT * aT;
	precision aL2 = aL * aL;
	precision common_prefactor = g * aT2 * aL * lambda4 / (4. * M_PI * M_PI);

	precision mbar2 = mbar * mbar;
	precision aT2_minus_aL2 = aT2 - aL2;

	// integration variables
	precision pbar, weight, pbar2, Ebar, exp_factor, total_weight;
	precision w, w2, w3, z;

	// put the arrays here for now (can I get away with 16 points?)
	const int pbar_pts = 32;

	double pbar_root_a2[pbar_pts] = {0.196943922146667,0.529487866050161,1.01026981913845,1.640616191672,2.42200673335506,3.35625823737525,4.44557319147359,5.69257570606939,7.10035048878373,8.67248915845674,10.413146435518,12.3271087558129,14.4198784243951,16.6977773650005,19.1680758788069,21.839153763432,24.7207039368187,27.823992811746,31.1621978174102,34.7508519173206,38.6084399084037,42.7572156420076,47.2243504952188,52.0435960848824,57.257778984273,62.9227106235616,69.1136582681551,75.9368320953467,83.5517824825995,92.221284870548,102.447989923982,115.52490220024};

	double pbar_weight_a2[pbar_pts] = {0.00825033790777967,0.0671033262747106,0.206386098255352,0.368179392999486,0.446389764546666,0.397211321904435,0.270703020914857,0.144937243765141,0.0619302157291065,0.0213227539141068,0.00594841159169929,0.00134795257769464,0.000248166548996264,3.7053223540482e-05,4.47057760459712e-06,4.33555258401213e-07,3.35571417159735e-08,2.05432200435071e-09,9.83646900727572e-11,3.63364388210833e-12,1.01834576904109e-13,2.12110313498633e-15,3.20100105319804e-17,3.39007439648141e-19,2.41904571899768e-21,1.10270714408855e-23,2.98827103874582e-26,4.34972188455989e-29,2.92108431650778e-32,7.0533942409897e-36,3.81617106981223e-40,1.39864930768275e-45};

	precision I_240 = 0;
	precision I_221 = 0;
	precision I_202 = 0;

	for(int i = 0; i < pbar_pts; i++)					// gauss integration loop by hand (can reduce time this way)
	{
		// compute n = 2 moments
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


		// then repeat for n = 4

	}


	I_240 *= aL2 * aL2 * common_prefactor;
	I_221 *= aT2 * aL2 * common_prefactor / 2.;
	I_202 *= aT2 * aT2 * common_prefactor / 8.;

	printf("I_240 = %lf\n", I_240);
	printf("I_221 = %lf\n", I_221);
	exit(-1);


// #if (NUMBER_OF_RESIDUAL_CURRENTS != 0)
// 	precision I_421 = prefactor * t_421 * Lambda2 * aL2 * 10.;
// 	precision I_402 = prefactor * t_402 * Lambda2 * 2.5;
// 	precision I_441 = prefactor * t_441 * Lambda2 * aL2 * 10.;
// 	precision I_422 = prefactor * t_422 * Lambda2 * 2.5;
// 	precision I_403 = prefactor * t_403 * Lambda2 * 5. / (12. * aL2);
// #endif

	// pl transport coefficients (there's also a quasiparticle contribution)
	zeta_LL = I_240  -  3. * pl;
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
	zeta_LT = I_221  -  pt;
	zeta_TT = 2. * (I_202 - pt);

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





