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





