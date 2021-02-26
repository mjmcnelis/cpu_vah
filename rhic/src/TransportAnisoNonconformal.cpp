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
		printf("compute_hypergeometric_functions_n_equals_0 error: z = %lf is out of bounds\n", z);
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
		printf("compute_hypergeometric_functions_n_equals_2 error: z = %lf is out of bounds\n", z);
		exit(-1);
	}
}




void aniso_transport_coefficients_nonconformal::compute_hypergeometric_functions_n_equals_4(precision z)
{
	precision z2 = z  * z;
	precision z3 = z2 * z;

	if(z > delta)
	{
		precision sqrtz = sqrt(z);
		precision t = atan(sqrtz) / sqrtz;

		t_421 = (3.  +  z  +  (1. + z) * (z - 3.) * t) / (4. * z2);
		t_402 = (3. * (z - 1.)  +  (z * (3.*z - 2.) + 3.) * t) / (4. * z2);
		t_441 = (-15.  -  13. * z  +  3. * (1. + z) * (5. + z) * t) / (4. * z3);
		t_422 = (15.  +  z +  (z * (z - 6.) - 15.)* t) / (4. * z3);
		t_403 = ((z - 3.) * (5. + 3.*z)  +  3. * (1. + z) * (z * (z - 2.) + 5.) * t) / (4. * z3 * (1. + z));
	}
	else if(z < -delta && z > -1.)
	{
		precision sqrtmz = sqrt(-z);
		precision t = atanh(sqrtmz) / sqrtmz;

		t_421 = (3.  +  z  +  (1. + z) * (z - 3.) * t) / (4. * z2);
		t_402 = (3. * (z - 1.)  +  (z * (3.*z - 2.) + 3.) * t) / (4. * z2);
		t_441 = (-15.  -  13. * z  +  3. * (1. + z) * (5. + z) * t) / (4. * z3);
		t_422 = (15.  +  z +  (z * (z - 6.) - 15.)* t) / (4. * z3);
		t_403 = ((z - 3.) * (5. + 3.*z)  +  3. * (1. + z) * (z * (z - 2.) + 5.) * t) / (4. * z3 * (1. + z));
	}
	else if(fabs(z) <= delta)
	{
		precision z4 = z3 * z;
		precision z5 = z4 * z;
		precision z6 = z5 * z;

		t_421 = 0.2666666666666666 - 0.0761904761904762*z + 0.0380952380952381*z2 - 0.023088023088023088*z3 + 0.015540015540015537*z4 - 0.011188811188811189*z5 + 0.00844645550527904*z6;
   		t_402 = 1.0666666666666667 - 0.4571428571428572*z + 0.3047619047619048*z2 - 0.23088023088023088*z3 + 0.1864801864801865*z4 - 0.15664335664335666*z5 + 0.13514328808446457*z6;
	   	t_441 = 0.1142857142857145 - 0.07619047619047613*z + 0.051948051948051896*z2 - 0.037296037296037254*z3 + 0.027972027972028003*z4 - 0.021719457013574667*z5 + 0.017337461300309585*z6;
   		t_422 = 0.15238095238095234 - 0.15238095238095234*z + 0.13852813852813856*z2 - 0.12432012432012436*z3 + 0.11188811188811187*z4 - 0.1013574660633484*z5 + 0.09246646026831787*z6;
   		t_403 = 0.9142857142857144 - 1.2190476190476192*z + 1.3852813852813852*z2 - 1.4918414918414917*z3 + 1.5664335664335665*z4 - 1.6217194570135747*z5 + 1.6643962848297214*z6;
	}
	else
	{
		printf("compute_hypergeometric_functions_n_equals_4 error: z = %lf is out of bounds\n", z);
		exit(-1);
	}
}


void aniso_transport_coefficients_nonconformal::compute_transport_coefficients(precision e, precision p, precision pl, precision pt, precision b, precision beq, precision lambda, precision aT, precision aL, precision mbar, precision mass, precision mdmde)
{
	// can I cross-check with the conformal case by setting mbar = 0?
	precision lambda2 = lambda * lambda;
	precision aT2 = aT * aT;
	precision aL2 = aL * aL;

	precision prefactor_n_equals_0 = g * aT2 * aL * lambda2 / (4. * M_PI * M_PI);
	precision prefactor_n_equals_2 = g * aT2 * aL * lambda2 * lambda2 / (4. * M_PI * M_PI);
	precision prefactor_n_equals_4 = g * aT2 * aL * lambda2 * lambda2 * lambda2 / (4. * M_PI * M_PI);

	precision mbar2 = mbar * mbar;
	precision aT2_minus_aL2 = aT2 - aL2;

	// integration variables
	precision pbar, weight, exp_factor, total_weight;
	precision w, w2, w3, z;

	precision I_020 = 0;								// n = 0 moments
	precision I_001 = 0;

	precision I_240 = 0;								// n = 2 moments
	precision I_221 = 0;
	precision I_202 = 0;

	precision I_421 = 0;								// n = 4 moments (think about which ones I really need to calculate)
	precision I_402 = 0;
	precision I_441 = 0;
	precision I_422 = 0;
	precision I_403 = 0;

	for(int i = 0; i < pbar_pts; i++)					// gauss integration loop by hand (can reduce time this way)
	{
		// compute n = 0 moments
		//::::::::::::::::::::::::::::::::::::::::::::
		pbar = pbar_root_a0[i];							// pbar roots / weights for a = n = 0
		weight = pbar_weight_a0[i];

		exp_factor = exp(pbar - sqrt(pbar * pbar  +  mbar2));
		total_weight = pbar * weight * exp_factor; 		// common weight for a = n = 0

		w = sqrt(aL2  +  mbar2 / (pbar * pbar));
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

		exp_factor = exp(pbar - sqrt(pbar * pbar  +  mbar2));
		total_weight = pbar * weight * exp_factor; 		// common weight for a = n = 2

		w = sqrt(aL2  +  mbar2 / (pbar * pbar));
		w2 = w * w;
		w3 = w2 * w;
		z = aT2_minus_aL2 / w2;

		compute_hypergeometric_functions_n_equals_2(z);

		I_240 += total_weight * t_240 / w3;
		I_221 += total_weight * t_221 / w3;
		I_202 += total_weight * t_202 / w3;
		//::::::::::::::::::::::::::::::::::::::::::::



		// compute n = 4 moments (I don't have anything to test against, maybe mathematica???)
		//::::::::::::::::::::::::::::::::::::::::::::
	#if (NUMBER_OF_RESIDUAL_CURRENTS != 0)
		pbar = pbar_root_a4[i];							// pbar roots / weights for a = n = 4
		weight = pbar_weight_a4[i];

		exp_factor = exp(pbar - sqrt(pbar * pbar  +  mbar2));
		total_weight = pbar * weight * exp_factor; 		// common weight for a = n = 4

		w = sqrt(aL2  +  mbar2 / (pbar * pbar));
		w2 = w * w;
		w3 = w2 * w;
		z = aT2_minus_aL2 / w2;

		compute_hypergeometric_functions_n_equals_4(z);

		I_421 += total_weight * t_421 / w;
		I_402 += total_weight * t_402 / w;
		I_441 += total_weight * t_441 / w3;
		I_422 += total_weight * t_422 / w3;
		I_403 += total_weight * t_403 / w3;
	#endif

	}

	I_020 *= aL2 * prefactor_n_equals_0;
	I_001 *= aT2 * prefactor_n_equals_0 / 2.;

	I_240 *= aL2 * aL2 * prefactor_n_equals_2;
	I_221 *= aT2 * aL2 * prefactor_n_equals_2 / 2.;
	I_202 *= aT2 * aT2 * prefactor_n_equals_2 / 8.;

#if (NUMBER_OF_RESIDUAL_CURRENTS != 0)
	I_421 *= aT2 * aL2 * prefactor_n_equals_4 / 2.;
	I_402 *= aT2 * aT2 * prefactor_n_equals_4 / 8.;
	I_441 *= aT2 * aL2 * aL2 * prefactor_n_equals_4 / 2.;
	I_422 *= aT2 * aT2 * aL2 * prefactor_n_equals_4 / 8.;
	I_403 *= aT2 * aT2 * aT2 * prefactor_n_equals_4 / 48.;
#endif


	// note: I don't remember how quasiparticle terms ~ mdmde enter in transport coefficients
	//		 (this is based off my paper, maybe there are notes somewhere...)

	// # unstable
	precision trace_term = (2.*pt + pl + 4.*b - e) / (mass * mass);

	// # option 1,2
	//precision trace_term = (2.*pt + pl + 4.*beq - e) / (mass * mass);

	// option #3
	//precision trace_term = (3.*p + 4.*beq - e) / (mass * mass);			// equilibrium part
	//precision Pi = (pl + 2.*pt)/3. - p;									// bulk pressure


	// option #2,3 need to gneralize to 3+1d



	// pl transport coefficients


	// oh this should have beq + db
	zeta_LL = I_240  -  3. * (pl + b)  +  mdmde * (e + pl) * (I_020 + trace_term);										// option #1 (or previous)
	//zeta_LL = I_240  -  3. * (pl + b)  +  mdmde * ((e + pl) * I_020  +  (e + p) * trace_term);							// option #2
	//zeta_LL = I_240  -  3. * (pl + b)  +  mdmde * ((e + p) * I_020  +  (e + p) * trace_term);							// option #2 (simplify mdot)


	//zeta_LL = I_240  -  3. * (pl + b)  +  mdmde * ((e + pl) * (I_020 + trace_term)  +  3. * (e + p) * Pi / (mass * mass));  // option #3



	zeta_TL = I_221  -  pl  -  b  +  mdmde * (e + pt) * (I_020 + trace_term);					// previous
	//zeta_TL = I_221  -  pl  -  b  +  mdmde * ((e + pt) * I_020  +  (e + p) * trace_term);		// option #2
	//zeta_TL = I_221  -  pl  -  b  +  mdmde * ((e + p) * I_020  +  (e + p) * trace_term);		// option #2  (simplify mdot)



#ifdef WTZMU
	lambda_WuL = I_441 / I_421  +  mdmde * (I_020 + trace_term);		// previous
	//lambda_WuL = I_441 / I_421  +  mdmde * I_020;						// option #2
	lambda_WTL = 1. - lambda_WuL;
#else
	lambda_WuL = 0;
	lambda_WTL = 0;
#endif
#ifdef PIMUNU
	lambda_piL = I_422 / I_402  +  mdmde * (I_020 + trace_term);		// changed to + trace on 5/14/20
	//lambda_piL = I_422 / I_402  +  mdmde * I_020;						// option #2
#else
	lambda_piL = 0;
#endif


	// pt transport coefficients
	zeta_LT = I_221  -  pt  -  b  +  mdmde * (e + pl) * (I_001  +  trace_term);										// option #1 (or previous)
	//zeta_LT = I_221  -  pt  -  b  +  mdmde * ((e + pl) * I_001  +  (e + p) * trace_term);								// option #2
	//zeta_LT = I_221  -  pt  -  b  +  mdmde * ((e + p) * I_001  +  (e + p) * trace_term);								// option #2 (simplify mdot)

	//zeta_LT = I_221  -  pt  -  b  +  mdmde * ((e + pl) * (I_001 + trace_term)  +  3. * (e + p) * Pi / (mass * mass));	// option #3



	zeta_TT = 2. * (I_202 - pt - b)  +  mdmde * (e + pt) * (I_001  +  trace_term);				// previous
	//zeta_TT = 2. * (I_202 - pt - b)  +  mdmde * ((e + pt) * I_001  +  (e + p) * trace_term);		// option #2
	//zeta_TT = 2. * (I_202 - pt - b)  +  mdmde * ((e + p) * I_001  +  (e + p) * trace_term);		// option #2 (simplify mdot)

#ifdef WTZMU
	lambda_WTT = 2. * I_422 / I_421  +  mdmde * (I_001 + trace_term);		// previous
	//lambda_WTT = 2. * I_422 / I_421  +  mdmde * I_001;						// option #2
	lambda_WuT = lambda_WTT  -  1.;
#else
	lambda_WTT = 0;
	lambda_WuT = 0;
#endif

#ifdef PIMUNU
	lambda_piT = 1.  -  3. * I_403 / I_402  -  mdmde * (I_001 + trace_term);	// previous
	//lambda_piT = 1.  -  3. * I_403 / I_402  -  mdmde * I_001;					// previous
#else
	lambda_piT = 0;
#endif


	// WTz transport coefficients
#ifdef WTZMU
	eta_uW = 0.5 * (pl + b - I_221);
	eta_TW = 0.5 * (pt + b - I_221);
	tau_zW = pl - pt;
	lambda_WTW = 2. * I_422 / I_421  -  1.;
	lambda_WuW = 2.  -  (I_441  +  mdmde * (e + pl) * I_221) / I_421;
	delta_WW = lambda_WTW  -  0.5  +  mdmde * (e + pt) * I_221 / I_421;
#ifdef PIMUNU
	lambda_piuW = I_422 / I_402;
	lambda_piTW = lambda_piuW - 1.;
#else
	lambda_piuW = 0;
	lambda_piTW = 0;
#endif
#endif


	// piT transport coefficients
#ifdef PIMUNU
	eta_T 		= pt + b - I_202;
	tau_pipi 	= 2.  -  4. * I_403 / I_402;
	delta_pipi	= (3. * tau_pipi  + 2.) / 4.  -  mdmde * (e + pt) * I_202 / I_402;
	lambda_pipi = I_422 / I_402  -  1.  +  mdmde * (e + pl) * I_202 / I_402;
#ifdef WTZMU
	lambda_Wupi = lambda_WTW  -  1.;
	// lambda_WTpi = lambda_WuW  +  2.;
	lambda_WTpi = lambda_Wupi  +  2.;		// fixed bug on 2/25/21 (haven't tested change to 3+1d nonconformal simulation)
#else
	lambda_Wupi = 0;
	lambda_WTpi = 0;
#endif
#endif

	//test_kinetic_solution(e, pl, pt, z, aL2, prefactor);


}





