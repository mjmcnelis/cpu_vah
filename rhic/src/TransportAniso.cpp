#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include "../include/Macros.h"
#include "../include/Hydrodynamics.h"
#include "../include/TransportAniso.h"
#include "../include/DynamicalVariables.h"
#include "../include/EquationOfState.h"
#include "../include/Precision.h"
#include "../include/Parameters.h"

double de_error = 1.e-7;
double dpl_error = 1.e-7;
double dpt_error = 1.e-7;


aniso_transport_coefficients::aniso_transport_coefficients()
{

}


aniso_transport_coefficients::~aniso_transport_coefficients()
{

}


void aniso_transport_coefficients::test_kinetic_solution(precision e, precision pl, precision pt, precision aL2, precision prefactor)
{
	precision e_a  = prefactor * t_200 * aL2;
	precision pl_a = prefactor * t_220 * aL2;
	precision pt_a = prefactor * t_201 / 2.;

	precision de  = fabs((e  - e_a)  / e);
	precision dpl = fabs((pl - pl_a) / pl);
	precision dpt = fabs((pt - pt_a) / pt);

	if(de > de_error || dpl > dpl_error || dpt > dpt_error)
	{
		de_error  = fmax(de,  de_error);
		dpl_error = fmax(dpl, dpl_error);
		dpt_error = fmax(dpt, dpt_error);

		printf("Transport coefficients error: |dF| = (%.6g, %.6g, %.6g)\n", de_error, dpl_error, dpt_error);
	}
}


void aniso_transport_coefficients::compute_hypergeometric_functions(precision z)
{
	precision z2 = z  * z;
	precision z3 = z2 * z;

	if(z > delta)
	{
		precision sqrtz = sqrt(z);
		precision t = atan(sqrtz) / sqrtz;

		t_200 = 1.  +  (1. + z) * t;
		t_220 = (-1.  +  (1. + z) * t) / z;
		t_201 = (1.  +  (z - 1.) * t) / z;
		t_240 = (3.  +  2. * z  -  3. * (1. + z) * t) / z2;
		t_221 = (-3.  +  (3. + z) * t) / z2;
		t_202 = (3.  +  z  +  (z - 3.) * (1. + z) * t) / (z2 * (1. + z));

	#if (NUMBER_OF_RESIDUAL_CURRENTS != 0)
		t_421 = (3.  +  z  +  (1. + z) * (z - 3.) * t) / (4. * z2);
		t_402 = (3. * (z - 1.)  +  (z * (3.*z - 2.) + 3.) * t) / (4. * z2);
		t_441 = (-15.  -  13. * z  +  3. * (1. + z) * (5. + z) * t) / (4. * z3);
		t_422 = (15.  +  z +  (z * (z - 6.) - 15.)* t) / (4. * z3);
		t_403 = ((z - 3.) * (5. + 3.*z)  +  3. * (1. + z) * (z * (z - 2.) + 5.) * t) / (4. * z3 * (1. + z));
	#endif
	}
	else if(z < -delta && z > -1.)
	{
		precision sqrtmz = sqrt(-z);
		precision t = atanh(sqrtmz) / sqrtmz;

		t_200 = 1.  +  (1. + z) * t;
		t_220 = (-1.  +  (1. + z) * t) / z;
		t_201 = (1.  +  (z - 1.) * t) / z;
		t_240 = (3.  +  2. * z  -  3. * (1. + z) * t) / z2;
		t_221 = (-3.  +  (3. + z) * t) / z2;
		t_202 = (3.  +  z  +  (z - 3.) * (1. + z) * t) / (z2 * (1. + z));

	#if (NUMBER_OF_RESIDUAL_CURRENTS != 0)
		t_421 = (3.  +  z  +  (1. + z) * (z - 3.) * t) / (4. * z2);
		t_402 = (3. * (z - 1.)  +  (z * (3.*z - 2.) + 3.) * t) / (4. * z2);
		t_441 = (-15.  -  13. * z  +  3. * (1. + z) * (5. + z) * t) / (4. * z3);
		t_422 = (15.  +  z +  (z * (z - 6.) - 15.)* t) / (4. * z3);
		t_403 = ((z - 3.) * (5. + 3.*z)  +  3. * (1. + z) * (z * (z - 2.) + 5.) * t) / (4. * z3 * (1. + z));
	#endif
	}
	else if(fabs(z) <= delta)
	{
		precision z4 = z3 * z;
		precision z5 = z4 * z;
		precision z6 = z5 * z;

		t_200 = 2. + 0.6666666666666667*z - 0.1333333333333333*z2 + 0.05714285714285716*z3 - 0.031746031746031744*z4 + 0.020202020202020193*z5 - 0.013986013986013984*z6;

		t_220 = 0.6666666666666667 - 0.1333333333333333*z + 0.05714285714285716*z2 - 0.031746031746031744*z3 + 0.020202020202020193*z4 - 0.013986013986013984*z5 + 0.010256410256410262*z6;

		t_201 = 1.3333333333333333 - 0.5333333333333333*z + 0.34285714285714286*z2 - 0.25396825396825395*z3 + 0.20202020202020202*z4 - 0.16783216783216784*z5 + 0.14358974358974358*z6;

		t_240 = 0.4 - 0.17142857142857149*z + 0.09523809523809523*z2 - 0.06060606060606058*z3 + 0.04195804195804195*z4 - 0.030769230769230785*z5 + 0.023529411764705882*z6;

		t_221 = 0.2666666666666668 - 0.22857142857142854*z + 0.19047619047619047*z2 - 0.1616161616161616*z3 + 0.13986013986013987*z4 - 0.12307692307692308*z5 + 0.10980392156862744*z6;

		t_202 = 1.0666666666666664 - 1.3714285714285712*z + 1.5238095238095237*z2 - 1.616161616161616*z3 + 1.6783216783216781*z4 - 1.7230769230769227*z5 + 1.756862745098039*z6;

	#if (NUMBER_OF_RESIDUAL_CURRENTS != 0)
   		t_421 = 0.2666666666666666 - 0.0761904761904762*z + 0.0380952380952381*z2 - 0.023088023088023088*z3 + 0.015540015540015537*z4 - 0.011188811188811189*z5 + 0.00844645550527904*z6;

   		t_402 = 1.0666666666666667 - 0.4571428571428572*z + 0.3047619047619048*z2 - 0.23088023088023088*z3 + 0.1864801864801865*z4 - 0.15664335664335666*z5 + 0.13514328808446457*z6;

	   	t_441 = 0.1142857142857145 - 0.07619047619047613*z + 0.051948051948051896*z2 - 0.037296037296037254*z3 + 0.027972027972028003*z4 - 0.021719457013574667*z5 + 0.017337461300309585*z6;

   		t_422 = 0.15238095238095234 - 0.15238095238095234*z + 0.13852813852813856*z2 - 0.12432012432012436*z3 + 0.11188811188811187*z4 - 0.1013574660633484*z5 + 0.09246646026831787*z6;

   		t_403 = 0.9142857142857144 - 1.2190476190476192*z + 1.3852813852813852*z2 - 1.4918414918414917*z3 + 1.5664335664335665*z4 - 1.6217194570135747*z5 + 1.6643962848297214*z6;
	#endif
	}
	else
	{
		printf("Error: z = %lf is out of bounds\n", z);
		exit(-1);
	}
}


void aniso_transport_coefficients::compute_transport_coefficients(precision e, precision pl, precision pt, precision conformal_eos_prefactor)
{
#ifdef CONFORMAL_EOS
	precision x   = pl / e;
	precision x2  = x   * x;			// perhaps an interpolation helps
	precision x3  = x2  * x;
	precision x4  = x3  * x;
	precision x5  = x4  * x;
	precision x6  = x5  * x;
	precision x7  = x6  * x;
	precision x8  = x7  * x;
	precision x9  = x8  * x;
	precision x10 = x9  * x;
	precision x11 = x10 * x;
	precision x12 = x11 * x;
	precision x13 = x12 * x;
	precision x14 = x13 * x;
	precision x15 = x14 * x;
	precision x16 = x15 * x;
	precision x17 = x16 * x;
	precision x18 = x17 * x;
	precision x19 = x18 * x;
	precision x20 = x19 * x;
	precision x21 = x20 * x;
	precision x22 = x21 * x;

	precision aL, aL2, z;

	if(x > 1./3.)
	{
		z = (33920.75424130814 - 72210.38155086128*x - 58150.373221605056*x2 - 123258.68865968155*x3 + 14269.18945991164*x4 + 170364.85584208343*x5 +
			169910.50665817957*x6 + 64799.56668758523*x7 + 262850.50962558796*x8 + 35449.81323782106*x9 - 248808.80620651352*x10 -
			288836.0950617432*x11 - 55525.59817904083*x12 - 249947.48251438234*x13 - 310420.8253593438*x14 - 48317.55700989516*x15 +
			522691.7302032236*x16 + 527504.8150488662*x17 + 219759.88337782127*x18 - 187603.57642353655*x19 - 199506.45061878706*x20 -
			611077.848257917*x21 + 432142.0012199023*x22)/
			(-192.35667843131404 + 43491.082537769165*x + 83354.47448892899*x2 + 45103.07343085356*x3 - 105414.36804542418*x4 + 140186.71296754244*x5 -
			531082.9994828509*x6 - 85658.91194364589*x7 + 377783.60198413196*x8 - 339045.0410056553*x9 - 95837.02795785779*x10 +
			284537.5663725089*x11 + 703062.1998023012*x12 + 223019.9316692852*x13 - 501784.5491947427*x14 - 145230.1534789184*x15 +
			55948.62853147295*x16 + 49679.34805386173*x17 - 641771.3022609851*x18 - 274804.41532698454*x19 + 726388.8998660464*x20 +
			350014.57800287893*x21 - 361748.9148710701*x22);

		if(!(z >= -0.99999999 && z <= 1.e-8))
		{
		#ifdef FLAGS
			printf("z = %lf is out of range (x = %lf). Putting z in bounds\n", z, x);
		#endif
			z = fmax(-0.99999999, fmin(z, 1.e-8));
		}
		aL  = 1. / sqrt(1. + z);
		aL2 = 1. / (1. + z);
	}
	else
	{
		aL = (2.372796737893896e-62 + 9.355496760751141e-54*x + 3.15985529218801e-46*x2 + 2.29804656071578e-39*x3 + 4.8654069671748624e-33*x4 +
			3.4686835009134695e-27*x5 + 9.052410236842743e-22*x6 + 9.132309729581051e-17*x7 + 3.705485165853083e-12*x8 + 6.240802836268058e-8*x9 +
			0.00044799689605487286*x10 + 1.4025011569370325*x11 + 1953.0793537979494*x12 + 1.229812256787706e6*x13 + 3.543561225712354e8*x14 +
			4.697330865356272e10*x15 + 2.8499566740003765e12*x16 + 7.731610782177606e13*x17 + 8.841791912264315e14*x18 + 3.673281425421166e15*x19 +
			3.042059896930142e15*x20 - 3.4368817938638095e15*x21 - 8.169907788507815e14*x22)/
			(7.281820681114894e-58 + 6.793723008169782e-50*x + 1.0073238134263982e-42*x2 + 3.8567133904345664e-36*x3 + 4.6749055427591935e-30*x4 +
			1.9992235460663164e-24*x5 + 3.2233724058457452e-19*x6 + 2.051207320369606e-14*x7 + 5.334552198382988e-10*x8 + 5.833728132253219e-6*x9 +
			0.02748383972843651*x10 + 56.954350298361284*x11 + 52824.406590310646*x12 + 2.2217655338084057e7*x13 + 4.267549397728813e9*x14 +
			3.733806109621652e11*x15 + 1.459513002063948e13*x16 + 2.4180382199020853e14*x17 + 1.4786509784350255e15*x18 + 1.8740406611426415e15*x19 -
			3.345323820802959e15*x20 - 1.2075997985771218e15*x21 + 1.136213305508547e15*x22);

   		if(!(aL >= 0.0001 && aL <= 1.0001))
		{
		#ifdef FLAGS
			printf("aL = %lf is out of range [%lf, %lf]. Putting aL in bounds\n", aL, 0.0001, 1.0001);
		#endif
			aL = fmax(0.0001, fmin(aL, 1.0001));
		}

		aL2 = aL * aL;

		z = 1. / aL2  -  1.;
	}
	// rational polynomial fit of aL as a function of pl / e
	// precision aL = (5.6098342562962155e-24 + 1.0056714201158781e-17*x + 8.574287549260127e-13*x2 + 8.639689853874967e-9*x3 + 0.000014337184308704522*x4 +
 //     0.0047402683487226555*x5 + 0.3461801244895056*x6 + 5.3061287395562*x7 + 3.7804213528647956*x8 - 55.646719325650224*x9 +
 //     71.68906037132133*x10 + 0.6485422288016947*x11 - 52.86438720903515*x12 + 32.635674688615836*x13 - 5.899614102635062*x14)/
 //   (1.2460117685059638e-20 + 3.9506205613753145e-15*x + 1.090135069930889e-10*x2 + 4.2931027828550746e-7*x3 + 0.00030704101799886117*x4 +
 //     0.04575504592470687*x5 + 1.4479634250149949*x6 + 6.077429142899631*x7 - 29.171395065126873*x8 + 13.501854646832847*x9 +
 //     65.98203155631907*x10 - 111.65365949648432*x11 + 71.83676912638525*x12 - 19.66184593458614*x13 + 1.5947903161928916*x14);


 //   	aL = max(0.001, min(aL, 20.0));

	compute_hypergeometric_functions(z);

	precision Lambda4 = 2. * e / (aL2 * conformal_eos_prefactor * t_200);
	precision Lambda2 = sqrt(Lambda4);

	precision prefactor = conformal_eos_prefactor * Lambda4 / 2.;

	// anisotropic functions
	precision I_240 = prefactor * t_240 * aL2;
	precision I_221 = prefactor * t_221 / 2.;
	precision I_202 = prefactor * t_202 / (8. * aL2);

#if (NUMBER_OF_RESIDUAL_CURRENTS != 0)
	precision I_421 = prefactor * t_421 * Lambda2 * aL2 * 10.;
	precision I_402 = prefactor * t_402 * Lambda2 * 2.5;
	precision I_441 = prefactor * t_441 * Lambda2 * aL2 * 10.;
	precision I_422 = prefactor * t_422 * Lambda2 * 2.5;
	precision I_403 = prefactor * t_403 * Lambda2 * 5. / (12. * aL2);
#endif


	// pl transport coefficients
	zeta_LL = I_240  -  3. * pl;
	zeta_TL = I_221 - pl;
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


	// pt transport coefficients
	zeta_LT = I_221  -  pt;
	zeta_TT = 2. * (I_202 - pt);
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
	// lambda_WTpi = lambda_WuW  +  2.;
	lambda_WTpi = lambda_Wupi  +  2.;		// fixed on 2/25/21 (have not tested change to 3+1d conformal simulation)
#else
	lambda_Wupi = 0;
	lambda_WTpi = 0;
#endif
#endif

	//test_kinetic_solution(e, pl, pt, z, aL2, prefactor);

#else
	// nonconformal formula
#endif
}





