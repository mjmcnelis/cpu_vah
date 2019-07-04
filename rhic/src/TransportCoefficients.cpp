
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h> // for math functions
#include <cmath>
#include "../include/TransportCoefficients.h"
#include "../include/DynamicalVariables.h"
#include "../include/EquationOfState.h"
#include "../include/Hypergeometric.h"


double de_error  = 5.e-7;
double dpl_error = 5.e-7;
double dpt_error = 5.e-7;

inline precision hyper(precision z)
{
	if(z > delta)
	{
		precision sqrtz = sqrt(z);
		return atan(sqrtz) / sqrtz;
	}
	else if(z < - delta)
	{
		precision sqrtminusz = sqrt(-z);
		return atanh(sqrtminusz) / sqrtminusz;
	}
	return 1.0;
}


transport_coefficients::transport_coefficients()
{

}


transport_coefficients::~transport_coefficients()
{

}


void transport_coefficients::compute_transport_coefficients(precision e, precision pl, precision pt)
{
	// if(std::isnan(e) || std::isnan(pl) || std::isnan(pt))
	// {
	// 	printf("Transport coefficients error: (e, pl, pt) = (%lf, %lf, %lf)\n", e, pl, pt);
	// 	exit(-1);
	// }

#ifdef CONFORMAL_EOS
	//precision x = 3.0 * pl / e;		// x = pl / peq
	precision x = pl / e;				// x = pl / e
	precision x2  = x   * x;
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

	// rational polynomial fit of aL as a function of pl / peq

	// what was the limit I imposed here?
	// aL = (2.307660683188896e-22 + 1.7179667824677117e-16*x + 7.2725449826862375e-12*x2 + 4.2846163672079405e-8*x3 + 0.00004757224421671691*x4 +
	//  0.011776118846199547*x5 + 0.7235583305942909*x6 + 11.582755440134724*x7 + 44.45243622597357*x8 + 12.673594148032494*x9 -
	//  33.75866652773691*x10 + 8.04299287188939*x11 + 1.462901772148128*x12 - 0.6320131889637761*x13 + 0.048528166213735346*x14)/
	// (5.595674409987461e-19 + 8.059757191879689e-14*x + 1.2033043382301483e-9*x2 + 2.9819348588423508e-6*x3 + 0.0015212379997299082*x4 +
	//  0.18185453852532632*x5 + 5.466199358534425*x6 + 40.1581708710626*x7 + 44.38310108782752*x8 - 55.213789667214364*x9 +
	// 1.5449108423263358*x10 + 11.636087951096759*x11 - 4.005934533735304*x12 + 0.4703844693488544*x13 - 0.014599143701745957*x14);

	precision aL = (5.6098342562962155e-24 + 1.0056714201158781e-17*x + 8.574287549260127e-13*x2 + 8.639689853874967e-9*x3 + 0.000014337184308704522*x4 +
     0.0047402683487226555*x5 + 0.3461801244895056*x6 + 5.3061287395562*x7 + 3.7804213528647956*x8 - 55.646719325650224*x9 +
     71.68906037132133*x10 + 0.6485422288016947*x11 - 52.86438720903515*x12 + 32.635674688615836*x13 - 5.899614102635062*x14)/
   (1.2460117685059638e-20 + 3.9506205613753145e-15*x + 1.090135069930889e-10*x2 + 4.2931027828550746e-7*x3 + 0.00030704101799886117*x4 +
     0.04575504592470687*x5 + 1.4479634250149949*x6 + 6.077429142899631*x7 - 29.171395065126873*x8 + 13.501854646832847*x9 +
     65.98203155631907*x10 - 111.65365949648432*x11 + 71.83676912638525*x12 - 19.66184593458614*x13 + 1.5947903161928916*x14);

   	if(!(aL >= 0.001 && aL <= 20.0))
   	{
   		//printf("aL = %lf is out of range\n", aL);
   		aL = max(0.001, min(aL, 20.0));
   	}



	// why does z go nan?

	precision aL2 = aL  * aL;
	//precision aL3 = aL2 * aL;
	//precision aL4 = aL3 * aL;
	//precision aL5 = aL4 * aL;

	precision z = 1.0 / aL2  -  1.0;	// z = xi (conformal limit)

	// if(z <= -1.0 || std::isnan(z))
	// {
	// 	printf("Transport coefficients error: (z, aL) = (%lf, %lf)\n", z, aL);
	// 	exit(-1);
	// }

	precision t = hyper(z);


	precision Lambda4 = 2.0 * e / (aL2 * EOS_FACTOR * t_200(z, t));
	//precision Lambda  = pow(Lambda4, 0.25);


	// test kinetic formulas
	precision e_a  = EOS_FACTOR * aL2 * Lambda4 * t_200(z, t) / 2.0;
	precision pl_a = EOS_FACTOR * aL2 * Lambda4 * t_220(z, t) / 2.0;
	precision pt_a = EOS_FACTOR       * Lambda4 * t_201(z, t) / 4.0;

	precision de  = fabs((e  - e_a)  / e);
	precision dpl = fabs((pl - pl_a) / pl);
	precision dpt = fabs((pt - pt_a) / pt);

	// check errors (largest error I see is 10^-7)
	if(de > de_error || dpl > dpl_error || dpt > dpt_error)
	{
		de_error  = fmax(de,  de_error);
		dpl_error = fmax(dpl, dpl_error);
		dpt_error = fmax(dpt, dpt_error);

		printf("Transport coefficients error: |dF| = (%.6g, %.6g, %.6g)\n", de_error, dpl_error, dpt_error);
	}






	I_240 = EOS_FACTOR * aL2 * Lambda4 * t_240(z, t) / 2.0;		// looks correct
	I_221 = EOS_FACTOR       * Lambda4 * t_221(z, t) / 4.0;

#if (PT_MATCHING == 1)
	I_202 = EOS_FACTOR       * Lambda4 * t_202(z, t) / 16.0 / aL2;
#endif

#else
	// nonconformal formula
#endif
}





