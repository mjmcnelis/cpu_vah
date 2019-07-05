
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h> // for math functions
#include <cmath>
#include "../include/TransportCoefficients.h"
#include "../include/DynamicalVariables.h"
#include "../include/EquationOfState.h"
#include "../include/Hypergeometric.h"


double de_error  = 1.e-7;
double dpl_error = 1.e-7;
double dpt_error = 1.e-7;


inline precision hyper(precision z)
{
	precision sqrtz = sqrt(fabs(z));

	if(z > delta)
	{
		return atan(sqrtz) / sqrtz;
	}

	return atanh(sqrtz) / sqrtz;
}


void test_kinetic_formulas(precision e, precision pl, precision pt, precision z, precision t, precision aL2, precision Lambda4)
{
	precision prefactor = EOS_FACTOR * Lambda4 / 2.0;

	precision e_a  = prefactor * t_200(z, t) * aL2;
	precision pl_a = prefactor * t_220(z, t) * aL2;
	precision pt_a = prefactor * t_201(z, t) / 2.0;

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


transport_coefficients::transport_coefficients()
{

}


transport_coefficients::~transport_coefficients()
{

}


void transport_coefficients::compute_transport_coefficients(precision e, precision pl, precision pt)
{

#ifdef CONFORMAL_EOS
	precision x = pl / e;		// x = pl / e
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

	// rational polynomial fit of aL as a function of pl / e
	precision aL = (5.6098342562962155e-24 + 1.0056714201158781e-17*x + 8.574287549260127e-13*x2 + 8.639689853874967e-9*x3 + 0.000014337184308704522*x4 +
     0.0047402683487226555*x5 + 0.3461801244895056*x6 + 5.3061287395562*x7 + 3.7804213528647956*x8 - 55.646719325650224*x9 +
     71.68906037132133*x10 + 0.6485422288016947*x11 - 52.86438720903515*x12 + 32.635674688615836*x13 - 5.899614102635062*x14)/
   (1.2460117685059638e-20 + 3.9506205613753145e-15*x + 1.090135069930889e-10*x2 + 4.2931027828550746e-7*x3 + 0.00030704101799886117*x4 +
     0.04575504592470687*x5 + 1.4479634250149949*x6 + 6.077429142899631*x7 - 29.171395065126873*x8 + 13.501854646832847*x9 +
     65.98203155631907*x10 - 111.65365949648432*x11 + 71.83676912638525*x12 - 19.66184593458614*x13 + 1.5947903161928916*x14);


   	aL = max(0.001, min(aL, 20.0));

	precision aL2 = aL  * aL;

	precision z = 1.0 / aL2  -  1.0;	// z = xi (conformal limit)
	precision t = hyper(z);


	precision Lambda4 = 2.0 * e / (aL2 * EOS_FACTOR * t_200(z, t));
	//precision Lambda  = pow(Lambda4, 0.25);

	//test_kinetic_formulas(e, pl, pt, z, t, aL2, Lambda4);

	precision prefactor = EOS_FACTOR * Lambda4 / 2.0;

	I_240 = prefactor * t_240(z, t) * aL2;
	I_221 = prefactor * t_221(z, t) / 2.0;

#if (PT_MATCHING == 1)
	I_202 = prefactor * t_202(z, t) / 8.0 / aL2;
#endif

#else
	// nonconformal formula
#endif
}





