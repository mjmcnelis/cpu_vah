
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h> // for math functions

#include "../include/TransportCoefficients.h"
#include "../include/EquationOfState.h"
#include "../include/Hypergeometric.h"


inline precision hyper(precision z)
{
	if(z > delta) 
	{
		precision sqrtz = sqrt(z);
		return atan(sqrtz) / sqrtz;
	}
	else if(z < - delta && z > -1.0)
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


void transport_coefficients::compute_anisotropic_parameters_conformal(precision e_scaled, precision x)
{	
	// x = pl / peq 
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
	aL = (2.307660683188896e-22 + 1.7179667824677117e-16*x + 7.2725449826862375e-12*x2 + 4.2846163672079405e-8*x3 + 0.00004757224421671691*x4 +
	 0.011776118846199547*x5 + 0.7235583305942909*x6 + 11.582755440134724*x7 + 44.45243622597357*x8 + 12.673594148032494*x9 -
	 33.75866652773691*x10 + 8.04299287188939*x11 + 1.462901772148128*x12 - 0.6320131889637761*x13 + 0.048528166213735346*x14)/
	(5.595674409987461e-19 + 8.059757191879689e-14*x + 1.2033043382301483e-9*x2 + 2.9819348588423508e-6*x3 + 0.0015212379997299082*x4 +
	 0.18185453852532632*x5 + 5.466199358534425*x6 + 40.1581708710626*x7 + 44.38310108782752*x8 - 55.213789667214364*x9 +
	1.5449108423263358*x10 + 11.636087951096759*x11 - 4.005934533735304*x12 + 0.4703844693488544*x13 - 0.014599143701745957*x14);

	z = 1.0 / (aL * aL)  -  1.0;	// z = xi (conformal limit)
	t = hyper(z);
	
	Lambda = pow(2.0 * e_scaled / (aL * aL * t_200(z, t)), 0.25);	
}


void transport_coefficients::compute_transport_coefficients(precision e, precision pl, precision pt)
{
#ifdef CONFORMAL_EOS
	compute_anisotropic_parameters_conformal(e / EOS_FACTOR, 3.0 * pl / e);

	precision aL2 = aL  * aL;
	precision aL3 = aL2 * aL;
	precision aL4 = aL3 * aL;
	precision aL5 = aL4 * aL;
	precision Lambda4 = Lambda * Lambda * Lambda * Lambda;

	I_240 = EOS_FACTOR * aL2 * Lambda4 * t_240(z, t) / 2.0;
	I_221 = EOS_FACTOR       * Lambda4 * t_221(z, t) / 4.0;
#else
	// nonconformal formula
#endif
}





