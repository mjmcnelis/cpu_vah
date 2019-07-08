
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
	else
	{
		return atanh(sqrtz) / sqrtz;
	}
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
	precision x   = pl / e;		
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
			printf("z = %lf is out of range\n", z);
			//exit(-1);
			z = max(-0.99999999, min(z, 1.e-8));
		}
		aL = 1.0 / sqrt(1.0 + z);
		aL2 = 1.0 / (1.0 + z);
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
			printf("aL = %lf is out of range [%lf, %lf]\n", aL, 0.0001, 1.0001);
			//exit(-1);
			aL = max(0.0001, min(aL, 1.0001));
		}

		aL2 = aL * aL;

		z = 1.0 / aL2  -  1.0;
	}
	// rational polynomial fit of aL as a function of pl / e
	// precision aL = (5.6098342562962155e-24 + 1.0056714201158781e-17*x + 8.574287549260127e-13*x2 + 8.639689853874967e-9*x3 + 0.000014337184308704522*x4 +
 //     0.0047402683487226555*x5 + 0.3461801244895056*x6 + 5.3061287395562*x7 + 3.7804213528647956*x8 - 55.646719325650224*x9 +
 //     71.68906037132133*x10 + 0.6485422288016947*x11 - 52.86438720903515*x12 + 32.635674688615836*x13 - 5.899614102635062*x14)/
 //   (1.2460117685059638e-20 + 3.9506205613753145e-15*x + 1.090135069930889e-10*x2 + 4.2931027828550746e-7*x3 + 0.00030704101799886117*x4 +
 //     0.04575504592470687*x5 + 1.4479634250149949*x6 + 6.077429142899631*x7 - 29.171395065126873*x8 + 13.501854646832847*x9 +
 //     65.98203155631907*x10 - 111.65365949648432*x11 + 71.83676912638525*x12 - 19.66184593458614*x13 + 1.5947903161928916*x14);


 //   	aL = max(0.001, min(aL, 20.0));

	precision t = hyper(z);


	//compute_aniso_functions()

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





