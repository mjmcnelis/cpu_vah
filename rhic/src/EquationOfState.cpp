
#include <math.h> // for math functions
#include <cmath>
#include "../include/Precision.h"
#include "../include/EquationOfState.h"


precision energy_density_cutoff(precision e)
{
	precision e_min = E_MIN;
	precision e_cutoff = fmax(0.0, e);

	// function approaches e_min exponentially as e -> 0 (avoids discontinuites in energy profile)
	return e_cutoff  +  e_min * exp(- e_cutoff / e_min);
	//return fmax(e_min, e);
}


// replaced powf with pow (7/2/19)
// it didn't change anything

precision equilibriumPressure(precision e)
{
#ifndef CONFORMAL_EOS
    // Equation of state from the Wuppertal-Budapest collaboration
    double e1 = (double)e;
    double e2 = e*e;
    double e3 = e2*e;
    double e4 = e3*e;
    double e5 = e4*e;
    double e6 = e5*e;
    double e7 = e6*e;
    double e8 = e7*e;
    double e9 = e8*e;
    double e10 = e9*e;
    double e11 = e10*e;
    double e12 = e11*e;

	double a0 = -0.25181736420168666;
	double a1 = 9737.845799644809;
	double a2 = 1.077580993288114e6;
	double a3 = 3.1729694865420084e6;
	double a4 = 1.6357487344679043e6;
	double a5 = 334334.4309240126;
	double a6 = 41913.439282708554;
	double a7 = 6340.448389300905;
	double a8 = 141.5073484468774;
	double a9 = 0.7158279081255019;
	double a10 = 0.0009417586777847889;
	double a11 = 3.1188455176941583e-7;
	double a12 = 1.9531729608963267e-11;
	precision a = (precision)fma(a12,e12,fma(a11,e11,fma(a10,e10,fma(a9,e9,fma(a8,e8,fma(a7,e7,fma(a6,e6,fma(a5,e5,fma(a4,e4,fma(a3,e3,fma(a2,e2,fma(a1,e1,a0))))))))))));

	double b0 = 45829.44617893836;
	double b1 = 4.0574329080826794e6;
	double b2 = 2.0931169138134286e7;
	double b3 = 1.3512402226067686e7;
	double b4 = 1.7851642641834426e6;
	double b5 = 278581.2989342773;
	double b6 = 26452.34905933697;
	double b7 = 499.04919730607065;
	double b8 = 2.3405487982094204;
	double b9 = 0.002962497695527404;
	double b10 = 9.601103399348206e-7;
	double b11 = 5.928138360995685e-11;
	double b12 = 3.2581066229887368e-18;
	precision b = (precision)fma(b12,e12,fma(b11,e11,fma(b10,e10,fma(b9,e9,fma(b8,e8,fma(b7,e7,fma(b6,e6,fma(b5,e5,fma(b4,e4,fma(b3,e3,fma(b2,e2,fma(b1,e1,b0))))))))))));

    return a/b;

#else
    return e / 3.0;
#endif
}

precision speedOfSoundSquared(precision e) {
#ifndef CONFORMAL_EOS
	// Speed of sound from the Wuppertal-Budapest collaboration
	double e1 = (double) e;
	double e2 = e * e1;
	double e3 = e2 * e1;
	double e4 = e3 * e1;
	double e5 = e4 * e1;
	double e6 = e5 * e1;
	double e7 = e6 * e1;
	double e8 = e7 * e1;
	double e9 = e8 * e1;
	double e10 = e9 * e1;
	double e11 = e10 * e1;
	double e12 = e11 * e1;
	double e13 = e12 * e1;
	return (5.191934309650155e-32 + 4.123605749683891e-23 * e
			+ 3.1955868410879504e-16 * e2 + 1.4170364808063119e-10 * e3
			+ 6.087136671592452e-6 * e4 + 0.02969737949090831 * e5
			+ 15.382615282179595 * e6 + 460.6487249985994 * e7
			+ 1612.4245252438795 * e8 + 275.0492627924299 * e9
			+ 58.60283714484669 * e10 + 6.504847576502024 * e11
			+ 0.03009027913262399 * e12 + 8.189430244031285e-6 * e13)
			/ (1.4637868900982493e-30 + 6.716598285341542e-22 * e
					+ 3.5477700458515908e-15 * e2 + 1.1225580509306008e-9 * e3
					+ 0.00003551782901018317 * e4 + 0.13653226327408863 * e5
					+ 60.85769171450653 * e6 + 1800.5461219450308 * e7
					+ 15190.225535036281 * e8 + 590.2572000057821 * e9
					+ 293.99144775704605 * e10 + 21.461303090563028 * e11
					+ 0.09301685073435291 * e12 + 0.000024810902623582917 * e13);
#else
	return 1.0/3.0;
#endif
}

precision effectiveTemperature(precision e) {
#ifndef CONFORMAL_EOS
	// Effective temperature from the Wuppertal-Budapest collaboration
	double e1 = (double) e;
	double e2 = e * e1;
	double e3 = e2 * e1;
	double e4 = e3 * e1;
	double e5 = e4 * e1;
	double e6 = e5 * e1;
	double e7 = e6 * e1;
	double e8 = e7 * e1;
	double e9 = e8 * e1;
	double e10 = e9 * e1;
	double e11 = e10 * e1;
	return (1.510073201405604e-29 + 8.014062800678687e-18 * e
			+ 2.4954778310451065e-10 * e2 + 0.000063810382643387 * e3
			+ 0.4873490574161924 * e4 + 207.48582344326206 * e5
			+ 6686.07424325115 * e6 + 14109.766109389702 * e7
			+ 1471.6180520527757 * e8 + 14.055788949565482 * e9
			+ 0.015421252394182246 * e10 + 1.5780479034557783e-6 * e11)
			/ (7.558667139355393e-28 + 1.3686372302041508e-16 * e
					+ 2.998130743142826e-9 * e2 + 0.0005036835870305458 * e3
					+ 2.316902328874072 * e4 + 578.0778724946719 * e5
					+ 11179.193315394154 * e6 + 17965.67607192861 * e7
					+ 1051.0730543534657 * e8 + 5.916312075925817 * e9
					+ 0.003778342768228011 * e10 + 1.8472801679382593e-7 * e11);
#else
	return pow(e / EOS_FACTOR, 0.25);
#endif
}

precision equilibriumEnergyDensity(precision T) {
#ifndef CONFORMAL_EOS
	// Effective temperature from the Wuppertal-Budapest collaboration
	double T1 = (double) T;
	double T2 = T1 * T1;
	double T3 = T2 * T1;
	double T4 = T3 * T1;
	double T5 = T4 * T1;
	double T6 = T5 * T1;
	double T7 = T6 * T1;
	double T8 = T7 * T1;
	double T9 = T8 * T1;
	double T10 = T9 * T1;
	double T11 = T10 * T1;
	double T12 = T11 * T1;
	double T13 = T12 * T1;
	double T14 = T13 * T1;
	double T15 = T14 * T1;
	double T16 = T15 * T1;
	double T17 = T16 * T1;
	double T18 = T17 * T1;
	double T19 = T18 * T1;
	double T20 = T19 * T1;
	double T21 = T20 * T1;
	double T22 = T21 * T1;
	double T23 = T22 * T1;
	return (-0.011958188410851651 + 119.89423098138208 * T
			- 3156.9475699248055 * T2 + 32732.86844374939 * T3
			- 187899.8994764422 * T4 + 712537.3610845465 * T5
			- 1.557049803609345e6 * T6 + 1.4852519861308339e6 * T7
			+ 532132.6079941876 * T8 - 1.963099445042592e6 * T9
			- 4484.44579242679 * T10 + 1.7984228830058286e6 * T11
			+ 119345.25619517374 * T12 - 1.3499773937058165e6 * T13
			- 207838.4995663606 * T14 + 654970.2138652403 * T15
			- 78643.00334616247 * T16 + 40274.00078068926 * T17
			+ 422619.58977657766 * T18 - 409688.07836393174 * T19
			- 62005.75915066359 * T20 + 46788.14270090656 * T21
			+ 40784.330477857235 * T22 - 12589.47744840392 * T23)
			/ (31630.074365558292 - 127100.88940643385 * T
					+ 173528.1225422275 * T2 - 39403.297956865215 * T3
					- 85582.57873541754 * T4 + 9320.560804233442 * T5
					+ 50882.74198960172 * T6 + 20335.926473421183 * T7
					- 14897.725710713818 * T8 - 23836.484117457 * T9
					- 13726.013896090335 * T10 + 4517.908673107615 * T11
					+ 18056.19917986404 * T12 + 14954.82860467155 * T13
					+ 2569.623976952738 * T14 - 9304.046211514986 * T15
					- 15606.429173842751 * T16 + 8383.710735812094 * T17
					+ 1591.3177623932843 * T18 - 678.748230997762 * T19
					- 33.58687934953277 * T20 + 3.2520554133126285 * T21
					- 0.19647288043440464 * T22 + 0.005443394551264717 * T23);
#else
	return EOS_FACTOR * pow(T, 4.0);
#endif
}


