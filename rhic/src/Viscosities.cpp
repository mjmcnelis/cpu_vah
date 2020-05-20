#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include "../include/Hydrodynamics.h"
#include "../include/Precision.h"
#include "../include/Parameters.h"

inline precision sign(precision x)
{
	if(x > 0.)
	{
		return 1.;
	}
	else if(x < 0.)
	{
		return -1.;
	}

	return 0.;
}

inline precision Theta(precision x)
{
	if(x > 0.)
	{
		return 1.;
	}

	return 0.;
}


precision eta_over_s(precision T, hydro_parameters hydro)
{
	precision etas_min = hydro.etas_min;

	if(hydro.temperature_etas == 1)						// temperature dependent piecewise parameterization
	{
		precision aL = hydro.etas_aL * hbarc;			// left slope [fm]
		precision aH = hydro.etas_aH * hbarc;			// right slope [fm]
		precision Tk = hydro.etas_Tk_GeV / hbarc;		// kink temperature fm^-1]
		precision etask = hydro.etas_etask;				// kink etas value

		return fmax(etas_min, etask  +  (T - Tk) * (aL * Theta(Tk - T)  +  aH * Theta(T - Tk)));
	}

	return fmax(etas_min, hydro.constant_etas);
}


precision zeta_over_s(precision T, hydro_parameters hydro)
{
	precision norm = hydro.zetas_normalization_factor;				// normalization factor
	precision Tpeak = hydro.zetas_peak_temperature_GeV / hbarc;		// peak temperature [fm^-1]
	precision width = hydro.zetas_width_GeV / hbarc;				// width [fm^-1]
	precision skew = hydro.zetas_skew;								// skew

	precision Lambda = width * (1.  +  skew * sign(T - Tpeak));

	return norm * Lambda * Lambda / (Lambda * Lambda  +  (T - Tpeak) * (T - Tpeak));
}


// these viscosities were used in 2018 vahydro paper

/*
precision eta_over_s(precision T, hydro_parameters hydro)
{
	if(hydro.temperature_etas == 1)
	{
		precision etas_min = 0.08;
		precision etas_slope = 0.167728;
		precision Tc = 0.154 / hbarc;

		if(T > Tc)
		{
			return etas_min + etas_slope * (T - Tc);
		}

		return etas_min;

	}

	return hydro.constant_etas;
}

precision zeta_over_s(precision T, hydro_parameters hydro)
{
	precision norm = 1.25;
	precision Tc = 0.154 / hbarc;

	precision a0 = -13.45;
	precision a1 = 27.55;
	precision a2 = -13.77;

	precision lambda1 = 0.9;
	precision lambda2 = 0.25;
	precision lambda3 = 0.9;
	precision lambda4 = 0.22;

	precision sigma1 = 0.025;
	precision sigma2 = 0.13;
	precision sigma3 = 0.0025;
	precision sigma4 = 0.022;

	precision x = T / Tc;

	if(x > 1.05)
	{
		return norm * (lambda1*exp(-(x-1.0)/sigma1) + lambda2*exp(-(x-1.0)/sigma2) + 0.001);
	}
	else if(x < 0.995)
	{
		return norm * (lambda3*exp((x-1.0)/sigma3)+ lambda4*exp((x-1.0)/sigma4) + 0.03);
	}

	return norm * (a0 + a1*x + a2*x*x);
}

*/












