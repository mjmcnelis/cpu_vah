
#ifndef TRANSPORTANISO_H_
#define TRANSPORTANISO_H_

#include "Precision.h"
#include "Macros.h"
#include "DynamicalVariables.h"
#include "Parameters.h"

// todo: possibly use cubic spline interpolation for the xi data (could be more accurate than a fit)
// I don't think this is really critical

const precision delta = 0.01;	// piecewise interval where hypergeometric functions are Taylor expanded

class aniso_transport_coefficients
{
	// conformal vahydro transport coefficients (nonconformal case will be a separate class)

	private:
		int alpha;			// # generalized powers (0 : alpha - 1)
		int points;			// # quadrature points
		double ** root;		// roots and weights for Gauss-Laguerre quadrature
    	double ** weight;

		// declare this once (in hydrodynamics) and pass it via pointer I guess
		// waste as little time opening and reading the files

	public:

		precision zeta_LL;			// pl coefficients
		precision zeta_TL;
		precision lambda_WuL;
		precision lambda_WTL;
		precision lambda_piL;

		precision zeta_LT;			// pt coefficients
		precision zeta_TT;
		precision lambda_WuT;
		precision lambda_WTT;
		precision lambda_piT;

	#ifdef WTZMU
		precision eta_uW;			// WTz coefficients
		precision eta_TW;
		precision tau_zW;
		precision delta_WW;
		precision lambda_WuW;
		precision lambda_WTW;
		precision lambda_piuW;
		precision lambda_piTW;
	#endif

	#ifdef PIMUNU 					// piT coefficients
		precision eta_T;
		precision delta_pipi;
		precision tau_pipi;
		precision lambda_pipi;
		precision lambda_Wupi;
		precision lambda_WTpi;
	#endif

		// hypergeometric functions needed
		precision t_200;
		precision t_220;
		precision t_201;

		precision t_240;
		precision t_221;
		precision t_202;

		// for now just assume WTz and pi both included (easier to sort organize)
	#if (NUMBER_OF_RESIDUAL_CURRENTS != 0)
		precision t_441;
		precision t_421;
		precision t_402;

		precision t_422;
		precision t_403;
	#endif


		aniso_transport_coefficients();
		~aniso_transport_coefficients();

		void test_kinetic_solution(precision e, precision pl, precision pt, precision aL2, precision prefactor);

		void compute_hypergeometric_functions(precision z);

		void compute_transport_coefficients(precision e, precision pl, precision pt, precision conformal_eos_prefactor);
};


#endif




