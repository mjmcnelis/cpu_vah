
#ifndef TRANSPORTCOEFFICIENTS_H_
#define TRANSPORTCOEFFICIENTS_H_

#include "Precision.h"
#include "Macros.h"
#include "DynamicalVariables.h"
#include "Parameters.h"

// todo: possibly use cubic spline interpolation for the xi data (could be more accurate than a fit)

const precision delta = 0.01;	// piecewise interval where hypergeometric functions are Taylor expanded

precision eta_over_s(precision T, hydro_parameters hydro);
precision zeta_over_s(precision T, hydro_parameters hydro);

class transport_coefficients
{
	private:
		int alpha;			// # generalized powers (0 : alpha - 1)
		int points;			// # quadrature points
		double ** root;		// roots and weights for Gauss-Laguerre quadrature
    	double ** weight;

		// declare this once (in hydrodynamics) and pass it via pointer I guess
		// waste as little time opening and reading the files

	public:

		// pl coefficients
		precision zeta_LL;
		precision zeta_TL;
		precision lambda_WuL;
		precision lambda_WTL;
		precision lambda_piL;

		// pt coefficients
	#if (PT_MATCHING == 1)
		precision zeta_LT;
		precision zeta_TT;
		precision lambda_WuT;
		precision lambda_WTT;
		precision lambda_piT;
	#endif

		// WTz coefficients
	#ifdef WTZMU
		precision eta_uW;
		precision eta_TW;
		precision tau_zW;
		precision delta_WW;
		precision lambda_WuW;
		precision lambda_WTW;
		precision lambda_piuW;
		precision lambda_piTW;
	#endif

		// piT coefficients
	#ifdef PIMUNU
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

	#if (PT_MATCHING == 1 || PIMUNU_COMPONENTS != 0)
		precision t_202;
	#endif

		// for now just assume WTz and pi both included (easier to sort organize)
	#if (NUMBER_OF_RESIDUAL_CURRENTS != 0)
		precision t_441;
		precision t_421;
		precision t_402;

		precision t_422;
		precision t_403;
	#endif


		transport_coefficients();
		~transport_coefficients();

		void test_kinetic_solution(precision e, precision pl, precision pt, precision aL2, precision prefactor);

		void compute_hypergeometric_functions(precision z);

		// gauss-laguerre data
		//void load_roots_and_weights();

		// have a root-solver here (it will be called in a root-finding kernel)

		void compute_transport_coefficients(precision e, precision pl, precision pt, precision conformal_eos_prefactor);
		// my idea last night was to use PL-PT matching (conformal EOS)
		// to propagate initial conditions instead of free-streaming
		// vahydro captures free-streaming (the problem is starting too early, may need change of variables)
};


#endif




