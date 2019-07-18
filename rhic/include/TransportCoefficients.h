
#ifndef TRANSPORTCOEFFICIENTS_H_
#define TRANSPORTCOEFFICIENTS_H_

#include "Precision.h"
#include "DynamicalVariables.h"

using namespace std;

const precision delta = 0.01;

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
	#ifdef WTZMU
		precision lambda_WuL;
		precision lambda_WTL;
	#endif

		// pt coefficients

	#if (PT_MATCHING == 1)
		precision zeta_LT;
		precision zeta_TT;
	#ifdef WTZMU
		precision lambda_WuT;
		precision lambda_WTT;
	#endif
	#endif
	
		precision t_200;
		precision t_240;
		precision t_221;

	#ifdef WTZMU
		precision t_441;
		precision t_421;
	#endif

	#if (PT_MATCHING == 1)
		precision t_202;
	#ifdef WTZMU
		precision t_422;
	#endif
	#endif
		
		

		transport_coefficients();
		~transport_coefficients();

		void compute_hypergeometric_functions(precision z);

		// gauss-laguerre data
		//void load_roots_and_weights();

		// have a root-solver here (it will be called in a root-finding kernel)

		void compute_transport_coefficients(precision e, precision pl, precision pt);
		// my idea last night was to use PL-PT matching (conformal EOS)
		// to propagate initial conditions instead of free-streaming
		// vahydro captures free-streaming (the problem is starting too early)
};


#endif




