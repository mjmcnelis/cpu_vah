
#ifndef TRANSPORTANISONONCONFORMAL_H_
#define TRANSPORTANISONONCONFORMAL_H_

#include "Precision.h"
#include "Macros.h"
#include "DynamicalVariables.h"
#include "Parameters.h"


class aniso_transport_coefficients_nonconformal
{
	// nonconformal vahydro transport coefficients 

	private:

	public:
		
		precision zeta_LL;			// pl coefficients
		precision zeta_TL;
		// precision lambda_WuL;
		// precision lambda_WTL;
		// precision lambda_piL;

		precision zeta_LT;			// pt coefficients
		precision zeta_TT;
		// precision lambda_WuT;
		// precision lambda_WTT;
		// precision lambda_piT;

	// #ifdef WTZMU
	// 	precision eta_uW;			// WTz coefficients
	// 	precision eta_TW;
	// 	precision tau_zW;
	// 	precision delta_WW;
	// 	precision lambda_WuW;
	// 	precision lambda_WTW;
	// 	precision lambda_piuW;
	// 	precision lambda_piTW;
	// #endif

	// #ifdef PIMUNU 					// piT coefficients
	// 	precision eta_T;
	// 	precision delta_pipi;
	// 	precision tau_pipi;
	// 	precision lambda_pipi;
	// 	precision lambda_Wupi;
	// 	precision lambda_WTpi;
	// #endif

		// hypergeometric functions needed
		// precision t_200;
		// precision t_220;
		// precision t_201;

		precision t_240;
		precision t_221;
		precision t_202;

	// 	// for now just assume WTz and pi both included (easier to sort organize)
	// #if (NUMBER_OF_RESIDUAL_CURRENTS != 0)
	// 	precision t_441;
	// 	precision t_421;
	// 	precision t_402;

	// 	precision t_422;
	// 	precision t_403;
	// #endif


		aniso_transport_coefficients_nonconformal();
		~aniso_transport_coefficients_nonconformal();

		void compute_hypergeometric_functions_n_equals_2(precision z);

		void compute_transport_coefficients(precision pl, precision pt, precision lambda, precision aT, precision aL, precision mbar);
		// my idea last night was to use conformal vahydro to propagate initial conditions instead of free-streaming
};


#endif




