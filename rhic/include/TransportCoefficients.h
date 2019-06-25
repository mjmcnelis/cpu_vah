
#ifndef TRANSPORTCOEFFICIENTS_H_
#define TRANSPORTCOEFFICIENTS_H_

#include "Precision.h"
using namespace std;


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
    	precision Lambda;
    	precision aT;
    	precision aL;

    	precision w;
    	precision z;
    	precision t;

		precision I_240;		// make a list of all the functions I need
		precision I_221;
		precision I_202;

		transport_coefficients();
		~transport_coefficients();

		// gauss-laguerre data
		//void load_roots_and_weights();

		// have a root-solver here (it will be called in a root-finding kernel)

		void compute_transport_coefficients(precision e, precision pl, precision pt);
		// my idea last night was to use PL-PT matching (conformal EOS)
		// to propagate initial conditions instead of free-streaming
		// vahydro captures free-streaming (the problem is starting too early)
};


#endif




