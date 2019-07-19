#ifndef PROJECTIONS_H_
#define PROJECTIONS_H_

#include "Precision.h"

using namespace std;


class double_transverse_projection
{
	private:
		precision Xitt;
		precision Xitx;
		precision Xity;
		precision Xitn;
		precision Xixx;
		precision Xixy;
		precision Xixn;

	public:

		// compute Xi and Xiuvab
		double_transverse_projection(precision Xitt_in, precision Xitx_in, precision Xity_in, precision Xitn_in, precision Xixx_in, precision Xixy_in, precision Xixn_in,
			precision Xiyy_in, precision Xiyn_in, precision Xinn_in);

		~double_transverse_projection();

		compute_double_transverse_projector();
};