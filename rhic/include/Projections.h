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
		precision Xiyy;
		precision Xiyn;
		precision Xinn;

	public:
		// start with this
		precision Xitt_tt;
		precision Xitt_tx;
		precision Xitt_ty;
		precision Xitt_tn;
		precision Xitt_xx;
		precision Xitt_xy;
		precision Xitt_xn;
		precision Xitt_yy;
		precision Xitt_yn;
		precision Xitt_nn;

		// compute Xi and Xiuvab
		double_transverse_projection(precision Xitt_in, precision Xitx_in, precision Xity_in, precision Xitn_in, precision Xixx_in, precision Xixy_in, precision Xixn_in,
			precision Xiyy_in, precision Xiyn_in, precision Xinn_in);

		~double_transverse_projection();

		compute_double_transverse_projector();
};

#endif