#ifndef PROJECTIONS_H_
#define PROJECTIONS_H_

#include "Precision.h"

using namespace std;

class transverse_projection
{
	private:
		precision ut;
		precision ux;
		precision uy;
		precision un;
		
		precision zt;
		precision zn;

		precision t2;
	public:
		precision Xitt;	// Xi^{\mu\nu}
		precision Xitx;
		precision Xity;
		precision Xitn;
		precision Xixx;
		precision Xixy;
		precision Xixn;
		precision Xiyy;
		precision Xiyn;
		precision Xinn;

		transverse_projection(precision ut_in, precision ux_in, precision uy_in, precision un_in, precision zt, precision zn, precision t2_in);

		~transverse_projection();

		void transverse_project_vector(precision & At, precision & Ax, precision & Ay, precision & An);

		void test_transverse_projector();
};


class double_transverse_projection
{
	private:
		precision t2;
		precision t4;

		precision Xitt;	// Xi^{\mu\nu}
		precision Xitx;
		precision Xity;
		precision Xitn;
		precision Xixx;
		precision Xixy;
		precision Xixn;
		precision Xiyy;
		precision Xiyn;
		precision Xinn;

		precision Xitt_tt;	// Xi^{\mu\nu\alpha\beta} 
		precision Xitt_tx;	// all upper indices (55 independent components)
		precision Xitt_ty;	
		precision Xitt_tn;
		precision Xitt_xx;
		precision Xitt_xy;
		precision Xitt_xn;
		precision Xitt_yy;
		precision Xitt_yn;
		precision Xitt_nn;

		precision Xitx_tx;
		precision Xitx_ty;
		precision Xitx_tn;
		precision Xitx_xx;
		precision Xitx_xy;
		precision Xitx_xn;
		precision Xitx_yy;
		precision Xitx_yn;
		precision Xitx_nn;

		precision Xity_ty;
		precision Xity_tn;
		precision Xity_xx;
		precision Xity_xy;
		precision Xity_xn;
		precision Xity_yy;
		precision Xity_yn;
		precision Xity_nn;

		precision Xitn_tn;
		precision Xitn_xx;
		precision Xitn_xy;
		precision Xitn_xn;
		precision Xitn_yy;
		precision Xitn_yn;
		precision Xitn_nn;

		precision Xixx_xx;
		precision Xixx_xy;
		precision Xixx_xn;
		precision Xixx_yy;
		precision Xixx_yn;
		precision Xixx_nn;

		precision Xixy_xy;
		precision Xixy_xn;
		precision Xixy_yy;
		precision Xixy_yn;
		precision Xixy_nn;

		precision Xixn_xn;
		precision Xixn_yy;
		precision Xixn_yn;
		precision Xixn_nn;

		precision Xiyy_yy;
		precision Xiyy_yn;
		precision Xiyy_nn;

		precision Xiyn_yn;
		precision Xiyn_nn;

		precision Xinn_nn;

	public:
	
		double_transverse_projection(transverse_projection Xi, precision t2_in, precision t4_in);

		~double_transverse_projection();

		void double_transverse_project_tensor(precision & Att, precision & Atx, precision & Aty, precision & Atn, precision & Axx, precision & Axy, precision & Axn, precision & Ayy, precision & Ayn, precision & Ann);

		void test_double_transverse_projector(precision ut, precision ux, precision uy, precision un, precision zt, precision zn);
};

#endif







