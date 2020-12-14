#ifndef PROJECTIONS_H_
#define PROJECTIONS_H_

#include "Precision.h"

class spatial_projection
{
	private:
		precision ut;
		precision ux;
		precision uy;
		precision un;
		precision t2;
	public:
		precision Dtt;	// \Delta^{\mu\nu}
		precision Dtx;
		precision Dty;
		precision Dtn;
		precision Dxx;
		precision Dxy;
		precision Dxn;
		precision Dyy;
		precision Dyn;
		precision Dnn;

		spatial_projection(precision ut_in, precision ux_in, precision uy_in, precision un_in, precision t2_in);
		~spatial_projection();

		void spatial_project_vector(precision & At, precision & Ax, precision & Ay, precision & An);

		void test_spatial_projector();
};

class double_spatial_projection
{
	private:
		precision t2;
		precision t4;

		precision Dtt;	// \Delta^{\mu\nu}
		precision Dtx;
		precision Dty;
		precision Dtn;
		precision Dxx;
		precision Dxy;
		precision Dxn;
		precision Dyy;
		precision Dyn;
		precision Dnn;

		precision Dtt_tt;	// \Delta^{\mu\nu\alpha\beta}
		precision Dtt_tx;	// all upper indices (55 independent components)
		precision Dtt_ty;
		precision Dtt_tn;
		precision Dtt_xx;
		precision Dtt_xy;
		precision Dtt_xn;
		precision Dtt_yy;
		precision Dtt_yn;
		precision Dtt_nn;

		precision Dtx_tx;
		precision Dtx_ty;
		precision Dtx_tn;
		precision Dtx_xx;
		precision Dtx_xy;
		precision Dtx_xn;
		precision Dtx_yy;
		precision Dtx_yn;
		precision Dtx_nn;

		precision Dty_ty;
		precision Dty_tn;
		precision Dty_xx;
		precision Dty_xy;
		precision Dty_xn;
		precision Dty_yy;
		precision Dty_yn;
		precision Dty_nn;

		precision Dtn_tn;
		precision Dtn_xx;
		precision Dtn_xy;
		precision Dtn_xn;
		precision Dtn_yy;
		precision Dtn_yn;
		precision Dtn_nn;

		precision Dxx_xx;
		precision Dxx_xy;
		precision Dxx_xn;
		precision Dxx_yy;
		precision Dxx_yn;
		precision Dxx_nn;

		precision Dxy_xy;
		precision Dxy_xn;
		precision Dxy_yy;
		precision Dxy_yn;
		precision Dxy_nn;

		precision Dxn_xn;
		precision Dxn_yy;
		precision Dxn_yn;
		precision Dxn_nn;

		precision Dyy_yy;
		precision Dyy_yn;
		precision Dyy_nn;

		precision Dyn_yn;
		precision Dyn_nn;

		precision Dnn_nn;

	public:

		double_spatial_projection(spatial_projection Delta, precision t2_in, precision t4_in);
		~double_spatial_projection();

		void double_spatial_project_tensor(precision & Att, precision & Atx, precision & Aty, precision & Atn, precision & Axx, precision & Axy, precision & Axn, precision & Ayy, precision & Ayn, precision & Ann);

		void test_double_spatial_projector(precision ut, precision ux, precision uy, precision un);
};

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







