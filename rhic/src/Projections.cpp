#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <random>
#include <chrono>
#include "../include/Precision.h"
#include "../include/Projections.h"
#include "../include/Macros.h"
using namespace std;

inline double canonical(default_random_engine & generator)
{
  // random number between [0,1)
  return generate_canonical<double, numeric_limits<double>::digits>(generator);
}


//////////////////////////////////////////////////////////////////////


spatial_projection::spatial_projection(precision ut_in, precision ux_in, precision uy_in, precision un_in, precision t2_in)
{
	ut = ut_in;
	ux = ux_in;
	uy = uy_in;
	un = un_in;
	t2 = t2_in;

	Dtt = 1.  -  ut * ut;
	Dtx = - ut * ux;
	Dty = - ut * uy;
	Dtn = - ut * un;
	Dxx = - 1.  -  ux * ux;
	Dxy = - ux * uy;
	Dxn = - ux * un;
	Dyy = - 1.  -  uy * uy;
	Dyn = - uy * un;
	Dnn = - 1. / t2  -  un * un;
}


spatial_projection::~spatial_projection()
{

}


void spatial_projection::spatial_project_vector(precision & At, precision & Ax, precision & Ay, precision & An)
{
	precision At_pro = Dtt * At  -  Dtx * Ax  -  Dty * Ay  -  t2 * Dtn * An;
	precision Ax_pro = Dtx * At  -  Dxx * Ax  -  Dxy * Ay  -  t2 * Dxn * An;
	precision Ay_pro = Dty * At  -  Dxy * Ax  -  Dyy * Ay  -  t2 * Dyn * An;
	precision An_pro = Dtn * At  -  Dxn * Ax  -  Dyn * Ay  -  t2 * Dnn * An;

	At = At_pro;
	Ax = Ax_pro;
	Ay = Ay_pro;
	An = An_pro;
}

void spatial_projection::test_spatial_projector()
{
	// test spatial projector by projecting a random vector
	unsigned long int seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);

	precision At = canonical(generator);
	precision Ax = canonical(generator);
	precision Ay = canonical(generator);
	precision An = canonical(generator);

	spatial_project_vector(At, Ax, Ay, An);

	precision A_mag = sqrt(fabs(At * At  -  Ax * Ax  -  Ay * Ay  -  t2 * An * An));

	precision Au = fabs(At * ut  -  Ax * ux  -  Ay * uy  -  t2 * An * un) / (ut * A_mag);

#ifdef FLAGS
    if(Au > 1.e-12)  printf("test_spatial_projector flag: A is not orthogonal to u (%.6g)\n", Au);
#endif
}


//////////////////////////////////////////////////////////////////////



double_spatial_projection::double_spatial_projection(spatial_projection Delta, precision t2_in, precision t4_in)
{
	t2 = t2_in;
	t4 = t4_in;

	Dtt = Delta.Dtt;
	Dtx = Delta.Dtx;
	Dty = Delta.Dty;
	Dtn = Delta.Dtn;
	Dxx = Delta.Dxx;
	Dxy = Delta.Dxy;
	Dxn = Delta.Dxn;
	Dyy = Delta.Dyy;
	Dyn = Delta.Dyn;
	Dnn = Delta.Dnn;

	Dtt_tt =  2./3. * Dtt * Dtt;
	Dtt_tx =  2./3. * Dtt * Dtx;
	Dtt_ty =  2./3. * Dtt * Dty;
	Dtt_tn =  2./3. * Dtt * Dtn;
	Dtt_xx =  Dtx * Dtx  -  Dtt * Dxx / 3.;
	Dtt_xy =  Dtx * Dty  -  Dtt * Dxy / 3.;
	Dtt_xn =  Dtx * Dtn  -  Dtt * Dxn / 3.;
	Dtt_yy =  Dty * Dty  -  Dtt * Dyy / 3.;
	Dtt_yn =  Dty * Dtn  -  Dtt * Dyn / 3.;
	Dtt_nn =  Dtn * Dtn  -  Dtt * Dnn / 3.;

	Dtx_tx = (Dtx * Dtx  +  3. * Dtt * Dxx) / 6.;
	Dtx_ty = (Dtx * Dty  +  3. * Dtt * Dxy) / 6.;
	Dtx_tn = (Dtx * Dtn  +  3. * Dtt * Dxn) / 6.;
	Dtx_xx = 2./3. * Dtx * Dxx;
	Dtx_xy = (Dtx * Dxy  +  3. * Dty * Dxx) / 6.;
	Dtx_xn = (Dtx * Dxn  +  3. * Dtn * Dxx) / 6.;
	Dtx_yy = Dty * Dxy  -  Dtx * Dyy / 3.;
	Dtx_yn = (3. * Dty * Dxn  +  3. * Dtn * Dxy  -  2. * Dtx * Dyn) / 6.;
	Dtx_nn = Dtn * Dxn  -  Dtx * Dnn / 3.;

	Dty_ty = (Dty * Dty  +  3. * Dtt * Dyy) / 6.;
	Dty_tn = (Dty * Dtn  +  3. * Dtt * Dyn) / 6.;
	Dty_xx = Dtx * Dxy  -  Dty * Dxx / 3.;
	Dty_xy = (Dty * Dxy  +  3. * Dtx * Dyy) / 6.;
	Dty_xn = (3. * Dtx * Dyn  +  3. * Dtn * Dxy  -  2. * Dty * Dxn) / 6.;
	Dty_yy = 2./3. * Dty * Dyy;
	Dty_yn = (Dty * Dyn  +  3. * Dtn * Dyy) / 6.;
	Dty_nn = Dtn * Dyn  -  Dnn * Dty / 3.;

	Dtn_tn = (Dtn * Dtn  +  3. * Dnn * Dtt) / 6.;
	Dtn_xx = Dtx * Dxn  -  Dtn * Dxx / 3.;
	Dtn_xy = (3. * Dty * Dxn  +  3. * Dtx * Dyn  -  2. * Dtn * Dxy) / 6.;
	Dtn_xn = (Dtn * Dxn  +  3. * Dnn * Dtx) / 6.;
	Dtn_yy = Dty * Dyn  -  Dtn * Dyy / 3.;
	Dtn_yn = (Dtn * Dyn  +  3. * Dnn * Dty) / 6.;
	Dtn_nn = 2./3. * Dnn * Dtn;

	Dxx_xx = 2./3. * Dxx * Dxx;
	Dxx_xy = 2./3. * Dxx * Dxy;
	Dxx_xn = 2./3. * Dxx * Dxn;
	Dxx_yy = Dxy * Dxy  -  Dxx * Dyy / 3.;
	Dxx_yn = Dxy * Dxn  -  Dxx * Dyn / 3.;
	Dxx_nn = Dxn * Dxn  -  Dxx * Dnn / 3.;

	Dxy_xy = (Dxy * Dxy  +  3. * Dxx * Dyy) / 6.;
	Dxy_xn = (Dxn * Dxy  +  3. * Dxx * Dyn) / 6.;
	Dxy_yy = 2./3. * Dyy * Dxy;
	Dxy_yn = (Dxy * Dyn  +  3. * Dxn * Dyy) / 6.;
	Dxy_nn = Dxn * Dyn  -  Dxy * Dnn / 3.;

	Dxn_xn = (Dxn * Dxn  +  3. * Dnn * Dxx) / 6.;
	Dxn_yy = Dxy * Dyn  -  Dyy * Dxn / 3.;
	Dxn_yn = (Dxn * Dyn  +  3. * Dnn * Dxy) / 6.;
	Dxn_nn = 2./3. * Dxn * Dnn;

	Dyy_yy = 2./3. * Dyy * Dyy;
	Dyy_yn = 2./3. * Dyy * Dyn;
	Dyy_nn = Dyn * Dyn  -  Dyy * Dnn / 3.;

	Dyn_yn = (Dyn * Dyn  +  3. * Dnn * Dyy) / 6.;
	Dyn_nn = 2./3. * Dyn * Dnn;

	Dnn_nn = 2./3. * Dnn * Dnn;
}


double_spatial_projection::~double_spatial_projection()
{

}


void double_spatial_projection::double_spatial_project_tensor(precision & Att, precision & Atx, precision & Aty, precision & Atn, precision & Axx, precision & Axy, precision & Axn, precision & Ayy, precision & Ayn, precision & Ann)
{
	// A_pro^{\mu\nu} = \Delta^{\mu\nu\alpha\beta} . A_{\alpha\beta}
	precision Att_pro = Dtt_tt * Att  +  Dtt_xx * Axx  +  Dtt_yy * Ayy  +  t4 * Dtt_nn * Ann  -  2. * (Dtt_tx * Atx  +  Dtt_ty * Aty  -  Dtt_xy * Axy  +  t2 * (Dtt_tn * Atn   -  Dtt_xn * Axn  -  Dtt_yn * Ayn));
	precision Atx_pro = Dtt_tx * Att  +  Dtx_xx * Axx  +  Dtx_yy * Ayy  +  t4 * Dtx_nn * Ann  -  2. * (Dtx_tx * Atx  +  Dtx_ty * Aty  -  Dtx_xy * Axy  +  t2 * (Dtx_tn * Atn   -  Dtx_xn * Axn  -  Dtx_yn * Ayn));
	precision Aty_pro = Dtt_ty * Att  +  Dty_xx * Axx  +  Dty_yy * Ayy  +  t4 * Dty_nn * Ann  -  2. * (Dtx_ty * Atx  +  Dty_ty * Aty  -  Dty_xy * Axy  +  t2 * (Dty_tn * Atn   -  Dty_xn * Axn  -  Dty_yn * Ayn));
	precision Atn_pro = Dtt_tn * Att  +  Dtn_xx * Axx  +  Dtn_yy * Ayy  +  t4 * Dtn_nn * Ann  -  2. * (Dtx_tn * Atx  +  Dty_tn * Aty  -  Dtn_xy * Axy  +  t2 * (Dtn_tn * Atn   -  Dtn_xn * Axn  -  Dtn_yn * Ayn));
	precision Axx_pro = Dtt_xx * Att  +  Dxx_xx * Axx  +  Dxx_yy * Ayy  +  t4 * Dxx_nn * Ann  -  2. * (Dtx_xx * Atx  +  Dty_xx * Aty  -  Dxx_xy * Axy  +  t2 * (Dtn_xx * Atn   -  Dxx_xn * Axn  -  Dxx_yn * Ayn));
	precision Axy_pro = Dtt_xy * Att  +  Dxx_xy * Axx  +  Dxy_yy * Ayy  +  t4 * Dxy_nn * Ann  -  2. * (Dtx_xy * Atx  +  Dty_xy * Aty  -  Dxy_xy * Axy  +  t2 * (Dtn_xy * Atn   -  Dxy_xn * Axn  -  Dxy_yn * Ayn));
	precision Axn_pro = Dtt_xn * Att  +  Dxx_xn * Axx  +  Dxn_yy * Ayy  +  t4 * Dxn_nn * Ann  -  2. * (Dtx_xn * Atx  +  Dty_xn * Aty  -  Dxy_xn * Axy  +  t2 * (Dtn_xn * Atn   -  Dxn_xn * Axn  -  Dxn_yn * Ayn));
	precision Ayy_pro = Dtt_yy * Att  +  Dxx_yy * Axx  +  Dyy_yy * Ayy  +  t4 * Dyy_nn * Ann  -  2. * (Dtx_yy * Atx  +  Dty_yy * Aty  -  Dxy_yy * Axy  +  t2 * (Dtn_yy * Atn   -  Dxn_yy * Axn  -  Dyy_yn * Ayn));
	precision Ayn_pro = Dtt_yn * Att  +  Dxx_yn * Axx  +  Dyy_yn * Ayy  +  t4 * Dyn_nn * Ann  -  2. * (Dtx_yn * Atx  +  Dty_yn * Aty  -  Dxy_yn * Axy  +  t2 * (Dtn_yn * Atn   -  Dxn_yn * Axn  -  Dyn_yn * Ayn));
	precision Ann_pro = Dtt_nn * Att  +  Dxx_nn * Axx  +  Dyy_nn * Ayy  +  t4 * Dnn_nn * Ann  -  2. * (Dtx_nn * Atx  +  Dty_nn * Aty  -  Dxy_nn * Axy  +  t2 * (Dtn_nn * Atn   -  Dxn_nn * Axn  -  Dyn_nn * Ayn));

	Att = Att_pro;
	Atx = Atx_pro;
	Aty = Aty_pro;
	Atn = Atn_pro;
	Axx = Axx_pro;
	Axy = Axy_pro;
	Axn = Axn_pro;
	Ayy = Ayy_pro;
	Ayn = Ayn_pro;
	Ann = Ann_pro;
}


void double_spatial_projection::test_double_spatial_projector(precision ut, precision ux, precision uy, precision un)
{
	// test double projector by projecting a random symmetric matrix
	unsigned long int seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);

	precision Att = canonical(generator);
	precision Atx = canonical(generator);
	precision Aty = canonical(generator);
	precision Atn = canonical(generator);
	precision Axx = canonical(generator);
	precision Axy = canonical(generator);
	precision Axn = canonical(generator);
	precision Ayy = canonical(generator);
	precision Ayn = canonical(generator);
	precision Ann = canonical(generator);

	double_spatial_project_tensor(Att, Atx, Aty, Atn, Axx, Axy, Axn, Ayy, Ayn, Ann);

	precision A_mag = sqrt(fabs(Att * Att  +  Axx * Axx  +  Ayy * Ayy  +  t4 * Ann * Ann  -  2. * (Atx * Atx  +  Aty * Aty  -  Axy * Axy  +  t2 * (Atn * Atn  -  Axn * Axn  -  Ayn * Ayn))));

	precision Au0 = fabs(Att * ut  -  Atx * ux  -  Aty * uy  -  t2 * Atn * un);
	precision Au1 = fabs(Atx * ut  -  Axx * ux  -  Axy * uy  -  t2 * Axn * un);
	precision Au2 = fabs(Aty * ut  -  Axy * ux  -  Ayy * uy  -  t2 * Ayn * un);
	precision Au3 = fabs(Atn * ut  -  Axn * ux  -  Ayn * uy  -  t2 * Ann * un);
	precision Au = fmax(Au0, fmax(Au1, fmax(Au2, Au3))) / (ut * A_mag);

	precision trA = fabs(Att  -  Axx  -  Ayy  -  t2 * Ann) / A_mag;

#ifdef FLAGS
    if(Au > 1.e-12)  printf("test_double_spatial_projector flag: A is not orthogonal to u (%.6g)\n", Au);
    if(trA > 1.e-12) printf("test_double_spatial_projector flag: A is not traceless (%.6g)\n", trA);
 #endif
}


//////////////////////////////////////////////////////////////////////


transverse_projection::transverse_projection(precision ut_in, precision ux_in, precision uy_in, precision un_in, precision zt_in, precision zn_in, precision t2_in)
{
	ut = ut_in;		ux = ux_in;		uy = uy_in;		un = un_in;
	zt = zt_in;		zn = zn_in;
	t2 = t2_in;

	Xitt = 1.  -  ut * ut  +  zt * zt;
	Xitx = - ut * ux;
	Xity = - ut * uy;
	Xitn = - ut * un  +  zt * zn;
	Xixx = - 1.  -  ux * ux;
	Xixy = - ux * uy;
	Xixn = - ux * un;
	Xiyy = - 1.  -  uy * uy;
	Xiyn = - uy * un;
	Xinn = - 1. / t2  -  un * un  +  zn * zn;
}


transverse_projection::~transverse_projection()
{

}


void transverse_projection::transverse_project_vector(precision & At, precision & Ax, precision & Ay, precision & An)
{
	precision At_pro = Xitt * At  -  Xitx * Ax  -  Xity * Ay  -  t2 * Xitn * An;
	precision Ax_pro = Xitx * At  -  Xixx * Ax  -  Xixy * Ay  -  t2 * Xixn * An;
	precision Ay_pro = Xity * At  -  Xixy * Ax  -  Xiyy * Ay  -  t2 * Xiyn * An;
	precision An_pro = Xitn * At  -  Xixn * Ax  -  Xiyn * Ay  -  t2 * Xinn * An;

	At = At_pro;
	Ax = Ax_pro;
	Ay = Ay_pro;
	An = An_pro;
}


void transverse_projection::test_transverse_projector()
{
	// test transverse projector by projecting a random vector
	unsigned long int seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);

	precision At = canonical(generator);
	precision Ax = canonical(generator);
	precision Ay = canonical(generator);
	precision An = canonical(generator);

	transverse_project_vector(At, Ax, Ay, An);

	precision Au = fabs(At * ut  -  Ax * ux  -  Ay * uy  -  t2 * An * un);
	precision Az = fabs(At * zt  -  t2 * An * zn);

#ifdef FLAGS
    if(Au > 1.e-14)  printf("test_transverse_projector flag: A is not orthogonal to u (%.6g)\n", Au);
    if(Az > 1.e-14)  printf("test_transverse_projector flag: A is not orthogonal to z (%.6g)\n", Az);
#endif
}


//////////////////////////////////////////////////////////////////////

double_transverse_projection::double_transverse_projection(transverse_projection Xi, precision t2_in, precision t4_in)
{
	t2 = t2_in;
	t4 = t4_in;

	Xitt = Xi.Xitt;
	Xitx = Xi.Xitx;
	Xity = Xi.Xity;
	Xitn = Xi.Xitn;
	Xixx = Xi.Xixx;
	Xixy = Xi.Xixy;
	Xixn = Xi.Xixn;
	Xiyy = Xi.Xiyy;
	Xiyn = Xi.Xiyn;
	Xinn = Xi.Xinn;

	Xitt_tt =  0.5 * Xitt * Xitt;
	Xitt_tx =  0.5 * Xitt * Xitx;
	Xitt_ty =  0.5 * Xitt * Xity;
	Xitt_tn =  0.5 * Xitt * Xitn;
	Xitt_xx =  Xitx * Xitx  -  0.5 * Xitt * Xixx;
	Xitt_xy =  Xitx * Xity  -  0.5 * Xitt * Xixy;
	Xitt_xn =  Xitx * Xitn  -  0.5 * Xitt * Xixn;
	Xitt_yy =  Xity * Xity  -  0.5 * Xitt * Xiyy;
	Xitt_yn =  Xity * Xitn  -  0.5 * Xitt * Xiyn;
	Xitt_nn =  Xitn * Xitn  -  0.5 * Xitt * Xinn;

	Xitx_tx = 0.5 * Xitt * Xixx;
	Xitx_ty = 0.5 * Xitt * Xixy;
	Xitx_tn = 0.5 * Xitt * Xixn;
	Xitx_xx = 0.5 * Xitx * Xixx;
	Xitx_xy = 0.5 * Xity * Xixx;
	Xitx_xn = 0.5 * Xitn * Xixx;
	Xitx_yy = Xity * Xixy  -  0.5 * Xitx * Xiyy;
	Xitx_yn = 0.5 * (Xity * Xixn  +  Xitn * Xixy  -  Xitx * Xiyn);
	Xitx_nn = Xitn * Xixn  -  0.5 * Xitx * Xinn;

	Xity_ty = 0.5 * Xitt * Xiyy;
	Xity_tn = 0.5 * Xitt * Xiyn;
	Xity_xx = Xitx * Xixy  -  0.5 * Xity * Xixx;
	Xity_xy = 0.5 * Xitx * Xiyy;
	Xity_xn = 0.5 * (Xitx * Xiyn  +  Xitn * Xixy  -  Xity * Xixn);
	Xity_yy = 0.5 * Xity * Xiyy;
	Xity_yn = 0.5 * Xitn * Xiyy;
	Xity_nn = Xitn * Xiyn  -  0.5 * Xity * Xinn;

	Xitn_tn = 0.5 * Xitt * Xinn;
	Xitn_xx = Xitx * Xixn  -  0.5 * Xitn * Xixx;
	Xitn_xy = 0.5 * (Xity * Xixn  +  Xitx * Xiyn  -  Xitn * Xixy);
	Xitn_xn = 0.5 * Xitx * Xinn;
	Xitn_yy = Xity * Xiyn  -  0.5 * Xitn * Xiyy;
	Xitn_yn = 0.5 * Xity * Xinn;
	Xitn_nn = 0.5 * Xitn * Xinn;

	Xixx_xx = 0.5 * Xixx * Xixx;
	Xixx_xy = 0.5 * Xixx * Xixy;
	Xixx_xn = 0.5 * Xixx * Xixn;
	Xixx_yy = Xixy * Xixy  -  0.5 * Xixx * Xiyy;
	Xixx_yn = Xixy * Xixn  -  0.5 * Xixx * Xiyn;
	Xixx_nn = Xixn * Xixn  -  0.5 * Xixx * Xinn;

	Xixy_xy = 0.5 * Xixx * Xiyy;
	Xixy_xn = 0.5 * Xixx * Xiyn;
	Xixy_yy = 0.5 * Xiyy * Xixy;
	Xixy_yn = 0.5 * Xiyy * Xixn;
	Xixy_nn = Xixn * Xiyn  -  0.5 * Xixy * Xinn;

	Xixn_xn = 0.5 * Xixx * Xinn;
	Xixn_yy = Xixy * Xiyn  -  0.5 * Xiyy * Xixn;
	Xixn_yn = 0.5 * Xixy * Xinn;
	Xixn_nn = 0.5 * Xixn * Xinn;

	Xiyy_yy = 0.5 * Xiyy * Xiyy;
	Xiyy_yn = 0.5 * Xiyy * Xiyn;
	Xiyy_nn = Xiyn * Xiyn  -  0.5 * Xiyy * Xinn;

	Xiyn_yn = 0.5 * Xiyy * Xinn;
	Xiyn_nn = 0.5 * Xiyn * Xinn;

	Xinn_nn = 0.5 * Xinn * Xinn;
}


double_transverse_projection::~double_transverse_projection()
{

}


void double_transverse_projection::double_transverse_project_tensor(precision & Att, precision & Atx, precision & Aty, precision & Atn, precision & Axx, precision & Axy, precision & Axn, precision & Ayy, precision & Ayn, precision & Ann)
{
	// A_pro^{\mu\nu} = Xi^{\mu\nu\alpha\beta} . A_{\alpha\beta}
	precision Att_pro = Xitt_tt * Att  +  Xitt_xx * Axx  +  Xitt_yy * Ayy  +  t4 * Xitt_nn * Ann  -  2. * (Xitt_tx * Atx  +  Xitt_ty * Aty  -  Xitt_xy * Axy  +  t2 * (Xitt_tn * Atn   -  Xitt_xn * Axn  -  Xitt_yn * Ayn));
	precision Atx_pro = Xitt_tx * Att  +  Xitx_xx * Axx  +  Xitx_yy * Ayy  +  t4 * Xitx_nn * Ann  -  2. * (Xitx_tx * Atx  +  Xitx_ty * Aty  -  Xitx_xy * Axy  +  t2 * (Xitx_tn * Atn   -  Xitx_xn * Axn  -  Xitx_yn * Ayn));
	precision Aty_pro = Xitt_ty * Att  +  Xity_xx * Axx  +  Xity_yy * Ayy  +  t4 * Xity_nn * Ann  -  2. * (Xitx_ty * Atx  +  Xity_ty * Aty  -  Xity_xy * Axy  +  t2 * (Xity_tn * Atn   -  Xity_xn * Axn  -  Xity_yn * Ayn));
	precision Atn_pro = Xitt_tn * Att  +  Xitn_xx * Axx  +  Xitn_yy * Ayy  +  t4 * Xitn_nn * Ann  -  2. * (Xitx_tn * Atx  +  Xity_tn * Aty  -  Xitn_xy * Axy  +  t2 * (Xitn_tn * Atn   -  Xitn_xn * Axn  -  Xitn_yn * Ayn));
	precision Axx_pro = Xitt_xx * Att  +  Xixx_xx * Axx  +  Xixx_yy * Ayy  +  t4 * Xixx_nn * Ann  -  2. * (Xitx_xx * Atx  +  Xity_xx * Aty  -  Xixx_xy * Axy  +  t2 * (Xitn_xx * Atn   -  Xixx_xn * Axn  -  Xixx_yn * Ayn));
	precision Axy_pro = Xitt_xy * Att  +  Xixx_xy * Axx  +  Xixy_yy * Ayy  +  t4 * Xixy_nn * Ann  -  2. * (Xitx_xy * Atx  +  Xity_xy * Aty  -  Xixy_xy * Axy  +  t2 * (Xitn_xy * Atn   -  Xixy_xn * Axn  -  Xixy_yn * Ayn));
	precision Axn_pro = Xitt_xn * Att  +  Xixx_xn * Axx  +  Xixn_yy * Ayy  +  t4 * Xixn_nn * Ann  -  2. * (Xitx_xn * Atx  +  Xity_xn * Aty  -  Xixy_xn * Axy  +  t2 * (Xitn_xn * Atn   -  Xixn_xn * Axn  -  Xixn_yn * Ayn));
	precision Ayy_pro = Xitt_yy * Att  +  Xixx_yy * Axx  +  Xiyy_yy * Ayy  +  t4 * Xiyy_nn * Ann  -  2. * (Xitx_yy * Atx  +  Xity_yy * Aty  -  Xixy_yy * Axy  +  t2 * (Xitn_yy * Atn   -  Xixn_yy * Axn  -  Xiyy_yn * Ayn));
	precision Ayn_pro = Xitt_yn * Att  +  Xixx_yn * Axx  +  Xiyy_yn * Ayy  +  t4 * Xiyn_nn * Ann  -  2. * (Xitx_yn * Atx  +  Xity_yn * Aty  -  Xixy_yn * Axy  +  t2 * (Xitn_yn * Atn   -  Xixn_yn * Axn  -  Xiyn_yn * Ayn));
	precision Ann_pro = Xitt_nn * Att  +  Xixx_nn * Axx  +  Xiyy_nn * Ayy  +  t4 * Xinn_nn * Ann  -  2. * (Xitx_nn * Atx  +  Xity_nn * Aty  -  Xixy_nn * Axy  +  t2 * (Xitn_nn * Atn   -  Xixn_nn * Axn  -  Xiyn_nn * Ayn));

	Att = Att_pro;
	Atx = Atx_pro;
	Aty = Aty_pro;
	Atn = Atn_pro;
	Axx = Axx_pro;
	Axy = Axy_pro;
	Axn = Axn_pro;
	Ayy = Ayy_pro;
	Ayn = Ayn_pro;
	Ann = Ann_pro;
}


void double_transverse_projection::test_double_transverse_projector(precision ut, precision ux, precision uy, precision un, precision zt, precision zn)
{
	// test double projector by projecting a random symmetric matrix
	unsigned long int seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);

	precision Att = canonical(generator);
	precision Atx = canonical(generator);
	precision Aty = canonical(generator);
	precision Atn = canonical(generator);
	precision Axx = canonical(generator);
	precision Axy = canonical(generator);
	precision Axn = canonical(generator);
	precision Ayy = canonical(generator);
	precision Ayn = canonical(generator);
	precision Ann = canonical(generator);

	double_transverse_project_tensor(Att, Atx, Aty, Atn, Axx, Axy, Axn, Ayy, Ayn, Ann);

	precision Au0 = fabs(Att * ut  -  Atx * ux  -  Aty * uy  -  t2 * Atn * un);
	precision Au1 = fabs(Atx * ut  -  Axx * ux  -  Axy * uy  -  t2 * Axn * un);
	precision Au2 = fabs(Aty * ut  -  Axy * ux  -  Ayy * uy  -  t2 * Ayn * un);
	precision Au3 = fabs(Atn * ut  -  Axn * ux  -  Ayn * uy  -  t2 * Ann * un);

	precision Az0 = fabs(Att * zt  -  t2 * Atn * zn);
	precision Az1 = fabs(Atx * zt  -  t2 * Axn * zn);
	precision Az2 = fabs(Aty * zt  -  t2 * Ayn * zn);
	precision Az3 = fabs(Atn * zt  -  t2 * Ann * zn);

	precision trA = fabs(Att  -  Axx  -  Ayy  -  t2 * Ann);

	precision Au = fmax(Au0, fmax(Au1, fmax(Au2, Au3)));
	precision Az = fmax(Az0, fmax(Az1, fmax(Az2, Az3)));

    precision eps = 1.e-14;

#ifdef FLAGS
    if(Au > eps)  printf("test_double_transverse_projector flag: A is not orthogonal to u (%.6g)\n", Au);
    if(Az > eps)  printf("test_double_transverse_projector flag: A is not orthogonal to z (%.6g)\n", Az);
    if(trA > eps) printf("test_double_transverse_projector flag: A is not traceless (%.6g)\n", trA);
#endif
}









