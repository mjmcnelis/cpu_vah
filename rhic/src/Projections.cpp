#include <stdlib.h>
#include <stdio.h>
#include <random>
#include <chrono>
#include "../include/Precision.h"
#include "../include/Projections.h"
using namespace std;

double canonical(default_random_engine & generator)
{
  // random number between [0,1)
  return generate_canonical<double, numeric_limits<double>::digits>(generator);
}


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

    precision eps = 1.e-14;

    if(Au > eps)  printf("test_transverse_projector error: A is not orthogonal to u (%.6g)\n", Au);
    if(Az > eps)  printf("test_transverse_projector error: A is not orthogonal to z (%.6g)\n", Az);
}


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
	precision Att_pro = Xitt_tt * Att  +  Xitt_xx * Axx  +  Xitt_yy * Ayy  +  t4 * Xitt_nn * Ann  -  2.0 * (Xitt_tx * Atx  +  Xitt_ty * Aty  -  Xitt_xy * Axy  +  t2 * (Xitt_tn * Atn   -  Xitt_xn * Axn  -  Xitt_yn * Ayn));
	precision Atx_pro = Xitt_tx * Att  +  Xitx_xx * Axx  +  Xitx_yy * Ayy  +  t4 * Xitx_nn * Ann  -  2.0 * (Xitx_tx * Atx  +  Xitx_ty * Aty  -  Xitx_xy * Axy  +  t2 * (Xitx_tn * Atn   -  Xitx_xn * Axn  -  Xitx_yn * Ayn));
	precision Aty_pro = Xitt_ty * Att  +  Xity_xx * Axx  +  Xity_yy * Ayy  +  t4 * Xity_nn * Ann  -  2.0 * (Xitx_ty * Atx  +  Xity_ty * Aty  -  Xity_xy * Axy  +  t2 * (Xity_tn * Atn   -  Xity_xn * Axn  -  Xity_yn * Ayn));
	precision Atn_pro = Xitt_tn * Att  +  Xitn_xx * Axx  +  Xitn_yy * Ayy  +  t4 * Xitn_nn * Ann  -  2.0 * (Xitx_tn * Atx  +  Xity_tn * Aty  -  Xitn_xy * Axy  +  t2 * (Xitn_tn * Atn   -  Xitn_xn * Axn  -  Xitn_yn * Ayn));
	precision Axx_pro = Xitt_xx * Att  +  Xixx_xx * Axx  +  Xixx_yy * Ayy  +  t4 * Xixx_nn * Ann  -  2.0 * (Xitx_xx * Atx  +  Xity_xx * Aty  -  Xixx_xy * Axy  +  t2 * (Xitn_xx * Atn   -  Xixx_xn * Axn  -  Xixx_yn * Ayn));
	precision Axy_pro = Xitt_xy * Att  +  Xixx_xy * Axx  +  Xixy_yy * Ayy  +  t4 * Xixy_nn * Ann  -  2.0 * (Xitx_xy * Atx  +  Xity_xy * Aty  -  Xixy_xy * Axy  +  t2 * (Xitn_xy * Atn   -  Xixy_xn * Axn  -  Xixy_yn * Ayn));
	precision Axn_pro = Xitt_xn * Att  +  Xixx_xn * Axx  +  Xixn_yy * Ayy  +  t4 * Xixn_nn * Ann  -  2.0 * (Xitx_xn * Atx  +  Xity_xn * Aty  -  Xixy_xn * Axy  +  t2 * (Xitn_xn * Atn   -  Xixn_xn * Axn  -  Xixn_yn * Ayn));
	precision Ayy_pro = Xitt_yy * Att  +  Xixx_yy * Axx  +  Xiyy_yy * Ayy  +  t4 * Xiyy_nn * Ann  -  2.0 * (Xitx_yy * Atx  +  Xity_yy * Aty  -  Xixy_yy * Axy  +  t2 * (Xitn_yy * Atn   -  Xixn_yy * Axn  -  Xiyy_yn * Ayn));
	precision Ayn_pro = Xitt_yn * Att  +  Xixx_yn * Axx  +  Xiyy_yn * Ayy  +  t4 * Xiyn_nn * Ann  -  2.0 * (Xitx_yn * Atx  +  Xity_yn * Aty  -  Xixy_yn * Axy  +  t2 * (Xitn_yn * Atn   -  Xixn_yn * Axn  -  Xiyn_yn * Ayn));
	precision Ann_pro = Xitt_nn * Att  +  Xixx_nn * Axx  +  Xiyy_nn * Ayy  +  t4 * Xinn_nn * Ann  -  2.0 * (Xitx_nn * Atx  +  Xity_nn * Aty  -  Xixy_nn * Axy  +  t2 * (Xitn_nn * Atn   -  Xixn_nn * Axn  -  Xiyn_nn * Ayn));

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

    if(Au > eps)  printf("test_double_transverse_projector error: A is not orthogonal to u (%.6g)\n", Au);
    if(Az > eps)  printf("test_double_transverse_projector error: A is not orthogonal to z (%.6g)\n", Az);
    if(trA > eps) printf("test_double_transverse_projector error: A is not traceless (%.6g)\n", trA);
}









