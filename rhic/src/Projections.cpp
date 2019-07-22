
#include <stdlib.h>
// #include <stdio.h>
// #include <iostream>
// #include <math.h> // for math functions
// #include <cmath>
#include "../include/Precision.h"
#include "../include/Projections.h"

using namespace std;

double_transverse_projection::double_transverse_projection(precision Xitt_in, precision Xitx_in, precision Xity_in, precision Xitn_in, precision Xixx_in, precision Xixy_in, precision Xixn_in, precision Xiyy_in, precision Xiyn_in, precision Xinn_in)
{
	Xitt = Xitt_in;
	Xitx = Xitx_in;
	Xity = Xity_in;
	Xitn = Xitn_in;
	Xixx = Xixx_in;
	Xixy = Xixy_in;
	Xixn = Xixn_in;
	Xiyy = Xiyy_in;
	Xiyn = Xiyn_in;
	Xinn = Xinn_in;
}


double_transverse_projection::~double_transverse_projection()
{
	// test the projection formula on a random matrix

	// Xi^{\mu\nu\alpha\beta} (all upper indices)

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

	//Xitx_tt = Xitt_tx;
	Xitx_tx = 0.5 * Xitt * Xixx;
	Xitx_ty = 0.5 * Xitt * Xixy;
	Xitx_tn = 0.5 * Xitt * Xixn;
	Xitx_xx = 0.5 * Xitx * Xixx;
	Xitx_xy = 0.5 * Xity * Xixx;
	Xitx_xn = 0.5 * Xitn * Xixx;
	Xitx_yy = Xity * Xixy  -  0.5 * Xitx * Xiyy;
	Xitx_yn = 0.5 * (Xity * Xixn  +  Xitn * Xixy  -  Xitx * Xiyn);
	Xitx_nn = Xitn * Xixn  -  0.5 * Xitx * Xinn;

	//Xity_tt = Xitt_ty;
	Xity_tx = 0.5 * Xitt * Xixy;
	//Xity_tx = Xitx_ty;
	Xity_ty = 0.5 * Xitt * Xiyy;
	Xity_tn = 0.5 * Xitt * Xiyn;
	Xity_xx = Xitx * Xixy  -  0.5 * Xity * Xixx;
	Xity_xy = 0.5 * Xitx * Xiyy;
	Xity_xn = 0.5 * (Xitx * Xiyn  +  Xitn * Xixy  -  Xity * Xixn);
	Xity_yy = 0.5 * Xity * Xiyy;
	Xity_yn = 0.5 * Xitn * Xiyy;
	Xity_nn = Xitn * Xiyn  -  0.5 * Xity * Xinn;

	Xitn_tt = Xitt_tn;






	precision Xixx_xx = 0.5 * Xixx * Xixx;
	precision Xiyy_yy = 0.5 * Xiyy * Xiyy;
	precision Xinn_nn = 0.5 * Xinn * Xinn;
}









