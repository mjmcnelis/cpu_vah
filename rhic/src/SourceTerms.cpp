#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include "../include/SourceTerms.h"
#include "../include/FluxTerms.h"
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"
#include "../include/EquationOfState.h" // for bulk terms
#include "../include/TransportCoefficients.h"
#include "../include/Projections.h"

using namespace std;

#define A_1 -13.77
#define A_2 27.55
#define A_3 13.45

#define LAMBDA_1 0.9
#define LAMBDA_2 0.25
#define LAMBDA_3 0.9
#define LAMBDA_4 0.22

#define SIGMA_1 0.025
#define SIGMA_2 0.13
#define SIGMA_3 0.0025
#define SIGMA_4 0.022


inline precision zetabar(precision T)
{
	precision x = T / 1.01355;

	if(x > 1.05)
	{
		return LAMBDA_1*exp(-(x-1)/SIGMA_1) + LAMBDA_2*exp(-(x-1)/SIGMA_2)+0.001;
	}
	else if(x < 0.995)
	{
		return LAMBDA_3*exp((x-1)/SIGMA_3)+ LAMBDA_4*exp((x-1)/SIGMA_4)+0.03;
	}
	else
	{
		return A_1*x*x + A_2*x - A_3;
	}
}


inline precision central_derivative(const precision * const __restrict__ f, int n, precision dx)
{
	return (f[n + 1] - f[n]) / (2. * dx);		// f[n] = fm  |	 f[n+1] = fp
}

/*
void source_terms(precision * const __restrict__ S, const precision * const __restrict__ q, precision e, precision t, const precision * const __restrict__ qi1, const precision * const __restrict__ qj1, const precision * const __restrict__ qk1, const precision * const __restrict__ e1, const precision * const __restrict__ ui1, const precision * const __restrict__ uj1, const precision * const __restrict__ uk1, precision ux, precision uy, precision un, precision ux_p, precision uy_p, precision un_p, precision dt, precision dx, precision dy, precision dn, precision etabar)
{
	precision t2 = t * t;
	precision t4 = t2 * t2;

	precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t2 * un * un);
	precision vx = ux / ut;				// spatial velocity
	precision vy = uy / ut;
	precision vn = un / ut;

	precision utperp  = sqrt(1.  +  ux * ux  +  uy * uy);
	precision zt = t * un / utperp;		// longitudinal basis vector
	precision zn = ut / t / utperp;

	precision zt2  = zt * zt;			// other useful expressions
	precision ztzn = zt * zn;
	precision zn2  = zn * zn;
	precision t2zn = t2 * zn;
	precision ut2  = ut * ut;
	precision un2  = un * un;
	precision utun = ut * un;
	precision tun  = t  * un;
	precision unzn = un * zn;
	precision utperp2 = utperp * utperp;

	precision T = effectiveTemperature(e);
	precision p = equilibriumPressure(e);

	// conserved variables
	precision ttt = q[0];
	precision ttx = q[1];
	precision tty = q[2];
	precision ttn = q[3];
	precision pl  = q[4];

	int a = 5;	// counter

#if (PT_MATCHING == 1)
	precision pt  = q[a];	a++;
#else
	precision pt  = (e - pl) / 2.;
#endif
#ifdef PIMUNU
	precision pitt = q[a];	a++;
	precision pitx = q[a];	a++;
	precision pity = q[a];	a++;
	precision pitn = q[a];	a++;
	precision pixx = q[a];	a++;
	precision pixy = q[a];	a++;
	precision pixn = q[a];	a++;
	precision piyy = q[a];	a++;
	precision piyn = q[a];	a++;
	precision pinn = q[a];	a++;
#else
	precision pitt = 0, pitx = 0, pity = 0, pitn = 0, pixx = 0, pixy = 0, pixn = 0, piyy = 0, piyn = 0, pinn = 0;
#endif
#ifdef WTZMU
	precision WtTz = q[a];	a++;
	precision WxTz = q[a];	a++;
	precision WyTz = q[a];	a++;
	precision WnTz = q[a];
#else
	precision WtTz = 0, WxTz = 0, WyTz = 0, WnTz = 0;
#endif

// relaxation times and transport coefficients
//-------------------------------------------------
#ifdef CONFORMAL_EOS
	precision taupiInv = T / (5. * etabar);
	precision taubulkInv = 0;
#else
	precision taupiInv = 1.;	// fill in quasiparticle later
	precision taubulkInv = 1.;
#endif

	// anisotropic transport coefficients
	transport_coefficients aniso;
	aniso.compute_transport_coefficients(e, pl, pt);

	// pl transport coefficients
	precision zeta_LL = aniso.zeta_LL;
	precision zeta_TL = aniso.zeta_TL;
	precision lambda_WuL = aniso.lambda_WuL;
	precision lambda_WTL = aniso.lambda_WTL;
	precision lambda_piL = aniso.lambda_piL;

	// pt transport coefficients
#if (PT_MATCHING == 1)
	precision zeta_LT = aniso.zeta_LT;
	precision zeta_TT = aniso.zeta_TT;
	precision lambda_WuT = aniso.lambda_WuT;
	precision lambda_WTT = aniso.lambda_WTT;
	precision lambda_piT = aniso.lambda_piT;
#endif

	// primary variable derivatives
#if (PT_MATCHING == 0)
	precision de_dx = (e1[1] - e1[0]) / (2. * dx);
	precision de_dy = (e1[3] - e1[2]) / (2. * dy);
	precision de_dn = (e1[5] - e1[4]) / (2. * dn);
#endif
	// pl derivatives
	int n = 8;
	precision dpl_dx = central_derivative(qi1, n, dx);
	precision dpl_dy = central_derivative(qj1, n, dy);
	precision dpl_dn = central_derivative(qk1, n, dn);		n += 2;

#if (PT_MATCHING == 1)	// pt derivatives
	precision dpt_dx = central_derivative(qi1, n, dx);
	precision dpt_dy = central_derivative(qj1, n, dy);
	precision dpt_dn = central_derivative(qk1, n, dn);		n += 2;
#else
	precision dpt_dx = (de_dx  -  dpl_dx) / 2.;
	precision dpt_dy = (de_dy  -  dpl_dy) / 2.;
	precision dpt_dn = (de_dn  -  dpl_dn) / 2.;
#endif

#ifdef PIMUNU	// pimunu derivatives
	precision dpitt_dx = central_derivative(qi1, n, dx);
	precision dpitt_dy = central_derivative(qj1, n, dy);
	precision dpitt_dn = central_derivative(qk1, n, dn);	n += 2;

	precision dpitx_dx = central_derivative(qi1, n, dx);
	precision dpitx_dy = central_derivative(qj1, n, dy);
	precision dpitx_dn = central_derivative(qk1, n, dn);	n += 2;

	precision dpity_dx = central_derivative(qi1, n, dx);
	precision dpity_dy = central_derivative(qj1, n, dy);
	precision dpity_dn = central_derivative(qk1, n, dn);	n += 2;

	precision dpitn_dx = central_derivative(qi1, n, dx);
	precision dpitn_dy = central_derivative(qj1, n, dy);
	precision dpitn_dn = central_derivative(qk1, n, dn);	n += 2;

	precision dpixx_dx = central_derivative(qi1, n, dx);	n += 2;

	precision dpixy_dx = central_derivative(qi1, n, dx);
	precision dpixy_dy = central_derivative(qj1, n, dy);	n += 2;

	precision dpixn_dx = central_derivative(qi1, n, dx);
	precision dpixn_dn = central_derivative(qk1, n, dn);	n += 2;

	precision dpiyy_dy = central_derivative(qj1, n, dy);	n += 2;

	precision dpiyn_dy = central_derivative(qj1, n, dy);
	precision dpiyn_dn = central_derivative(qk1, n, dn);	n += 2;

	precision dpinn_dn = central_derivative(qk1, n, dn);	n += 2;

#else
	precision dpitt_dx = 0, dpitt_dy = 0, dpitt_dn = 0, dpitx_dx = 0, dpitx_dy = 0, dpitx_dn = 0, dpity_dx = 0, dpity_dy = 0, dpity_dn = 0, dpitn_dx = 0, dpitn_dy = 0, dpitn_dn = 0, dpixx_dx = 0, dpixy_dx = 0, dpixy_dy = 0, dpixn_dx = 0, dpixn_dn = 0, dpiyy_dy = 0, dpiyn_dy = 0, dpiyn_dn = 0, dpinn_dn = 0;
#endif
#ifdef WTZMU
	precision dWtTz_dx = central_derivative(qi1, n, dx);
	precision dWtTz_dy = central_derivative(qj1, n, dy);
	precision dWtTz_dn = central_derivative(qk1, n, dn);	n += 2;

	precision dWxTz_dx = central_derivative(qi1, n, dx);
	precision dWxTz_dy = central_derivative(qj1, n, dy);
	precision dWxTz_dn = central_derivative(qk1, n, dn);	n += 2;

	precision dWyTz_dx = central_derivative(qi1, n, dx);
	precision dWyTz_dy = central_derivative(qj1, n, dy);
	precision dWyTz_dn = central_derivative(qk1, n, dn);	n += 2;

	precision dWnTz_dx = central_derivative(qi1, n, dx);
	precision dWnTz_dy = central_derivative(qj1, n, dy);
	precision dWnTz_dn = central_derivative(qk1, n, dn);

#else
	precision dWtTz_dx = 0, dWtTz_dy = 0, dWtTz_dn = 0, dWxTz_dx = 0, dWxTz_dy = 0, dWxTz_dn = 0, dWyTz_dx = 0, dWyTz_dy = 0, dWyTz_dn = 0, dWnTz_dx = 0, dWnTz_dy = 0, dWnTz_dn = 0;
#endif

	// fluid velocity derivatives
	n = 0;
	precision dux_dt = (ux - ux_p) / dt;
	precision dux_dx = central_derivative(ui1, n, dx);
	precision dux_dy = central_derivative(uj1, n, dy);
	precision dux_dn = central_derivative(uk1, n, dn);	n += 2;

	precision duy_dt = (uy - uy_p) / dt;
	precision duy_dx = central_derivative(ui1, n, dx);
	precision duy_dy = central_derivative(uj1, n, dy);
	precision duy_dn = central_derivative(uk1, n, dn);	n += 2;

	precision dun_dt = (un - un_p) / dt;
	precision dun_dx = central_derivative(ui1, n, dx);
	precision dun_dy = central_derivative(uj1, n, dy);
	precision dun_dn = central_derivative(uk1, n, dn);


	// use chain rule for ut derivatives
	precision dut_dt = vx * dux_dt  +  vy * duy_dt  +  t2 * vn * dun_dt  +  t * vn * un;
	precision dut_dx = vx * dux_dx  +  vy * duy_dx  +  t2 * vn * dun_dx;
	precision dut_dy = vx * dux_dy  +  vy * duy_dy  +  t2 * vn * dun_dy;
	precision dut_dn = vx * dux_dn  +  vy * duy_dn  +  t2 * vn * dun_dn;

	// radial velocity derivatives
	precision duT2_dt = ux * dux_dt  +  uy * duy_dt;
	precision duT2_dx = ux * dux_dx  +  uy * duy_dx;
	precision duT2_dy = ux * dux_dy  +  uy * duy_dy;
	precision duT2_dn = ux * dux_dn  +  uy * duy_dn;

	// scalar, longitudinal and transverse expansion rates: theta = D_\mu u^\mu, thetaL = z_\mu Dz u^\mu, thetaT = NablaT_\mu u^\mu
	precision theta = dut_dt  +  dux_dx  +  duy_dy  +  dun_dn  +  ut / t;
	precision thetaL = - zt2 * dut_dt  +  t2 * zn2 * dun_dn  +  ztzn * (t2 * dun_dt  -  dut_dn)  +  t * zn2 * ut;
	precision thetaT = theta  -  thetaL;

	// spatial velocity derivatives and divergence
	precision dvx_dx = (dux_dx  -  vx * dut_dx) / ut;
	precision dvy_dy = (duy_dy  -  vy * dut_dy) / ut;
	precision dvn_dn = (dun_dn  -  vn * dut_dn) / ut;

	precision div_v = dvx_dx  +  dvy_dy  +  dvn_dn;

	// other spatial velocity derivatives (get rid of chain rule later)
	precision dvx_dn = (dux_dn  -  vx * dut_dn) / ut;
	precision dvy_dn = (duy_dn  -  vy * dut_dn) / ut;

	// longitudinal vector derivatives
	precision dzt_dt = t * (dun_dt  -  un * duT2_dt / utperp2) / utperp  +  zt / t;
	precision dzt_dx = t * (dun_dx  -  un * duT2_dx / utperp2) / utperp;
	precision dzt_dy = t * (dun_dy  -  un * duT2_dy / utperp2) / utperp;
	precision dzt_dn = t * (dun_dn  -  un * duT2_dn / utperp2) / utperp;

	precision dzn_dt = (dut_dt  -  ut * duT2_dt / utperp2) / (t * utperp)  -  zn / t;
	precision dzn_dx = (dut_dx  -  ut * duT2_dx / utperp2) / (t * utperp);
	precision dzn_dy = (dut_dy  -  ut * duT2_dy / utperp2) / (t * utperp);
	precision dzn_dn = (dut_dn  -  ut * duT2_dn / utperp2) / (t * utperp);


#if (NUMBER_OF_RESIDUAL_CURRENTS != 0)
	transverse_projection Xi(ut, ux, uy, un, zt, zn, t2);	// Xi^{\mu\nu}
	double_transverse_projection Xi_2(Xi, t2, t4);			// Xi^{\mu\nu\alpha\beta}

	//double_transverse_projection Xi_2_sigma(Xi, t2, t4);

	// project residual shear stress that appear in source terms
#ifdef PIMUNU
	Xi_2.double_transverse_project_tensor(pitt, pitx, pity, pitn, pixx, pixy, pixn, piyy, piyn, pinn);
#endif
#ifdef WTZMU
	Xi.transverse_project_vector(WtTz, WxTz, WyTz, WnTz);
#endif

	// acceleration = D u^\mu
	precision at = ut * dut_dt  +  ux * dut_dx  +  uy * dut_dy  +  un * dut_dn  +  t * un2;
	precision ax = ut * dux_dt  +  ux * dux_dx  +  uy * dux_dy  +  un * dux_dn;
	precision ay = ut * duy_dt  +  ux * duy_dx  +  uy * duy_dy  +  un * duy_dn;
	precision an = ut * dun_dt  +  ux * dun_dx  +  uy * dun_dy  +  un * dun_dn  +  2. * utun / t;

	// covariant time derivatives of z = D z^\mu
	precision D_zt = ut * dzt_dt  +  ux * dzt_dx  +  uy * dzt_dy  +  un * dzt_dn  +  t * unzn;
	precision D_zn = ut * dzn_dt  +  ux * dzn_dx  +  uy * dzn_dy  +  un * dzn_dn  +  (ut * zn  +  un * zt) / t;

	// longitudinal covariant derivative of u = Dz u^\mu
	precision Dz_ut = - zt * dut_dt  -  zn * dut_dn  -  t * unzn;
	precision Dz_ux = - zt * dux_dt  -  zn * dux_dn;
	precision Dz_uy = - zt * duy_dt  -  zn * duy_dn;
	precision Dz_un = - zt * dun_dt  -  zn * dun_dn  -  (ut * zn  +  un * zt) / t;

	// covariant transverse derivative of u contracted with z = z_\nu NablaT^\mu u^\nu (first compute z_\nu D^\mu u^\nu )
	precision z_NabTt_u =   zt * dut_dt  -  t2zn * dun_dt  -  t * unzn;
	precision z_NabTx_u = -(zt * dut_dx  -  t2zn * dun_dx);
	precision z_NabTy_u = -(zt * dut_dy  -  t2zn * dun_dy);
	precision z_NabTn_u = -(zt * dut_dn  -  t2zn * dun_dn  +  t * (un * zt  -  ut * zn)) / t2;
	Xi.transverse_project_vector(z_NabTt_u, z_NabTx_u, z_NabTy_u, z_NabTn_u);

	// transverse shear velocity tensor = sigmaT^{\mu\nu} (first compute D^{(\mu) u^{\nu)} )
	precision sTtt = dut_dt;
	precision sTtx = (dux_dt  -  dut_dx) / 2.;
	precision sTty = (duy_dt  -  dut_dy) / 2.;
	precision sTtn = (dun_dt  +  un / t  -  (dut_dn  +  tun) / t2) / 2.;  // I don't see anything wrong here...
	precision sTxx = - dux_dx;
	precision sTxy = - (dux_dy  +  duy_dx) / 2.;
	precision sTxn = - (dun_dx  +  dux_dn / t2) / 2.;
	precision sTyy = - duy_dy;
	precision sTyn = - (dun_dy  +  duy_dn / t2) / 2.;
	precision sTnn = - (dun_dn  +  ut / t) / t2;
	Xi_2.double_transverse_project_tensor(sTtt, sTtx, sTty, sTtn, sTxx, sTxy, sTxn, sTyy, sTyn, sTnn);
#endif

	// L^munu components and derivatives
	precision dp  = pl - pt;
	precision Ltt = dp * zt2;
	precision Ltn = dp * ztzn;
	precision Lnn = dp * zn2;

	precision dLtt_dx = (dpl_dx - dpt_dx) * zt2  +  2. * dp * zt * dzt_dx;
	precision dLtt_dy = (dpl_dy - dpt_dy) * zt2  +  2. * dp * zt * dzt_dy;
	precision dLtt_dn = (dpl_dn - dpt_dn) * zt2  +  2. * dp * zt * dzt_dn;

	precision dLtn_dx = (dpl_dx - dpt_dx) * ztzn  +  dp * (dzt_dx * zn  +  dzn_dx * zt);
	precision dLtn_dy = (dpl_dy - dpt_dy) * ztzn  +  dp * (dzt_dy * zn  +  dzn_dy * zt);
	precision dLtn_dn = (dpl_dn - dpt_dn) * ztzn  +  dp * (dzt_dn * zn  +  dzn_dn * zt);

	precision dLnn_dn = (dpl_dn - dpt_dn) * zn2  +  2. * dp * zn * dzn_dn;

#ifdef PIMUNU
	// piT transport coefficients
	precision eta_T = aniso.eta_T;
	precision delta_pipi = aniso.delta_pipi;
	precision tau_pipi = aniso.tau_pipi;
	precision lambda_pipi = aniso.lambda_pipi;
	precision lambda_Wupi = aniso.lambda_Wupi;
	precision lambda_WTpi = aniso.lambda_WTpi;

	precision pi_sT = pitt * sTtt  +  pixx * sTxx  +  piyy * sTyy  +  t4 * pinn * sTnn  +  2. * (pixy * sTxy  -  pitx * sTtx  -  pity * sTty  +  t2 * (pixn * sTxn  +  piyn * sTyn  -  pitn * sTtn));

	// pi_T gradient terms (tensors are symmetric but (2,4,5,7,8) not double transverse projected until end

	// 2 . \eta_T . \sigma_T^{\mu\nu}
	precision Itt_pi1 = 2. * eta_T * sTtt;
	precision Itx_pi1 = 2. * eta_T * sTtx;
	precision Ity_pi1 = 2. * eta_T * sTty;
	precision Itn_pi1 = 2. * eta_T * sTtn;
	precision Ixx_pi1 = 2. * eta_T * sTxx;
	precision Ixy_pi1 = 2. * eta_T * sTxy;
	precision Ixn_pi1 = 2. * eta_T * sTxn;
	precision Iyy_pi1 = 2. * eta_T * sTyy;
	precision Iyn_pi1 = 2. * eta_T * sTyn;
	precision Inn_pi1 = 2. * eta_T * sTnn;

	// 2 . WTz^{(\mu} . \dot{z}^{\nu)}
#ifdef WTZMU
	precision Itt_pi2 = 2. * WtTz * D_zt;
	precision Itx_pi2 = WxTz * D_zt;
	precision Ity_pi2 = WyTz * D_zt;
	precision Itn_pi2 = WtTz * D_zn  +  WnTz * D_zt;
	precision Ixn_pi2 = WxTz * D_zn;
	precision Iyn_pi2 = WyTz * D_zn;
	precision Inn_pi2 = 2. * WnTz * D_zn;
#else
	precision Itt_pi2 = 0, Itx_pi2 = 0, Ity_pi2 = 0, Itn_pi2 = 0, Ixn_pi2 = 0, Iyn_pi2 = 0, Inn_pi2 = 0;
#endif

	// \delta^\pi_\pi . pi_T^{\mu\nu} . \theta_T
	precision Itt_pi3 = delta_pipi * pitt * thetaT;
	precision Itx_pi3 = delta_pipi * pitx * thetaT;
	precision Ity_pi3 = delta_pipi * pity * thetaT;
	precision Itn_pi3 = delta_pipi * pitn * thetaT;
	precision Ixx_pi3 = delta_pipi * pixx * thetaT;
	precision Ixy_pi3 = delta_pipi * pixy * thetaT;
	precision Ixn_pi3 = delta_pipi * pixn * thetaT;
	precision Iyy_pi3 = delta_pipi * piyy * thetaT;
	precision Iyn_pi3 = delta_pipi * piyn * thetaT;
	precision Inn_pi3 = delta_pipi * pinn * thetaT;

	// \tau^\pi_\pi . \pi_T^{\alpha (\mu} . \sigma_T^{\nu)}_\alpha
	precision Itt_pi4 = tau_pipi * (pitt * sTtt  -  pitx * sTtx  -  pity * sTty  -  t2 * pitn * sTtn);
	precision Itx_pi4 = tau_pipi * (pitt * sTtx  +  pitx * sTtt  -  pitx * sTxx  -  pixx * sTtx  -  pity * sTxy  -  pixy * sTty  -  t2 * (pitn * sTxn  +  pixn * sTtn)) / 2.;
	precision Ity_pi4 = tau_pipi * (pitt * sTty  +  pity * sTtt  -  pitx * sTxy  -  pixy * sTtx  -  pity * sTyy  -  piyy * sTty  -  t2 * (pitn * sTyn  +  piyn * sTtn)) / 2.;
	precision Itn_pi4 = tau_pipi * (pitt * sTtn  +  pitn * sTtt  -  pitx * sTxn  -  pixn * sTtx  -  pity * sTyn  -  piyn * sTty  -  t2 * (pitn * sTnn  +  pinn * sTtn)) / 2.;
	precision Ixx_pi4 = tau_pipi * (pitx * sTtx  -  pixx * sTxx  -  pixy * sTxy  -  t2 * pixn * sTxn);
	precision Ixy_pi4 = tau_pipi * (pitx * sTty  +  pity * sTtx  -  pixx * sTxy  -  pixy * sTxx  -  pixy * sTyy  -  piyy * sTxy  -  t2 * (pixn * sTyn  +  piyn * sTxn)) / 2.;
	precision Ixn_pi4 = tau_pipi * (pitx * sTtn  +  pitn * sTtx  -  pixx * sTxn  -  pixn * sTxx  -  pixy * sTyn  -  piyn * sTxy  -  t2 * (pixn * sTnn  +  pinn * sTxn)) / 2.;
	precision Iyy_pi4 = tau_pipi * (pity * sTty  -  pixy * sTxy  -  piyy * sTyy  -  t2 * piyn * sTyn);
	precision Iyn_pi4 = tau_pipi * (pity * sTtn  +  pitn * sTty  -  pixy * sTxn  -  pixn * sTxy  -  piyy * sTyn  -  piyn * sTyy  -  t2 * (piyn * sTnn  +  pinn * sTyn)) / 2.;
	precision Inn_pi4 = tau_pipi * (pitn * sTtn  -  pixn * sTxn  -  piyn * sTyn  -  t2 * pinn * sTnn);

	// ignore the vorticity term for now
	precision Itt_pi5 = 0;
	precision Itx_pi5 = 0;
	precision Ity_pi5 = 0;
	precision Itn_pi5 = 0;
	precision Ixx_pi5 = 0;
	precision Ixy_pi5 = 0;
	precision Ixn_pi5 = 0;
	precision Iyy_pi5 = 0;
	precision Iyn_pi5 = 0;
	precision Inn_pi5 = 0;

	// \lambda^\pi_\pi . \pi_T^{\mu\nu} . \theta_L
	precision Itt_pi6 = lambda_pipi * pitt * thetaL;
	precision Itx_pi6 = lambda_pipi * pitx * thetaL;
	precision Ity_pi6 = lambda_pipi * pity * thetaL;
	precision Itn_pi6 = lambda_pipi * pitn * thetaL;
	precision Ixx_pi6 = lambda_pipi * pixx * thetaL;
	precision Ixy_pi6 = lambda_pipi * pixy * thetaL;
	precision Ixn_pi6 = lambda_pipi * pixn * thetaL;
	precision Iyy_pi6 = lambda_pipi * piyy * thetaL;
	precision Iyn_pi6 = lambda_pipi * piyn * thetaL;
	precision Inn_pi6 = lambda_pipi * pinn * thetaL;

#ifdef WTZMU
	// \lambda_Wu^\pi . WTz^{(\mu} . Dz u^{\nu)}
	precision Itt_pi7 = lambda_Wupi * (WtTz * Dz_ut);
	precision Itx_pi7 = lambda_Wupi * (WtTz * Dz_ux  +  WxTz * Dz_ut) / 2.;
	precision Ity_pi7 = lambda_Wupi * (WtTz * Dz_uy  +  WyTz * Dz_ut) / 2.;
	precision Itn_pi7 = lambda_Wupi * (WtTz * Dz_un  +  WnTz * Dz_ut) / 2.;
	precision Ixx_pi7 = lambda_Wupi * (WxTz * Dz_ux);
	precision Ixy_pi7 = lambda_Wupi * (WxTz * Dz_uy  +  WyTz * Dz_ux) / 2.;
	precision Ixn_pi7 = lambda_Wupi * (WxTz * Dz_un  +  WnTz * Dz_ux) / 2.;
	precision Iyy_pi7 = lambda_Wupi * (WyTz * Dz_uy);
	precision Iyn_pi7 = lambda_Wupi * (WyTz * Dz_un  +  WnTz * Dz_uy) / 2.;
	precision Inn_pi7 = lambda_Wupi * (WnTz * Dz_un);

	// \lambda_WT^\pi . WTz^{(\mu} . z_\alpha . \Nabla_T^{\nu)} . u^\alpha
	precision Itt_pi8 = lambda_WTpi * (WtTz * z_NabTt_u);
	precision Itx_pi8 = lambda_WTpi * (WtTz * z_NabTx_u  +  WxTz * z_NabTt_u) / 2.;
	precision Ity_pi8 = lambda_WTpi * (WtTz * z_NabTy_u  +  WyTz * z_NabTt_u) / 2.;
	precision Itn_pi8 = lambda_WTpi * (WtTz * z_NabTn_u  +  WnTz * z_NabTt_u) / 2.;
	precision Ixx_pi8 = lambda_WTpi * (WxTz * Dz_ux);
	precision Ixy_pi8 = lambda_WTpi * (WxTz * z_NabTy_u  +  WyTz * z_NabTx_u) / 2.;
	precision Ixn_pi8 = lambda_WTpi * (WxTz * z_NabTn_u  +  WnTz * z_NabTx_u) / 2.;
	precision Iyy_pi8 = lambda_WTpi * (WyTz * z_NabTy_u);
	precision Iyn_pi8 = lambda_WTpi * (WyTz * z_NabTn_u  +  WnTz * z_NabTy_u) / 2.;
	precision Inn_pi8 = lambda_WTpi * (WnTz * z_NabTn_u);
#else
	precision Itt_pi7 = 0, Itx_pi7 = 0, Ity_pi7 = 0, Itn_pi7 = 0, Ixx_pi7 = 0, Ixy_pi7 = 0, Ixn_pi7 = 0, Iyy_pi7 = 0, Iyn_pi7 = 0, Inn_pi7 = 0;
	precision Itt_pi8 = 0, Itx_pi8 = 0, Ity_pi8 = 0, Itn_pi8 = 0, Ixx_pi8 = 0, Ixy_pi8 = 0, Ixn_pi8 = 0, Iyy_pi8 = 0, Iyn_pi8 = 0, Inn_pi8 = 0;
#endif

	precision Itt_pi = Itt_pi1  -  Itt_pi2  -  Itt_pi3  -  Itt_pi4  +  Itt_pi5  +  Itt_pi6  -  Itt_pi7  +  Itt_pi8;
	precision Itx_pi = Itx_pi1  -  Itx_pi2  -  Itx_pi3  -  Itx_pi4  +  Itx_pi5  +  Itx_pi6  -  Itx_pi7  +  Itx_pi8;
	precision Ity_pi = Ity_pi1  -  Ity_pi2  -  Ity_pi3  -  Ity_pi4  +  Ity_pi5  +  Ity_pi6  -  Ity_pi7  +  Ity_pi8;
	precision Itn_pi = Itn_pi1  -  Itn_pi2  -  Itn_pi3  -  Itn_pi4  +  Itn_pi5  +  Itn_pi6  -  Itn_pi7  +  Itn_pi8;
	precision Ixx_pi = Ixx_pi1  		    -  Ixx_pi3  -  Ixx_pi4  +  Ixx_pi5  +  Ixx_pi6  -  Ixx_pi7  +  Ixx_pi8;
	precision Ixy_pi = Ixy_pi1  			-  Ixy_pi3  -  Ixy_pi4  +  Ixy_pi5  +  Ixy_pi6  -  Ixy_pi7  +  Ixy_pi8;
	precision Ixn_pi = Ixn_pi1  -  Ixn_pi2  -  Ixn_pi3  -  Ixn_pi4  +  Ixn_pi5  +  Ixn_pi6  -  Ixn_pi7  +  Ixn_pi8;
	precision Iyy_pi = Iyy_pi1              -  Iyy_pi3  -  Iyy_pi4  +  Iyy_pi5  +  Iyy_pi6  -  Iyy_pi7  +  Iyy_pi8;
	precision Iyn_pi = Iyn_pi1  -  Iyn_pi2  -  Iyn_pi3  -  Iyn_pi4  +  Iyn_pi5  +  Iyn_pi6  -  Iyn_pi7  +  Iyn_pi8;
	precision Inn_pi = Inn_pi1  -  Inn_pi2  -  Inn_pi3  -  Inn_pi4  +  Inn_pi5  +  Inn_pi6  -  Inn_pi7  +  Inn_pi8;
	Xi_2.double_transverse_project_tensor(Itt_pi, Itx_pi, Ity_pi, Itn_pi, Ixx_pi, Ixy_pi, Ixn_pi, Iyy_pi, Iyn_pi, Inn_pi);

	// Christofel terms: G_\pi^{\mu\nu} = 2 . u^\alpha . \Gamma^{(\mu}_{\alpha\beta} . \pi_T^{\beta\nu)}
	precision Gtt_pi = 2. * tun * pitn;
	precision Gtx_pi = tun * pixn;
	precision Gty_pi = tun * piyn;
	precision Gtn_pi = tun * pinn  +  (ut * pitn  +  un * pitt) / t;
	precision Gxn_pi = (ut * pixn  +  un * pitx) / t;
	precision Gyn_pi = (ut * piyn  +  un * pity) / t;
	precision Gnn_pi = 2. * (ut * pinn  +  un * pitn) / t;		// fixed factor of 2 bug on 7/23

	// pi_T^{\mu\alpha} . a_\alpha
	precision piat = pitt * at  -  pitx * ax  -  pity * ay  -  t2 * pitn * an;
	precision piax = pitx * at  -  pixx * ax  -  pixy * ay  -  t2 * pixn * an;
	precision piay = pity * at  -  pixy * ax  -  piyy * ay  -  t2 * piyn * an;
	precision pian = pitn * at  -  pixn * ax  -  piyn * ay  -  t2 * pinn * an;

	// pi_T^{\mu\alpha} . \dot{z}_\alpha
	precision piDzt = pitt * D_zt  -  t2 * pitn * D_zn;
	precision piDzx = pitx * D_zt  -  t2 * pixn * D_zn;
	precision piDzy = pity * D_zt  -  t2 * piyn * D_zn;
	precision piDzn = pitn * D_zt  -  t2 * pinn * D_zn;

	// Product rule terms: P^{\mu\nu} = 2.(- u^{(\mu} . \pi_{\nu)\alpha} . a_\alpha  +  z^{(\mu} . \pi_{\nu)\alpha} . \dot{z}_\alpha)
	precision Ptt_pi = 2. * (- ut * piat  +  zt * piDzt);
	precision Ptx_pi = - ut * piax  -  ux * piat  +  zt * piDzx;
	precision Pty_pi = - ut * piay  -  uy * piat  +  zt * piDzy;
	precision Ptn_pi = - ut * pian  -  un * piat  +  zt * piDzn  +  zn * piDzt;
	precision Pxx_pi = - 2. * ux * piax;
	precision Pxy_pi = - ux * piay  -  uy * piax;
	precision Pxn_pi = - ux * pian  -  un * piax  +  zn * piDzx;	// don't see anything wrong with this
	precision Pyy_pi = - 2. * uy * piay;
	precision Pyn_pi = - uy * pian  -  un * piay  +  zn * piDzy;
	precision Pnn_pi = 2. * (- un * pian  +  zn * piDzn);

#else
	precision pi_sT = 0;
#endif


#ifdef WTZMU
	// W^munu components and derivatives
	precision Wtt = 2. * WtTz * zt;
	precision Wtx = WxTz * zt;
	precision Wty = WyTz * zt;
	precision Wtn = WtTz * zn  +  WnTz * zt;
	precision Wxn = WxTz * zn;
	precision Wyn = WyTz * zn;
	precision Wnn = 2. * WnTz * zn;

	precision dWtt_dx = 2. * (dWtTz_dx * zt  +  WtTz * dzt_dx);
	precision dWtt_dy = 2. * (dWtTz_dy * zt  +  WtTz * dzt_dy);
	precision dWtt_dn = 2. * (dWtTz_dn * zt  +  WtTz * dzt_dn);

	precision dWtx_dx = dWxTz_dx * zt  +  WxTz * dzt_dx;
	precision dWtx_dy = dWxTz_dy * zt  +  WxTz * dzt_dy;
	precision dWtx_dn = dWxTz_dn * zt  +  WxTz * dzt_dn;

	precision dWty_dx = dWyTz_dx * zt  +  WyTz * dzt_dx;
	precision dWty_dy = dWyTz_dy * zt  +  WyTz * dzt_dy;
	precision dWty_dn = dWyTz_dn * zt  +  WyTz * dzt_dn;

	precision dWtn_dx = dWtTz_dx * zn  +  WtTz * dzn_dx  +  dWnTz_dx * zt  +  WnTz * dzt_dx;
	precision dWtn_dy = dWtTz_dy * zn  +  WtTz * dzn_dy  +  dWnTz_dy * zt  +  WnTz * dzt_dy;
	precision dWtn_dn = dWtTz_dn * zn  +  WtTz * dzn_dn  +  dWnTz_dn * zt  +  WnTz * dzt_dn;

	precision dWxn_dx = dWxTz_dx * zn  +  WxTz * dzn_dx;
	precision dWxn_dn = dWxTz_dn * zn  +  WxTz * dzn_dn;

	precision dWyn_dy = dWyTz_dy * zn  +  WyTz * dzn_dy;
	precision dWyn_dn = dWyTz_dn * zn  +  WyTz * dzn_dn;

	precision dWnn_dn = 2. * (dWnTz_dn * zn  +  WnTz * dzn_dn);

	// 2nd order gradient terms in dpl and dpt
	precision WTz_D_z = WtTz * D_zt  -  t2 * WnTz * D_zn;
	precision WTz_Dz_u = WtTz * Dz_ut  -  WxTz * Dz_ux  -  WyTz * Dz_uy  -  t2 * WnTz * Dz_un;
	precision WTz_z_NabT_u = WtTz * z_NabTt_u  -  WxTz * z_NabTx_u  -  WyTz * z_NabTy_u  -  t2 * WnTz * z_NabTn_u;

	precision IplW = - 2. * WTz_D_z  +  lambda_WuL * WTz_Dz_u  +  lambda_WTL * WTz_z_NabT_u;

#if (PT_MATCHING == 1)
	precision IptW = WTz_D_z  +  lambda_WuT * WTz_Dz_u  -  lambda_WTT * WTz_z_NabT_u;
#endif

	// WTz transport coefficients
	precision eta_uW = aniso.eta_uW;
	precision eta_TW = aniso.eta_TW;
	precision tau_zW = aniso.tau_zW;
	precision delta_WW = aniso.delta_WW;
	precision lambda_WuW = aniso.lambda_WuW;
	precision lambda_WTW = aniso.lambda_WTW;
	precision lambda_piuW = aniso.lambda_piuW;
	precision lambda_piTW = aniso.lambda_piTW;

	// gradient terms in dWTz (some are unprojected)
	precision It_W1 = 2. * eta_uW * Dz_ut;
	precision Ix_W1 = 2. * eta_uW * Dz_ux;
	precision Iy_W1 = 2. * eta_uW * Dz_uy;
	precision In_W1 = 2. * eta_uW * Dz_un;

	precision It_W2 = 2. * eta_TW * z_NabTt_u;
	precision Ix_W2 = 2. * eta_TW * z_NabTx_u;
	precision Iy_W2 = 2. * eta_TW * z_NabTy_u;
	precision In_W2 = 2. * eta_TW * z_NabTn_u;

	precision It_W3 = tau_zW * D_zt;
	precision In_W3 = tau_zW * D_zn;

#ifdef PIMUNU
	precision It_W4 = pitt * D_zt  -  t2 * pitn * D_zn;
	precision In_W4 = pitn * D_zt  -  t2 * pinn * D_zn;
#else
	precision It_W4 = 0, In_W4 = 0;
#endif

	precision It_W5 = delta_WW * WtTz * thetaT;
	precision Ix_W5 = delta_WW * WxTz * thetaT;
	precision Iy_W5 = delta_WW * WyTz * thetaT;
	precision In_W5 = delta_WW * WnTz * thetaT;

	precision It_W6 = lambda_WuW * WtTz * thetaL;
	precision Ix_W6 = lambda_WuW * WxTz * thetaL;
	precision Iy_W6 = lambda_WuW * WyTz * thetaL;
	precision In_W6 = lambda_WuW * WnTz * thetaL;

	precision It_W7 = lambda_WTW * (sTtt * WtTz  -  sTtx * WxTz  -  sTty * WyTz  - t2 * sTtn * WnTz);
	precision Ix_W7 = lambda_WTW * (sTtx * WtTz  -  sTxx * WxTz  -  sTxy * WyTz  - t2 * sTxn * WnTz);
	precision Iy_W7 = lambda_WTW * (sTty * WtTz  -  sTxy * WxTz  -  sTyy * WyTz  - t2 * sTyn * WnTz);
	precision In_W7 = lambda_WTW * (sTtn * WtTz  -  sTxn * WxTz  -  sTyn * WyTz  - t2 * sTnn * WnTz);

	precision It_W8 = 0;		// keep vorticity term zero for now..
	precision Ix_W8 = 0;
	precision Iy_W8 = 0;
	precision In_W8 = 0;

#ifdef PIMUNU
	precision It_W9 = lambda_piuW * (pitt * Dz_ut  -  pitx * Dz_ux  -  pity * Dz_uy  -  t2 * pitn * Dz_un);
	precision Ix_W9 = lambda_piuW * (pitx * Dz_ut  -  pixx * Dz_ux  -  pixy * Dz_uy  -  t2 * pixn * Dz_un);
	precision Iy_W9 = lambda_piuW * (pity * Dz_ut  -  pixy * Dz_ux  -  piyy * Dz_uy  -  t2 * piyn * Dz_un);
	precision In_W9 = lambda_piuW * (pitn * Dz_ut  -  pixn * Dz_ux  -  piyn * Dz_uy  -  t2 * pinn * Dz_un);

	precision It_W10 = lambda_piTW * (pitt * z_NabTt_u  -  pitx * z_NabTx_u  -  pity * z_NabTy_u  -  t2 * pitn * z_NabTn_u);
	precision Ix_W10 = lambda_piTW * (pitx * z_NabTt_u  -  pixx * z_NabTx_u  -  pixy * z_NabTy_u  -  t2 * pixn * z_NabTn_u);
	precision Iy_W10 = lambda_piTW * (pity * z_NabTt_u  -  pixy * z_NabTx_u  -  piyy * z_NabTy_u  -  t2 * piyn * z_NabTn_u);
	precision In_W10 = lambda_piTW * (pitn * z_NabTt_u  -  pixn * z_NabTx_u  -  piyn * z_NabTy_u  -  t2 * pinn * z_NabTn_u);
#else
	precision It_W9 = 0, Ix_W9 = 0, Iy_W9 = 0, In_W9 = 0;
	precision It_W10 = 0, Ix_W10 = 0, Iy_W10 = 0, In_W10 = 0;
#endif

	precision It_W = It_W1  -  It_W2  -  It_W3  -  It_W4  +  It_W5  -  It_W6  +  It_W7  +  It_W8  +  It_W9  -  It_W10;
	precision Ix_W = Ix_W1  -  Ix_W2  			   		  +  Ix_W5  -  Ix_W6  +  Ix_W7  +  Ix_W8  +  Ix_W9  -  Ix_W10;
	precision Iy_W = Iy_W1  -  Iy_W2  			 		  +  Iy_W5  -  Iy_W6  +  Iy_W7  +  Iy_W8  +  Iy_W9  -  Iy_W10;
	precision In_W = In_W1  -  In_W2  -  In_W3  -  In_W4  +  In_W5  -  In_W6  +  In_W7  +  In_W8  +  In_W9  -  In_W10;
	Xi.transverse_project_vector(It_W, Ix_W, Iy_W, In_W);

	// Christofel terms (G_W^\mu  =  u^\alpha . \Gamma^\mu_{\alpha\beta} . WTz^\mu)
	precision Gt_W = tun * WnTz;
	precision Gn_W = (ut * WnTz  +  un * WtTz) / t;

	// product rule terms (P_W^\mu  =  - u^\mu . Wtz^\nu . a_\nu  +  z^\mu . Wtz^\nu . Dz_\nu)
	precision WTza  = WtTz * at  -  WxTz * ax  -  WyTz * ay  -  t2 * WnTz * an;
	precision WTzDz = WtTz * D_zt  -  t2 * WnTz * D_zn;

	precision Pt_W = - ut * WTza  +  zt * WTzDz;
	precision Px_W = - ux * WTza;
	precision Py_W = - uy * WTza;
	precision Pn_W = - un * WTza  +  zn * WTzDz;
#else
	precision Wtt = 0, Wtx = 0, Wty = 0, Wtn = 0, Wxn = 0, Wyn = 0, Wnn = 0;
	precision dWtt_dx = 0, dWtt_dy = 0, dWtt_dn = 0, dWtx_dx = 0, dWtx_dy = 0, dWtx_dn = 0, dWty_dx = 0, dWty_dy = 0, dWty_dn = 0, dWtn_dx = 0, dWtn_dy = 0, dWtn_dn = 0, dWxn_dx = 0,dWxn_dn = 0, dWyn_dy = 0, dWyn_dn = 0, dWnn_dn = 0;

	precision IplW = 0;

#if (PT_MATCHING == 1)
	precision IptW = 0;
#endif
#endif

	// conservation laws (check again except last line comments)
	precision tnn = (e + pt) * un2  +  pt / t2  +  Lnn  +  Wnn  +  pinn;

	S[0] =	- (ttt / t  +  t * tnn)  +  div_v * (Ltt  +  Wtt  +  pitt  -  pt)
			+  vx * (dLtt_dx  +  dWtt_dx  +  dpitt_dx  -  dpt_dx)  -  dWtx_dx  -  dpitx_dx
			+  vy * (dLtt_dy  +  dWtt_dy  +  dpitt_dy  -  dpt_dy)  -  dWty_dy  -  dpity_dy
			+  vn * (dLtt_dn  +  dWtt_dn  +  dpitt_dn  -  dpt_dn)  -  dWtn_dn  -  dpitn_dn  -  dLtn_dn;




	S[1] =	- ttx / t  -  dpt_dx  +  div_v * (Wtx  +  pitx)
			+  vx * (dWtx_dx  +  dpitx_dx)  -  dpixx_dx

			+  vy * (dWtx_dy  +  dpitx_dy)  -  dpixy_dy

			+  vn * (dWtx_dn  +  dpitx_dn)  -  dpixn_dn  -  dWxn_dn
			-  0.5 * (vx * dLtn_dn  -  Ltn * dvx_dn);	// go over this line again





	S[2] =	- tty / t  -  dpt_dy  +  div_v * (Wty  +  pity)
			+  vx * (dWty_dx  +  dpity_dx)  -  dpixy_dx
			+  vy * (dWty_dy  +  dpity_dy)  -  dpiyy_dy
			+  vn * (dWty_dn  +  dpity_dn)  -  dpiyn_dn  -  dWyn_dn
			-  0.5 * (vy * dLtn_dn  -  Ltn * dvy_dn);  // go over this line again

	S[3] =	- 3. * ttn / t  -  dpt_dn / t2  +  div_v * (Ltn  +  Wtn  +  pitn)
			+  vx * (dLtn_dx  +  dWtn_dx  +  dpitn_dx)  -  dWxn_dx  -  dpixn_dx
			+  vy * (dLtn_dy  +  dWtn_dy  +  dpitn_dy)  -  dWyn_dy  -  dpiyn_dy
			+  vn * (dLtn_dn  +  dWtn_dn  +  dpitn_dn)  -  dWnn_dn  -  dpinn_dn  -  dLnn_dn;

	// pl relaxation equation
	precision dpl = - dp * taupiInv / 1.5  +  zeta_LL * thetaL  +  zeta_TL * thetaT  +  IplW  -  lambda_piL * pi_sT;
	S[4] =	dpl / ut  +  div_v * pl;

	a = 5;		// reset counter
	// pt relaxation equation
#if (PT_MATCHING == 1)
	precision dpt =	dp * taupiInv / 3.  +  zeta_LT * thetaL  +  zeta_TT * thetaT  +  IptW  +  lambda_piT * pi_sT;
	S[a] = dpt / ut  +  div_v * pt;		a++;
#endif

	// piT relaxation equation (checked it, looks ok)
#ifdef PIMUNU
	precision dpitt = - pitt * taupiInv  +  Itt_pi  +  Ptt_pi  -  Gtt_pi;
	precision dpitx = - pitx * taupiInv  +  Itx_pi  +  Ptx_pi  -  Gtx_pi;
	precision dpity = - pity * taupiInv  +  Ity_pi  +  Pty_pi  -  Gty_pi;
	precision dpitn = - pitn * taupiInv  +  Itn_pi  +  Ptn_pi  -  Gtn_pi;
	precision dpixx = - pixx * taupiInv  +  Ixx_pi  +  Pxx_pi;
	precision dpixy = - pixy * taupiInv  +  Ixy_pi  +  Pxy_pi;
	precision dpixn = - pixn * taupiInv  +  Ixn_pi  +  Pxn_pi  -  Gxn_pi;
	precision dpiyy = - piyy * taupiInv  +  Iyy_pi  +  Pyy_pi;
	precision dpiyn = - piyn * taupiInv  +  Iyn_pi  +  Pyn_pi  -  Gyn_pi;
	precision dpinn = - pinn * taupiInv  +  Inn_pi  +  Pnn_pi  -  Gnn_pi;

	S[a] = dpitt / ut  +  div_v * pitt;		a++;
	S[a] = dpitx / ut  +  div_v * pitx;		a++;
	S[a] = dpity / ut  +  div_v * pity;		a++;
	S[a] = dpitn / ut  +  div_v * pitn;		a++;
	S[a] = dpixx / ut  +  div_v * pixx;		a++;
	S[a] = dpixy / ut  +  div_v * pixy;		a++;
	S[a] = dpixn / ut  +  div_v * pixn;		a++;
	S[a] = dpiyy / ut  +  div_v * piyy;		a++;
	S[a] = dpiyn / ut  +  div_v * piyn;		a++;
	S[a] = dpinn / ut  +  div_v * pinn;		a++;
#endif

	// WTz relaxation equation
#ifdef WTZMU
	precision dWtTz = - WtTz * taupiInv  +  It_W  +  Pt_W  -  Gt_W;
	precision dWxTz = - WxTz * taupiInv  +  Ix_W  +  Px_W;
	precision dWyTz = - WyTz * taupiInv  +  Iy_W  +  Py_W;
	precision dWnTz = - WnTz * taupiInv  +  In_W  +  Pn_W  -  Gn_W;

	S[a] = dWtTz / ut  +  div_v * WtTz;		a++;
	S[a] = dWxTz / ut  +  div_v * WxTz;		a++;
	S[a] = dWyTz / ut  +  div_v * WyTz;		a++;
	S[a] = dWnTz / ut  +  div_v * WnTz;
#endif

}
*/







void source_terms_viscous_hydro(precision * const __restrict__ S, const precision * const __restrict__ q, precision e, precision t, const precision * const __restrict__ qi1, const precision * const __restrict__ qj1, const precision * const __restrict__ qk1, const precision * const __restrict__ e1, const precision * const __restrict__ ui1, const precision * const __restrict__ uj1, const precision * const __restrict__ uk1, precision ux, precision uy, precision un, precision ux_p, precision uy_p, precision un_p, precision dt, precision dx, precision dy, precision dn, precision etabar)
{
	precision t2 = t * t;
	precision t4 = t2 * t2;
	precision tun  = t * un;

	precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  tun * tun);
	precision vx = ux / ut;
	precision vy = uy / ut;
	precision vn = un / ut;

	precision T = effectiveTemperature(e);
	precision p = equilibriumPressure(e);
	precision cs2 = speedOfSoundSquared(e);

	// conserved variables
	precision ttt = q[0];
	precision ttx = q[1];
	precision tty = q[2];
	precision ttn = q[3];

	int a = 4;	// counter
#ifdef PIMUNU
	precision pitt = q[a];	a++;
	precision pitx = q[a];	a++;
	precision pity = q[a];	a++;
	precision pitn = q[a];	a++;
	precision pixx = q[a];	a++;
	precision pixy = q[a];	a++;
	precision pixn = q[a];	a++;
	precision piyy = q[a];	a++;
	precision piyn = q[a];	a++;
	precision pinn = q[a];	a++;

	// shear transport coefficients
	precision betapi = (e + p) / 5.;
	precision taupiInv = T / (5. * etabar);
	precision delta_pipi = 4./3. * taupiInv;
	precision tau_pipi = 10./7 * taupiInv;
	precision lambda_pibulkPi = 1.2 * taupiInv;
#else
	precision pitt = 0, pitx = 0, pity = 0, pitn = 0, pixx = 0, pixy = 0, pixn = 0, piyy = 0, piyn = 0, pinn = 0;
#endif

#ifdef PI
	precision Pi = q[a];

	// bulk transport coefficients
	precision betabulk = 15. * (1./3. - cs2) * (1./3. - cs2) * (e + p);
	precision taubulkInv = 15. * (1./3. - cs2) * (1./3. - cs2) * T / zetabar(T);
	precision lambda_bulkPipi = 1.6 * (1./3. - cs2) * taubulkInv;
	precision delta_bulkPibulkPi = 2./3. * taubulkInv;
#else
	precision Pi = 0;
#endif

	// energy density and equilibrium pressure derivatives
	precision de_dx = central_derivative(e1, 0, dx);
	precision de_dy = central_derivative(e1, 2, dy);
	precision de_dn = central_derivative(e1, 4, dn);
	precision dp_dx = cs2 * de_dx;
	precision dp_dy = cs2 * de_dy;
	precision dp_dn = cs2 * de_dn;

	int n = 8;
#ifdef PIMUNU	// pi^{\mu\nu} derivatives
	precision dpitt_dx = central_derivative(qi1, n, dx);
	precision dpitt_dy = central_derivative(qj1, n, dy);
	precision dpitt_dn = central_derivative(qk1, n, dn);	n += 2;

	precision dpitx_dx = central_derivative(qi1, n, dx);
	precision dpitx_dy = central_derivative(qj1, n, dy);
	precision dpitx_dn = central_derivative(qk1, n, dn);	n += 2;

	precision dpity_dx = central_derivative(qi1, n, dx);
	precision dpity_dy = central_derivative(qj1, n, dy);
	precision dpity_dn = central_derivative(qk1, n, dn);	n += 2;

	precision dpitn_dx = central_derivative(qi1, n, dx);
	precision dpitn_dy = central_derivative(qj1, n, dy);
	precision dpitn_dn = central_derivative(qk1, n, dn);	n += 2;

	precision dpixx_dx = central_derivative(qi1, n, dx);	n += 2;

	precision dpixy_dx = central_derivative(qi1, n, dx);
	precision dpixy_dy = central_derivative(qj1, n, dy);	n += 2;

	precision dpixn_dx = central_derivative(qi1, n, dx);
	precision dpixn_dn = central_derivative(qk1, n, dn);	n += 2;

	precision dpiyy_dy = central_derivative(qj1, n, dy);	n += 2;

	precision dpiyn_dy = central_derivative(qj1, n, dy);
	precision dpiyn_dn = central_derivative(qk1, n, dn);	n += 2;

	precision dpinn_dn = central_derivative(qk1, n, dn);	n += 2;

#else
	precision dpitt_dx = 0, dpitt_dy = 0, dpitt_dn = 0, dpitx_dx = 0, dpitx_dy = 0, dpitx_dn = 0, dpity_dx = 0, dpity_dy = 0, dpity_dn = 0, dpitn_dx = 0, dpitn_dy = 0, dpitn_dn = 0, dpixx_dx = 0, dpixy_dx = 0, dpixy_dy = 0, dpixn_dx = 0, dpixn_dn = 0, dpiyy_dy = 0, dpiyn_dy = 0, dpiyn_dn = 0, dpinn_dn = 0;
#endif
#ifdef PI 	// Pi derivatives
	precision dPi_dx = central_derivative(qi1, n, dx);
	precision dPi_dy = central_derivative(qj1, n, dy);
	precision dPi_dn = central_derivative(qk1, n, dn);
#else
	precision dPi_dx = 0, dPi_dy = 0, dPi_dn = 0;
#endif

	// fluid velocity derivatives
	precision dux_dt = (ux - ux_p) / dt;
	precision dux_dx = central_derivative(ui1, 0, dx);
	precision dux_dy = central_derivative(uj1, 0, dy);
	precision dux_dn = central_derivative(uk1, 0, dn);

	precision duy_dt = (uy - uy_p) / dt;
	precision duy_dx = central_derivative(ui1, 2, dx);
	precision duy_dy = central_derivative(uj1, 2, dy);
	precision duy_dn = central_derivative(uk1, 2, dn);

	precision dun_dt = (un - un_p) / dt;
	precision dun_dx = central_derivative(ui1, 4, dx);
	precision dun_dy = central_derivative(uj1, 4, dy);
	precision dun_dn = central_derivative(uk1, 4, dn);

	precision dut_dt = vx * dux_dt  +  vy * duy_dt  +  t2 * vn * dun_dt  +  tun * vn;
	precision dut_dx = vx * dux_dx  +  vy * duy_dx  +  t2 * vn * dun_dx;
	precision dut_dy = vx * dux_dy  +  vy * duy_dy  +  t2 * vn * dun_dy;
	precision dut_dn = vx * dux_dn  +  vy * duy_dn  +  t2 * vn * dun_dn;

	// scalar expansion rate: theta = D_\mu u^\mu
	precision theta = dut_dt  +  dux_dx  +  duy_dy  +  dun_dn  +  ut / t;

	// spatial velocity derivatives and divergence
	precision dvx_dx = (dux_dx  -  vx * dut_dx) / ut;
	precision dvy_dy = (duy_dy  -  vy * dut_dy) / ut;
	precision dvn_dn = (dun_dn  -  vn * dut_dn) / ut;
	precision div_v = dvx_dx  +  dvy_dy  +  dvn_dn;

#if (NUMBER_OF_VISCOUS_CURRENTS != 0)
	spatial_projection Delta(ut, ux, uy, un, t2);		// \Delta^{\mu\nu}
	double_spatial_projection Delta_2(Delta, t2, t4);	// \Delta^{\mu\nu\alpha\beta}

#ifdef PIMUNU
	Delta_2.double_spatial_project_tensor(pitt, pitx, pity, pitn, pixx, pixy, pixn, piyy, piyn, pinn);

	// acceleration = D u^\mu
	precision at = ut * dut_dt  +  ux * dut_dx  +  uy * dut_dy  +  un * dut_dn  +  tun * un;
	precision ax = ut * dux_dt  +  ux * dux_dx  +  uy * dux_dy  +  un * dux_dn;
	precision ay = ut * duy_dt  +  ux * duy_dx  +  uy * duy_dy  +  un * duy_dn;
	precision an = ut * dun_dt  +  ux * dun_dx  +  uy * dun_dy  +  un * dun_dn  +  2. * ut * un / t;
#endif

	// shear velocity tensor = sigma^{\mu\nu} (first compute D^{(\mu) u^{\nu)})
	precision stt = dut_dt;
	precision stx = (dux_dt  -  dut_dx) / 2.;
	precision sty = (duy_dt  -  dut_dy) / 2.;
	precision stn = (dun_dt  +  un / t  -  (dut_dn  +  tun) / t2) / 2.;
	precision sxx = - dux_dx;
	precision sxy = - (dux_dy  +  duy_dx) / 2.;
	precision sxn = - (dun_dx  +  dux_dn / t2) / 2.;
	precision syy = - duy_dy;
	precision syn = - (dun_dy  +  duy_dn / t2) / 2.;
	precision snn = - (dun_dn  +  ut / t) / t2;
	Delta_2.double_spatial_project_tensor(stt, stx, sty, stn, sxx, sxy, sxn, syy, syn, snn);
#endif

#ifdef PIMUNU
	precision pi_s = pitt * stt  +  pixx * sxx  +  piyy * syy  +  t4 * pinn * snn  +  2. * (pixy * sxy  -  pitx * stx  -  pity * sty  +  t2 * (pixn * sxn  +  piyn * syn  -  pitn * stn));

	// pi gradient terms (tensors are symmetric but (2,4) not double spatial projected until end)

	// 2 . (eta/tau_pi) . \sigma_^{\mu\nu}
	precision Itt_pi1 = 2. * betapi * stt;
	precision Itx_pi1 = 2. * betapi * stx;
	precision Ity_pi1 = 2. * betapi * sty;
	precision Itn_pi1 = 2. * betapi * stn;
	precision Ixx_pi1 = 2. * betapi * sxx;
	precision Ixy_pi1 = 2. * betapi * sxy;
	precision Ixn_pi1 = 2. * betapi * sxn;
	precision Iyy_pi1 = 2. * betapi * syy;
	precision Iyn_pi1 = 2. * betapi * syn;
	precision Inn_pi1 = 2. * betapi * snn;

	// ignore the vorticity term for now
	precision Itt_pi2 = 0;
	precision Itx_pi2 = 0;
	precision Ity_pi2 = 0;
	precision Itn_pi2 = 0;
	precision Ixx_pi2 = 0;
	precision Ixy_pi2 = 0;
	precision Ixn_pi2 = 0;
	precision Iyy_pi2 = 0;
	precision Iyn_pi2 = 0;
	precision Inn_pi2 = 0;

	// \delta_{\pi\pi} . pi^{\mu\nu} . \theta
	precision Itt_pi3 = delta_pipi * pitt * theta;
	precision Itx_pi3 = delta_pipi * pitx * theta;
	precision Ity_pi3 = delta_pipi * pity * theta;
	precision Itn_pi3 = delta_pipi * pitn * theta;
	precision Ixx_pi3 = delta_pipi * pixx * theta;
	precision Ixy_pi3 = delta_pipi * pixy * theta;
	precision Ixn_pi3 = delta_pipi * pixn * theta;
	precision Iyy_pi3 = delta_pipi * piyy * theta;
	precision Iyn_pi3 = delta_pipi * piyn * theta;
	precision Inn_pi3 = delta_pipi * pinn * theta;

	// \tau_{\pi\pi} . \pi^{\lambda (\mu} . \sigma_T^{\nu)}_\lambda
	precision Itt_pi4 = tau_pipi * (pitt * stt  -  pitx * stx  -  pity * sty  -  t2 * pitn * stn);
	precision Itx_pi4 = tau_pipi * (pitt * stx  +  pitx * stt  -  pitx * sxx  -  pixx * stx  -  pity * sxy  -  pixy * sty  -  t2 * (pitn * sxn  +  pixn * stn)) / 2.;
	precision Ity_pi4 = tau_pipi * (pitt * sty  +  pity * stt  -  pitx * sxy  -  pixy * stx  -  pity * syy  -  piyy * sty  -  t2 * (pitn * syn  +  piyn * stn)) / 2.;
	precision Itn_pi4 = tau_pipi * (pitt * stn  +  pitn * stt  -  pitx * sxn  -  pixn * stx  -  pity * syn  -  piyn * sty  -  t2 * (pitn * snn  +  pinn * stn)) / 2.;
	precision Ixx_pi4 = tau_pipi * (pitx * stx  -  pixx * sxx  -  pixy * sxy  -  t2 * pixn * sxn);
	precision Ixy_pi4 = tau_pipi * (pitx * sty  +  pity * stx  -  pixx * sxy  -  pixy * sxx  -  pixy * syy  -  piyy * sxy  -  t2 * (pixn * syn  +  piyn * sxn)) / 2.;
	precision Ixn_pi4 = tau_pipi * (pitx * stn  +  pitn * stx  -  pixx * sxn  -  pixn * sxx  -  pixy * syn  -  piyn * sxy  -  t2 * (pixn * snn  +  pinn * sxn)) / 2.;
	precision Iyy_pi4 = tau_pipi * (pity * sty  -  pixy * sxy  -  piyy * syy  -  t2 * piyn * syn);
	precision Iyn_pi4 = tau_pipi * (pity * stn  +  pitn * sty  -  pixy * sxn  -  pixn * sxy  -  piyy * syn  -  piyn * syy  -  t2 * (piyn * snn  +  pinn * syn)) / 2.;
	precision Inn_pi4 = tau_pipi * (pitn * stn  -  pixn * sxn  -  piyn * syn  -  t2 * pinn * snn);

	// \lambda^{\pi\Pi} . \sigma^{\mu\nu} . \Pi
	precision Itt_pi5 = lambda_pibulkPi * Pi * stt;
	precision Itx_pi5 = lambda_pibulkPi * Pi * stx;
	precision Ity_pi5 = lambda_pibulkPi * Pi * sty;
	precision Itn_pi5 = lambda_pibulkPi * Pi * stn;
	precision Ixx_pi5 = lambda_pibulkPi * Pi * sxx;
	precision Ixy_pi5 = lambda_pibulkPi * Pi * sxy;
	precision Ixn_pi5 = lambda_pibulkPi * Pi * sxn;
	precision Iyy_pi5 = lambda_pibulkPi * Pi * syy;
	precision Iyn_pi5 = lambda_pibulkPi * Pi * syn;
	precision Inn_pi5 = lambda_pibulkPi * Pi * snn;

	precision Itt_pi = Itt_pi1  +  Itt_pi2  -  Itt_pi3  -  Itt_pi4  +  Itt_pi5;
	precision Itx_pi = Itx_pi1  +  Itx_pi2  -  Itx_pi3  -  Itx_pi4  +  Itx_pi5;
	precision Ity_pi = Ity_pi1  +  Ity_pi2  -  Ity_pi3  -  Ity_pi4  +  Ity_pi5;
	precision Itn_pi = Itn_pi1  +  Itn_pi2  -  Itn_pi3  -  Itn_pi4  +  Itn_pi5;
	precision Ixx_pi = Ixx_pi1  +  Ixx_pi2  -  Ixx_pi3  -  Ixx_pi4  +  Ixx_pi5;
	precision Ixy_pi = Ixy_pi1  +  Ixy_pi2	-  Ixy_pi3  -  Ixy_pi4  +  Ixy_pi5;
	precision Ixn_pi = Ixn_pi1  +  Ixn_pi2  -  Ixn_pi3  -  Ixn_pi4  +  Ixn_pi5;
	precision Iyy_pi = Iyy_pi1  +  Iyy_pi2  -  Iyy_pi3  -  Iyy_pi4  +  Iyy_pi5;
	precision Iyn_pi = Iyn_pi1  +  Iyn_pi2  -  Iyn_pi3  -  Iyn_pi4  +  Iyn_pi5;
	precision Inn_pi = Inn_pi1  +  Inn_pi2  -  Inn_pi3  -  Inn_pi4  +  Inn_pi5;
	Delta_2.double_spatial_project_tensor(Itt_pi, Itx_pi, Ity_pi, Itn_pi, Ixx_pi, Ixy_pi, Ixn_pi, Iyy_pi, Iyn_pi, Inn_pi);

	// Christofel terms: G_\pi^{\mu\nu} = 2 . u^\alpha . \Gamma^{(\mu}_{\alpha\beta} . \pi^{\beta\nu)}
	precision Gtt_pi = 2. * tun * pitn;
	precision Gtx_pi = tun * pixn;
	precision Gty_pi = tun * piyn;
	precision Gtn_pi = tun * pinn  +  (ut * pitn  +  un * pitt) / t;
	precision Gxn_pi = (ut * pixn  +  un * pitx) / t;
	precision Gyn_pi = (ut * piyn  +  un * pity) / t;
	precision Gnn_pi = 2. * (ut * pinn  +  un * pitn) / t;

	// pi^{\mu\alpha} . a_\alpha
	precision piat = pitt * at  -  pitx * ax  -  pity * ay  -  t2 * pitn * an;
	precision piax = pitx * at  -  pixx * ax  -  pixy * ay  -  t2 * pixn * an;
	precision piay = pity * at  -  pixy * ax  -  piyy * ay  -  t2 * piyn * an;
	precision pian = pitn * at  -  pixn * ax  -  piyn * ay  -  t2 * pinn * an;

	// product rule terms: P^{\mu\nu} = 2 . u^{(\mu} . \pi_{\nu)\lambda} . a_\lambda
	precision Ptt_pi = 2. * ut * piat;
	precision Ptx_pi = ut * piax  -  ux * piat;
	precision Pty_pi = ut * piay  -  uy * piat;
	precision Ptn_pi = ut * pian  -  un * piat;
	precision Pxx_pi = 2. * ux * piax;
	precision Pxy_pi = ux * piay  -  uy * piax;
	precision Pxn_pi = ux * pian  -  un * piax;
	precision Pyy_pi = 2. * uy * piay;
	precision Pyn_pi = uy * pian  -  un * piay;
	precision Pnn_pi = 2. * un * pian;
#else
	precision pi_s = 0;
#endif

	// conservation laws
	precision tnn = (e + p) * un * un  +  (p + Pi) / t2  +  pinn;

	S[0] =	- (ttt / t  +  t * tnn)  +  div_v * (pitt  -  p  -  Pi)
			+  vx * (dpitt_dx  -  dp_dx  -  dPi_dx)  -  dpitx_dx
			+  vy * (dpitt_dy  -  dp_dy  -  dPi_dy)  -  dpity_dy
			+  vn * (dpitt_dn  -  dp_dn  -  dPi_dn)  -  dpitn_dn;

	S[1] =	- ttx / t  -  dp_dx  -  dPi_dx  +  div_v * pitx
			+  vx * dpitx_dx  -  dpixx_dx
			+  vy * dpitx_dy  -  dpixy_dy
			+  vn * dpitx_dn  -  dpixn_dn;

	S[2] =	- tty / t  -  dp_dy  -  dPi_dy  +  div_v * pity
			+  vx * dpity_dx  -  dpixy_dx
			+  vy * dpity_dy  -  dpiyy_dy
			+  vn * dpity_dn  -  dpiyn_dn;

	S[3] =	- 3. * ttn / t  -  (dp_dn + dPi_dn) / t2  +  div_v * pitn
			+  vx * dpitn_dx  -  dpixn_dx
			+  vy * dpitn_dy  -  dpiyn_dy
			+  vn * dpitn_dn  -  dpinn_dn;

	a = 4;	// reset counter

	// pimunu relaxation equation
#ifdef PIMUNU
	precision dpitt = - pitt * taupiInv  +  Itt_pi  -  Ptt_pi  -  Gtt_pi;
	precision dpitx = - pitx * taupiInv  +  Itx_pi  -  Ptx_pi  -  Gtx_pi;
	precision dpity = - pity * taupiInv  +  Ity_pi  -  Pty_pi  -  Gty_pi;
	precision dpitn = - pitn * taupiInv  +  Itn_pi  -  Ptn_pi  -  Gtn_pi;
	precision dpixx = - pixx * taupiInv  +  Ixx_pi  -  Pxx_pi;
	precision dpixy = - pixy * taupiInv  +  Ixy_pi  -  Pxy_pi;
	precision dpixn = - pixn * taupiInv  +  Ixn_pi  -  Pxn_pi  -  Gxn_pi;
	precision dpiyy = - piyy * taupiInv  +  Iyy_pi  -  Pyy_pi;
	precision dpiyn = - piyn * taupiInv  +  Iyn_pi  -  Pyn_pi  -  Gyn_pi;
	precision dpinn = - pinn * taupiInv  +  Inn_pi  -  Pnn_pi  -  Gnn_pi;

	S[a] = dpitt / ut  +  div_v * pitt;		a++;
	S[a] = dpitx / ut  +  div_v * pitx;		a++;
	S[a] = dpity / ut  +  div_v * pity;		a++;
	S[a] = dpitn / ut  +  div_v * pitn;		a++;
	S[a] = dpixx / ut  +  div_v * pixx;		a++;
	S[a] = dpixy / ut  +  div_v * pixy;		a++;
	S[a] = dpixn / ut  +  div_v * pixn;		a++;
	S[a] = dpiyy / ut  +  div_v * piyy;		a++;
	S[a] = dpiyn / ut  +  div_v * piyn;		a++;
	S[a] = dpinn / ut  +  div_v * pinn;		a++;
#endif

	// bulkPi relaxation equation
#ifdef PI
	precision dPi = - Pi * taubulkInv  -  betabulk * theta  -  delta_bulkPibulkPi * Pi * theta  +  lambda_bulkPipi * pi_s;

	S[a] = dPi / ut  +  div_v * Pi;
#endif

}








