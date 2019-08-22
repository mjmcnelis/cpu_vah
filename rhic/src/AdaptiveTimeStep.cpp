#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../include/Precision.h"
#include "../include/FluxTerms.h"
#include "../include/DynamicalVariables.h"
#include "../include/SourceTerms.h"
#include "../include/NeighborCells.h"
#include "../include/EquationOfState.h"
#include "../include/TransportCoefficients.h"

const precision fraction = 0.1;

// set the adaptive time step
precision set_time_step(hydro_time_scales dt_hydro, precision dt_min)
{
	precision dt_CFL 	= dt_hydro.dt_CFL;
	precision dt_micro 	= dt_hydro.dt_micro;
	precision dt_rate	= dt_hydro.dt_rate;

	printf("%lf\t%lf\t%lf\n", dt_CFL, dt_micro, dt_rate);

	precision dt = fmin(dt_CFL, fraction * fmin(dt_micro, dt_rate));

	if(dt < dt_min)
	{
		printf("set_time_step error: dt = %.6f < %.3f\n", dt, dt_min);
	}

	return dt_min * fmax(1., floor(dt / dt_min));	// round down dt to numerical precision
}


inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}

inline precision central_derivative(const precision * const __restrict__ f, int n, precision dx)
{
	return (f[n + 1] - f[n]) / (2. * dx);		// f[n] = fm  |	 f[n+1] = fp
}


precision compute_CFL_time_scale(const precision * const __restrict__ vxi, const precision * const __restrict__ vyj, const precision * const __restrict__ vnk, precision ut, precision ux, precision uy, precision un, precision dx, precision dy, precision dn)
{
	precision vx = ux / ut;
	precision vy = uy / ut;
	precision vn = un / ut;

	precision ax = compute_max_local_propagation_speed(vxi, vx);
	precision ay = compute_max_local_propagation_speed(vyj, vy);
	precision an = compute_max_local_propagation_speed(vnk, vn);

	return 0.125 / fmax(fabs(ax / dx), fmax(fabs(ay / dy), fabs(an / dn)));
	//return  0.125 / fmax(fabs(vx / dx), fmax(fabs(vy / dy), fabs(vn / dn)));
}


precision compute_microscopic_time_scale(precision e, precision etabar_const)
{
	precision T = effectiveTemperature(e);
	precision etabar = eta_over_s(T, etabar_const);

#ifdef CONFORMAL_EOS
	precision tau_pi = 5. * etabar / T;
	return tau_pi;
#else
	precision tau_pi = 1.;		// fill in quasiparticle later
	precision tau_bulk = 1.;	// is the bulk relaxation time ever incredibly small?
	return fmin(tau_pi, tau_bulk);
#endif
}


precision compute_comoving_evolution_rate_aniso_hydro(hydro_variables q, precision e, precision t, const precision * const __restrict__ ui1, const precision * const __restrict__ uj1, const precision * const __restrict__ uk1, precision ut, precision ux, precision uy, precision un, precision ux_p, precision uy_p, precision un_p, precision dt_prev, precision dx, precision dy, precision dn, precision etabar_const, precision dt_min)
{
	// useful expressions
	precision t2 = t * t;
	precision t4 = t2 * t2;

	// spatial velocity components
	precision vx = ux / ut;
	precision vy = uy / ut;
	precision vn = un / ut;

	// longitudinal basis vector
	precision utperp  = sqrt(1.  +  ux * ux  +  uy * uy);
	precision zt = t * un / utperp;
	precision zn = ut / t / utperp;

	// other useful expressions
	precision zt2  = zt * zt;
	precision ztzn = zt * zn;
	precision zn2  = zn * zn;
	precision t2zn = t2 * zn;
	precision ut2  = ut * ut;
	precision un2  = un * un;
	precision utun = ut * un;
	precision tun  = t  * un;
	precision unzn = un * zn;
	precision utperp2 = utperp * utperp;

	// hydro variables
	precision p = equilibriumPressure(e);
	precision cs2 = speedOfSoundSquared(e);
	precision T = effectiveTemperature(e);
	precision etabar = eta_over_s(T, etabar_const);

#ifdef CONFORMAL_EOS
	precision taupiInv = T / (5. * etabar);
#else
	precision taupiInv = 1.;	// fill in quasiparticle later
	precision taubulkInv = 1.;
#endif

	precision pl  = q.pl;

#if (PT_MATCHING == 1)
	precision pt  = q.pt;
#else
	precision pt  = (e - pl) / 2.;
#endif
#ifdef PIMUNU
	precision pitt = q.pitt;
	precision pitx = q.pitx;
	precision pity = q.pity;
	precision pitn = q.pitn;
	precision pixx = q.pixx;
	precision pixy = q.pixy;
	precision pixn = q.pixn;
	precision piyy = q.piyn;
	precision piyn = q.piyn;
	precision pinn = q.pinn;
#endif
#ifdef WTZMU
	precision WtTz = q.WtTz;
	precision WxTz = q.WxTz;
	precision WyTz = q.WyTz;
	precision WnTz = q.WnTz;
#endif

	// anisotropic transport coefficients
	transport_coefficients aniso;
	aniso.compute_transport_coefficients(e, pl, pt);

	// pl coefficients
	precision zeta_LL = aniso.zeta_LL;
	precision zeta_TL = aniso.zeta_TL;
	precision lambda_WuL = aniso.lambda_WuL;
	precision lambda_WTL = aniso.lambda_WTL;
	precision lambda_piL = aniso.lambda_piL;

	// pt coefficients
#if (PT_MATCHING == 1)
	precision zeta_LT = aniso.zeta_LT;
	precision zeta_TT = aniso.zeta_TT;
	precision lambda_WuT = aniso.lambda_WuT;
	precision lambda_WTT = aniso.lambda_WTT;
	precision lambda_piT = aniso.lambda_piT;
#endif

	// fluid velocity derivatives
	precision dux_dt = (ux - ux_p) / dt_prev;
	precision dux_dx = central_derivative(ui1, 0, dx);
	precision dux_dy = central_derivative(uj1, 0, dy);
	precision dux_dn = central_derivative(uk1, 0, dn);

	precision duy_dt = (uy - uy_p) / dt_prev;
	precision duy_dx = central_derivative(ui1, 2, dx);
	precision duy_dy = central_derivative(uj1, 2, dy);
	precision duy_dn = central_derivative(uk1, 2, dn);

	precision dun_dt = (un - un_p) / dt_prev;
	precision dun_dx = central_derivative(ui1, 4, dx);
	precision dun_dy = central_derivative(uj1, 4, dy);
	precision dun_dn = central_derivative(uk1, 4, dn);

	// chain rule for ut derivatives
	precision dut_dt = vx * dux_dt  +  vy * duy_dt  +  t2 * vn * dun_dt  +  tun * vn;
	precision dut_dx = vx * dux_dx  +  vy * duy_dx  +  t2 * vn * dun_dx;
	precision dut_dy = vx * dux_dy  +  vy * duy_dy  +  t2 * vn * dun_dy;
	precision dut_dn = vx * dux_dn  +  vy * duy_dn  +  t2 * vn * dun_dn;

	// scalar, longitudinal and transverse expansion rates: theta = D_\mu u^\mu, thetaL = z_\mu Dz u^\mu, thetaT = NablaT_\mu u^\mu
	precision theta = dut_dt  +  dux_dx  +  duy_dy  +  dun_dn  +  ut / t;
	precision thetaL = - zt2 * dut_dt  +  t2 * zn2 * dun_dn  +  ztzn * (t2 * dun_dt  -  dut_dn)  +  t * zn2 * ut;
	precision thetaT = theta  -  thetaL;

	// longitudinal vector derivatives
	precision dzt_dt = t * (dun_dt  -  un * (ux * dux_dt  +  uy * duy_dt) / utperp2) / utperp  +  zt / t;
	precision dzt_dx = t * (dun_dx  -  un * (ux * dux_dx  +  uy * duy_dx) / utperp2) / utperp;
	precision dzt_dy = t * (dun_dy  -  un * (ux * dux_dy  +  uy * duy_dy) / utperp2) / utperp;
	precision dzt_dn = t * (dun_dn  -  un * (ux * dux_dn  +  uy * duy_dn) / utperp2) / utperp;

	precision dzn_dt = (dut_dt  -  ut * (ux * dux_dt  +  uy * duy_dt) / utperp2) / (t * utperp)  -  zn / t;
	precision dzn_dx = (dut_dx  -  ut * (ux * dux_dx  +  uy * duy_dx) / utperp2) / (t * utperp);
	precision dzn_dy = (dut_dy  -  ut * (ux * dux_dy  +  uy * duy_dy) / utperp2) / (t * utperp);
	precision dzn_dn = (dut_dn  -  ut * (ux * dux_dn  +  uy * duy_dn) / utperp2) / (t * utperp);

#if (NUMBER_OF_RESIDUAL_CURRENTS != 0)
	transverse_projection Xi(ut, ux, uy, un, zt, zn, t2);	// Xi^{\mu\nu}
	double_transverse_projection Xi_2(Xi, t2, t4);			// Xi^{\mu\nu\alpha\beta}

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

	// piT source terms
#ifdef PIMUNU
	// piT transport coefficients
	precision eta_T = aniso.eta_T;
	precision delta_pipi = aniso.delta_pipi;
	precision tau_pipi = aniso.tau_pipi;
	precision lambda_pipi = aniso.lambda_pipi;
	precision lambda_Wupi = aniso.lambda_Wupi;
	precision lambda_WTpi = aniso.lambda_WTpi;

	precision pi_sT = pitt * sTtt  +  pixx * sTxx  +  piyy * sTyy  +  t4 * pinn * sTnn  +  2. * (pixy * sTxy  -  pitx * sTtx  -  pity * sTty  +  t2 * (pixn * sTxn  +  piyn * sTyn  -  pitn * sTtn));

	// gradient source terms
	// 2 . \eta_T . \sigma_T^{\mu\nu}
	precision Itt = 2. * eta_T * sTtt;
	precision Itx = 2. * eta_T * sTtx;
	precision Ity = 2. * eta_T * sTty;
	precision Itn = 2. * eta_T * sTtn;
	precision Ixx = 2. * eta_T * sTxx;
	precision Ixy = 2. * eta_T * sTxy;
	precision Ixn = 2. * eta_T * sTxn;
	precision Iyy = 2. * eta_T * sTyy;
	precision Iyn = 2. * eta_T * sTyn;
	precision Inn = 2. * eta_T * sTnn;

	// 2 . WTz^{(\mu} . \dot{z}^{\nu)}
#ifdef WTZMU
	precision Itt -= 2. * WtTz * D_zt;
	precision Itx -= WxTz * D_zt;
	precision Ity -= WyTz * D_zt;
	precision Itn -= WtTz * D_zn  +  WnTz * D_zt;
	precision Ixn -= WxTz * D_zn;
	precision Iyn -= WyTz * D_zn;
	precision Inn -= 2. * WnTz * D_zn;
#endif

	// \delta^\pi_\pi . pi_T^{\mu\nu} . \theta_T
	precision Itt -= delta_pipi * pitt * thetaT;
	precision Itx -= delta_pipi * pitx * thetaT;
	precision Ity -= delta_pipi * pity * thetaT;
	precision Itn -= delta_pipi * pitn * thetaT;
	precision Ixx -= delta_pipi * pixx * thetaT;
	precision Ixy -= delta_pipi * pixy * thetaT;
	precision Ixn -= delta_pipi * pixn * thetaT;
	precision Iyy -= delta_pipi * piyy * thetaT;
	precision Iyn -= delta_pipi * piyn * thetaT;
	precision Inn -= delta_pipi * pinn * thetaT;

	// \tau^\pi_\pi . \pi_T^{\alpha (\mu} . \sigma_T^{\nu)}_\alpha
	precision Itt -= tau_pipi * (pitt * sTtt  -  pitx * sTtx  -  pity * sTty  -  t2 * pitn * sTtn);
	precision Itx -= tau_pipi * (pitt * sTtx  +  pitx * sTtt  -  pitx * sTxx  -  pixx * sTtx  -  pity * sTxy  -  pixy * sTty  -  t2 * (pitn * sTxn  +  pixn * sTtn)) / 2.;
	precision Ity -= tau_pipi * (pitt * sTty  +  pity * sTtt  -  pitx * sTxy  -  pixy * sTtx  -  pity * sTyy  -  piyy * sTty  -  t2 * (pitn * sTyn  +  piyn * sTtn)) / 2.;
	precision Itn -= tau_pipi * (pitt * sTtn  +  pitn * sTtt  -  pitx * sTxn  -  pixn * sTtx  -  pity * sTyn  -  piyn * sTty  -  t2 * (pitn * sTnn  +  pinn * sTtn)) / 2.;
	precision Ixx -= tau_pipi * (pitx * sTtx  -  pixx * sTxx  -  pixy * sTxy  -  t2 * pixn * sTxn);
	precision Ixy -= tau_pipi * (pitx * sTty  +  pity * sTtx  -  pixx * sTxy  -  pixy * sTxx  -  pixy * sTyy  -  piyy * sTxy  -  t2 * (pixn * sTyn  +  piyn * sTxn)) / 2.;
	precision Ixn -= tau_pipi * (pitx * sTtn  +  pitn * sTtx  -  pixx * sTxn  -  pixn * sTxx  -  pixy * sTyn  -  piyn * sTxy  -  t2 * (pixn * sTnn  +  pinn * sTxn)) / 2.;
	precision Iyy -= tau_pipi * (pity * sTty  -  pixy * sTxy  -  piyy * sTyy  -  t2 * piyn * sTyn);
	precision Iyn -= tau_pipi * (pity * sTtn  +  pitn * sTty  -  pixy * sTxn  -  pixn * sTxy  -  piyy * sTyn  -  piyn * sTyy  -  t2 * (piyn * sTnn  +  pinn * sTyn)) / 2.;
	precision Inn -= tau_pipi * (pitn * sTtn  -  pixn * sTxn  -  piyn * sTyn  -  t2 * pinn * sTnn);

	// fill in vorticity source terms later
#ifdef VORTICITY
	precision Itt += 0;
	precision Itx += 0;
	precision Ity += 0;
	precision Itn += 0;
	precision Ixx += 0;
	precision Ixy += 0;
	precision Ixn += 0;
	precision Iyy += 0;
	precision Iyn += 0;
	precision Inn += 0;
#endif

	// \lambda^\pi_\pi . \pi_T^{\mu\nu} . \theta_L
	precision Itt += lambda_pipi * pitt * thetaL;
	precision Itx += lambda_pipi * pitx * thetaL;
	precision Ity += lambda_pipi * pity * thetaL;
	precision Itn += lambda_pipi * pitn * thetaL;
	precision Ixx += lambda_pipi * pixx * thetaL;
	precision Ixy += lambda_pipi * pixy * thetaL;
	precision Ixn += lambda_pipi * pixn * thetaL;
	precision Iyy += lambda_pipi * piyy * thetaL;
	precision Iyn += lambda_pipi * piyn * thetaL;
	precision Inn += lambda_pipi * pinn * thetaL;

#ifdef WTZMU
	// \lambda_Wu^\pi . WTz^{(\mu} . Dz u^{\nu)}
	precision Itt -= lambda_Wupi * (WtTz * Dz_ut);
	precision Itx -= lambda_Wupi * (WtTz * Dz_ux  +  WxTz * Dz_ut) / 2.;
	precision Ity -= lambda_Wupi * (WtTz * Dz_uy  +  WyTz * Dz_ut) / 2.;
	precision Itn -= lambda_Wupi * (WtTz * Dz_un  +  WnTz * Dz_ut) / 2.;
	precision Ixx -= lambda_Wupi * (WxTz * Dz_ux);
	precision Ixy -= lambda_Wupi * (WxTz * Dz_uy  +  WyTz * Dz_ux) / 2.;
	precision Ixn -= lambda_Wupi * (WxTz * Dz_un  +  WnTz * Dz_ux) / 2.;
	precision Iyy -= lambda_Wupi * (WyTz * Dz_uy);
	precision Iyn -= lambda_Wupi * (WyTz * Dz_un  +  WnTz * Dz_uy) / 2.;
	precision Inn -= lambda_Wupi * (WnTz * Dz_un);

	// \lambda_WT^\pi . WTz^{(\mu} . z_\alpha . \Nabla_T^{\nu)} . u^\alpha
	precision Itt += lambda_WTpi * (WtTz * z_NabTt_u);
	precision Itx += lambda_WTpi * (WtTz * z_NabTx_u  +  WxTz * z_NabTt_u) / 2.;
	precision Ity += lambda_WTpi * (WtTz * z_NabTy_u  +  WyTz * z_NabTt_u) / 2.;
	precision Itn += lambda_WTpi * (WtTz * z_NabTn_u  +  WnTz * z_NabTt_u) / 2.;
	precision Ixx += lambda_WTpi * (WxTz * Dz_ux);
	precision Ixy += lambda_WTpi * (WxTz * z_NabTy_u  +  WyTz * z_NabTx_u) / 2.;
	precision Ixn += lambda_WTpi * (WxTz * z_NabTn_u  +  WnTz * z_NabTx_u) / 2.;
	precision Iyy += lambda_WTpi * (WyTz * z_NabTy_u);
	precision Iyn += lambda_WTpi * (WyTz * z_NabTn_u  +  WnTz * z_NabTy_u) / 2.;
	precision Inn += lambda_WTpi * (WnTz * z_NabTn_u);
#endif

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
	precision It = 2. * eta_uW * Dz_ut;
	precision Ix = 2. * eta_uW * Dz_ux;
	precision Iy = 2. * eta_uW * Dz_uy;
	precision In = 2. * eta_uW * Dz_un;

	precision It -= 2. * eta_TW * z_NabTt_u;
	precision Ix -= 2. * eta_TW * z_NabTx_u;
	precision Iy -= 2. * eta_TW * z_NabTy_u;
	precision In -= 2. * eta_TW * z_NabTn_u;

	precision It -= tau_zW * D_zt;
	precision In -= tau_zW * D_zn;

#ifdef PIMUNU
	precision It -= pitt * D_zt  -  t2 * pitn * D_zn;
	precision In -= pitn * D_zt  -  t2 * pinn * D_zn;
#endif

	precision It += delta_WW * WtTz * thetaT;
	precision Ix += delta_WW * WxTz * thetaT;
	precision Iy += delta_WW * WyTz * thetaT;
	precision In += delta_WW * WnTz * thetaT;

	precision It -= lambda_WuW * WtTz * thetaL;
	precision Ix -= lambda_WuW * WxTz * thetaL;
	precision Iy -= lambda_WuW * WyTz * thetaL;
	precision In -= lambda_WuW * WnTz * thetaL;

	precision It += lambda_WTW * (sTtt * WtTz  -  sTtx * WxTz  -  sTty * WyTz  - t2 * sTtn * WnTz);
	precision Ix += lambda_WTW * (sTtx * WtTz  -  sTxx * WxTz  -  sTxy * WyTz  - t2 * sTxn * WnTz);
	precision Iy += lambda_WTW * (sTty * WtTz  -  sTxy * WxTz  -  sTyy * WyTz  - t2 * sTyn * WnTz);
	precision In += lambda_WTW * (sTtn * WtTz  -  sTxn * WxTz  -  sTyn * WyTz  - t2 * sTnn * WnTz);

	// todo: write up vorticity terms
#ifdef VORTICITY
	precision It += 0;
	precision Ix += 0;
	precision Iy += 0;
	precision In += 0;
#endif

#ifdef PIMUNU
	precision It += lambda_piuW * (pitt * Dz_ut  -  pitx * Dz_ux  -  pity * Dz_uy  -  t2 * pitn * Dz_un);
	precision Ix += lambda_piuW * (pitx * Dz_ut  -  pixx * Dz_ux  -  pixy * Dz_uy  -  t2 * pixn * Dz_un);
	precision Iy += lambda_piuW * (pity * Dz_ut  -  pixy * Dz_ux  -  piyy * Dz_uy  -  t2 * piyn * Dz_un);
	precision In += lambda_piuW * (pitn * Dz_ut  -  pixn * Dz_ux  -  piyn * Dz_uy  -  t2 * pinn * Dz_un);

	precision It -= lambda_piTW * (pitt * z_NabTt_u  -  pitx * z_NabTx_u  -  pity * z_NabTy_u  -  t2 * pitn * z_NabTn_u);
	precision Ix -= lambda_piTW * (pitx * z_NabTt_u  -  pixx * z_NabTx_u  -  pixy * z_NabTy_u  -  t2 * pixn * z_NabTn_u);
	precision Iy -= lambda_piTW * (pity * z_NabTt_u  -  pixy * z_NabTx_u  -  piyy * z_NabTy_u  -  t2 * piyn * z_NabTn_u);
	precision In -= lambda_piTW * (pitn * z_NabTt_u  -  pixn * z_NabTx_u  -  piyn * z_NabTy_u  -  t2 * pinn * z_NabTn_u);
#endif
	Xi.transverse_project_vector(It_W, Ix_W, Iy_W, In_W);

	// christofel terms (G_W^\mu  =  u^\alpha . \Gamma^\mu_{\alpha\beta} . WTz^\mu)
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
	precision WTz_Dz_u = 0;
	precision WTz_z_NabT_u = 0;

	precision IplW = 0;

#if (PT_MATCHING == 1)
	precision IptW = 0;
#endif
#endif

	// compute comoving derivatives of "scalar" hydrodynamic quantities
	precision dt_rate;

	// energy density
	precision de = - (e + pt) * thetaT  -  (e + pl) * thetaL  -  WTz_Dz_u  +  WTz_z_NabT_u  +  pi_sT;

	dt_rate = fabs(e / de);

	// fluid velocity (temporary approximation)
	precision dut = ut * dut_dt  +  ux * dut_dx  +  uy * dut_dy  +  un * dut_dn  +  t * un2;
	dt_rate = fmin(dt_rate, fabs(ut / dut));

	// viscous part of longitudinal pressure
	precision dpl = - 0*(pl - pt) * taupiInv / 1.5  +  zeta_LL * thetaL  +  zeta_TL * thetaT  +  IplW  -  lambda_piL * pi_sT;

	precision dt_pl = fabs((pl - 0*p) / (dpl  -  cs2 * de));

	dt_rate = fmin(dt_rate, dt_pl);
	//if(fabs((pl - p) / p) > 0.1)	// get rid of bad cases
	// if(e > E_MIN && pl > P_MIN)
	// {
	// 	dt_rate = fmin(dt_rate, dt_pl);

	// 	if(fraction * dt_pl < dt_min)		// I don't understand why they don't add up
	// 	{
	// 		//printf("dt_pl error at %lf fm/c: dt_pl = %.3g / |%.6g - %.6g| = %.3g\n", t, fabs(pl - p), dpl, cs2 * de, dt_pl);
	// 		//printf("(e + pt)/3 = %.6g\nzeta_TL = %.6g\n", -(e + pt)/3., zeta_TL);
	// 		printf("(e + pl)/3 = %.6g\nzeta_LL = %.6g\n", -(e + pl)/3., zeta_LL);
	// 		exit(-1);
	// 	}
	// }



	// if(fraction * dt_pl > dt_min)	// get rid of bad cases
	// {
	// 	dt_rate = fmin(dt_rate, dt_pl);
	// }



#if (PT_MATCHING == 1)
	// viscous part of transverse pressure
	precision dpt = zeta_LT * thetaL  +  zeta_TT * thetaT  +  IptW  +  lambda_piT * pi_sT;
	dt_rate = fmin(dt_rate, fabs((pt - p) / (dpt - dp)));
#endif
#ifdef PIMUNU
	precision dpitt = Itt_pi  +  Ptt_pi  -  Gtt_pi;
	precision dpitx = Itx_pi  +  Ptx_pi  -  Gtx_pi;
	precision dpity = Ity_pi  +  Pty_pi  -  Gty_pi;
	precision dpitn = Itn_pi  +  Ptn_pi  -  Gtn_pi;
	precision dpixx = Ixx_pi  +  Pxx_pi;
	precision dpixy = Ixy_pi  +  Pxy_pi;
	precision dpixn = Ixn_pi  +  Pxn_pi  -  Gxn_pi;
	precision dpiyy = Iyy_pi  +  Pyy_pi;
	precision dpiyn = Iyn_pi  +  Pyn_pi  -  Gyn_pi;
	precision dpinn = Inn_pi  +  Pnn_pi  -  Gnn_pi;

	precision dpi_mag = 0;

	precision dt_pi = fabs( pt / dpi_mag);

#endif
#ifdef WTZMU
	precision dWtTz = It_W  +  Pt_W  -  Gt_W;
	precision dWxTz = Ix_W  +  Px_W;
	precision dWyTz = Iy_W  +  Py_W;
	precision dWnTz = In_W  +  Pn_W  -  Gn_W;

	precision dWTz_mag = 0;
#endif

	return dt_rate;
}


hydro_time_scales compute_hydro_time_scales(precision t, const hydro_variables * const __restrict__ q, const precision * const __restrict__ e, const fluid_velocity * const __restrict__ u, const fluid_velocity * const __restrict__ up, int nx, int ny, int nz, precision dt, precision dt_prev, precision dx, precision dy, precision dn, precision etabar_const, precision dt_min)
{
	hydro_time_scales dt_hydro;						// default value for hydro time scales
	dt_hydro.dt_CFL = 100;
	dt_hydro.dt_micro = 100;
	dt_hydro.dt_rate = 100;

	precision t2 = t * t;
	int stride_y = nx + 4;							// strides for neighbor cells along x, y, n (stride_x = 1)
	int stride_z = (nx + 4) * (ny + 4);				// stride formulas based from linear_column_index()

	precision ui1[6];								// fluid velocity of neighbor cells along x [i-1, i+1]
	precision uj1[6];								// fluid velocity of neighbor cells along y [j-1, j+1]
	precision uk1[6];								// fluid velocity of neighbor cells along n [k-1, k+1]

	precision vxi[4];								// vx of neighbor cells along x [i-2, i-1, i+1, i+2]
	precision vyj[4];								// vy of neighbor cells along y [j-2, j-1, j+1, j+2]
	precision vnk[4];								// vn of neighbor cells along n [k-2, k-1, k+1, k+2]

	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				int simm = s - 2;			// neighbor cell indices (x)
				int sim  = s - 1;
				int sip  = s + 1;
				int sipp = s + 2;

				int sjmm = s - 2*stride_y;	// neighbor cell indices (y)
				int sjm  = s - stride_y;
				int sjp  = s + stride_y;
				int sjpp = s + 2*stride_y;

				int skmm = s - 2*stride_z;	// neighbor cell indices (n)
				int skm  = s - stride_z;
				int skp  = s + stride_z;
				int skpp = s + 2*stride_z;

				precision e_s = e[s];

				precision ux = u[s].ux;		// current fluid velocity
				precision uy = u[s].uy;
				precision un = u[s].un;
				precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t2 * un * un);

				precision ux_p = up[s].ux;	// previous fluid velocity
				precision uy_p = up[s].uy;
				precision un_p = up[s].un;

				get_fluid_velocity_neighbor_cells(u[simm], u[sim], u[sip], u[sipp], u[sjmm], u[sjm], u[sjp], u[sjpp], u[skmm], u[skm], u[skp], u[skpp], ui1, uj1, uk1, vxi, vyj, vnk, t2);


				// compute the hydrodynamic time scales for fluid cell s
				precision dt_CFL = compute_CFL_time_scale(vxi, vyj, vnk, ut, ux, uy, un, dx, dy, dn);

				precision dt_micro = compute_microscopic_time_scale(e_s, etabar_const);

			#ifdef ANISO_HYDRO
				precision dt_rate = compute_comoving_evolution_rate_aniso_hydro(q[s], e_s, t, ui1, uj1, uk1, ut, ux, uy, un, ux_p, uy_p, un_p, dt_prev, dx, dy, dn, etabar_const, dt_min);
			#else
				printf("compute_hydro_time_scales error: nothing for viscous hydro rate\n");
			#endif

				// update the minimum hydro time scales of the entire grid
				dt_hydro.dt_CFL	= fmin(dt_CFL, dt_hydro.dt_CFL);
				dt_hydro.dt_micro = fmin(dt_micro, dt_hydro.dt_micro);
				dt_hydro.dt_rate = fmin(dt_rate, dt_hydro.dt_rate);
			}
		}
	}
	return dt_hydro;
}






