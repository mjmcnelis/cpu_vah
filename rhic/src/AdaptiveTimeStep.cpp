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

inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


inline precision sign(precision x)
{
	if(x > 0.) return 1.;
	else if(x < 0.) return -1.;
	else return 0.;
}


precision compute_adaptive_time_step(precision t, precision dt_CFL, precision dt_source, precision dt_min)
{
	precision dt = fmin(dt_CFL, dt_source);			// currently have this

	//if(dt < dt_min) printf("compute_adaptive_time_step error: dt = %.3g < %.2g\n", dt, dt_min);

	//printf("compute_adaptive_time_step: %.3g\t%.3g\n", dt_CFL, dt_source);

	//dt = 0.01 * dt_min * floor(100. * dt / dt_min);

	dt = 0.001 * dt_min * floor(1000. * dt / dt_min);	// round dt to numerical precision (dt_min / 100)

	dt = fmax(dt_min, dt);

	FILE * dt_adaptive;
	dt_adaptive = fopen("output/dt_adaptive.dat", "a");
	fprintf(dt_adaptive, "%.8f\t%.8f\n", t, dt);
	fclose(dt_adaptive);

	return dt;
}


precision compute_dt_CFL(precision t, lattice_parameters lattice, hydro_parameters hydro)
{
	precision dt_CFL = 1./0.;

	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision dx = lattice.lattice_spacing_x;
	precision dy = lattice.lattice_spacing_y;
	precision dn = lattice.lattice_spacing_eta;

	precision t2 = t * t;
	precision Theta = hydro.flux_limiter;

	int stride_y = nx + 4;							// strides for neighbor cells along x, y, n (stride_x = 1)
	int stride_z = (nx + 4) * (ny + 4);				// stride formulas based from linear_column_index()

	precision ui1[6], uj1[6], uk1[6];				// these are just filler arguments below

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

				precision ux = u[s].ux;		// current fluid velocity
				precision uy = u[s].uy;
				precision un = u[s].un;
				precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t2 * un * un);

				get_fluid_velocity_neighbor_cells(u[simm], u[sim], u[sip], u[sipp], u[sjmm], u[sjm], u[sjp], u[sjpp], u[skmm], u[skm], u[skp], u[skpp], ui1, uj1, uk1, vxi, vyj, vnk, t2);

				precision ax = compute_max_local_propagation_speed(vxi, ux / ut, Theta);
				precision ay = compute_max_local_propagation_speed(vyj, uy / ut, Theta);
				precision an = compute_max_local_propagation_speed(vnk, un / ut, Theta);

				dt_CFL = fmin(dt_CFL, 0.125 * fmin(dx / ax, fmin(dy / ay, dn / an)));	// take the minimum wave propagation time
			}
		}
	}

	FILE * dt_CFL_bound;
	dt_CFL_bound = fopen("output/dt_CFL.dat", "a");
	fprintf(dt_CFL_bound, "%.4f\t%.8f\n", t, dt_CFL);
	fclose(dt_CFL_bound);

	return dt_CFL;
}


hydro_variables compute_q_star(hydro_variables q, hydro_variables f, precision dt_old)
{
	hydro_variables q_star;		// add euler step using old time step

	q_star.ttt = q.ttt  +  dt_old * f.ttt;
	q_star.ttx = q.ttx  +  dt_old * f.ttx;
	q_star.tty = q.tty  +  dt_old * f.tty;
	q_star.ttn = q.ttn  +  dt_old * f.ttn;

#ifdef ANISO_HYDRO
	q_star.pl  = q.pl  +  dt_old * f.pl;
#if (PT_MATCHING == 1)
	q_star.pt  = q.pt  +  dt_old * f.pt;
#endif

	// should I look at pi_mag?

#endif
#ifdef PIMUNU
	q_star.pitt = q.pitt  +  dt_old * f.pitt;
	q_star.pitx = q.pitx  +  dt_old * f.pitx;
	q_star.pity = q.pity  +  dt_old * f.pity;
	q_star.pitn = q.pitn  +  dt_old * f.pitn;
	q_star.pixx = q.pixx  +  dt_old * f.pixx;
	q_star.pixy = q.pixy  +  dt_old * f.pixy;
	q_star.pixn = q.pixn  +  dt_old * f.pixn;
	q_star.piyy = q.piyy  +  dt_old * f.piyy;
	q_star.piyn = q.piyn  +  dt_old * f.piyn;
	q_star.pinn = q.pinn  +  dt_old * f.pinn;
#endif
#ifdef WTZMU
	q_star.WtTz = q.WtTz  +  dt_old * f.WtTz;
	q_star.WxTz = q.WxTz  +  dt_old * f.WxTz;
	q_star.WyTz = q.WyTz  +  dt_old * f.WyTz;
	q_star.WnTz = q.WnTz  +  dt_old * f.WnTz;
#endif
#ifdef PI
	q_star.Pi = q.Pi  +  dt_old * f.Pi;
#endif

	return q_star;
}

precision compute_q_star_norm(hydro_variables q_star)		// I should be adding things with the same dimension...(which time to use?)
{
	precision ttt = q_star.ttt;
	precision ttx = q_star.ttx;
	precision tty = q_star.tty;
	precision ttn = q_star.ttn;
	precision pl = q_star.pl;

	precision norm2 =  ttt * ttt  +  ttx * ttx  +  tty * tty  +  ttn * ttn  +  pl * pl;

#if (PT_MATCHING == 1)
	precision pt = q_star.pt;

	norm2 += (pt * pt);
#endif
#ifdef PIMUNU
	precision pitt = q_star.pitt;
	precision pitx = q_star.pitx;
	precision pity = q_star.pity;
	precision pitn = q_star.pitn;
	precision pixx = q_star.pixx;
	precision pixy = q_star.pixy;
	precision pixn = q_star.pixn;
	precision piyy = q_star.piyy;
	precision piyn = q_star.piyn;
	precision pinn = q_star.pinn;

	norm2 += (pitt * pitt  +  pitx * pitx  +  pity * pity  +  pitn * pitn  +  pixx * pixx  +  pixy * pixy  +  pixn * pixn  +  piyy * piyy  +  piyn * piyn  +  pinn * pinn);
#endif
#ifdef WTZMU
	precision WtTz = q_star.WtTz;
	precision WxTz = q_star.WxTz;
	precision WyTz = q_star.WyTz;
	precision WnTz = q_star.WnTz;

	norm2 += (WtTz * WtTz  +  WxTz * WxTz  +  WyTz * WyTz  +  WnTz * WnTz);
#endif
#ifdef PI
	precision Pi = q_star.Pi;

	norm2 += (Pi * Pi);
#endif

	return sqrt(norm2);
}


inline precision second_derivative_squared(precision q_prev, precision q, precision q_star)	
{
	precision second_derivative = q_prev  -  2. * q  +  q_star;		// left out factor of 2 / dt_prev^2

	return second_derivative * second_derivative;
}


precision compute_second_derivative_norm(hydro_variables q_prev, hydro_variables q, hydro_variables q_star)
{
	precision norm2 = second_derivative_squared(q_prev.ttt, q.ttt, q_star.ttt);

	norm2 += second_derivative_squared(q_prev.ttx, q.ttx, q_star.ttx);
	norm2 += second_derivative_squared(q_prev.tty, q.tty, q_star.tty);
	norm2 += second_derivative_squared(q_prev.ttn, q.ttn, q_star.ttn);

#ifdef ANISO_HYDRO
	norm2 += second_derivative_squared(q_prev.pl, q.pl, q_star.pl);
#if (PT_MATCHING == 1)
	norm2 += second_derivative_squared(q_prev.pt, q.pt, q_star.pt);
#endif
#endif
#ifdef PIMUNU
	norm2 += second_derivative_squared(q_prev.pitt, q.pitt, q_star.pitt);
	norm2 += second_derivative_squared(q_prev.pitx, q.pitx, q_star.pitx);
	norm2 += second_derivative_squared(q_prev.pity, q.pity, q_star.pity);
	norm2 += second_derivative_squared(q_prev.pitn, q.pitn, q_star.pitn);
	norm2 += second_derivative_squared(q_prev.pixx, q.pixx, q_star.pixx);
	norm2 += second_derivative_squared(q_prev.pixy, q.pixy, q_star.pixy);
	norm2 += second_derivative_squared(q_prev.pixn, q.pixn, q_star.pixn);
	norm2 += second_derivative_squared(q_prev.piyy, q.piyy, q_star.piyy);
	norm2 += second_derivative_squared(q_prev.piyn, q.piyn, q_star.piyn);
	norm2 += second_derivative_squared(q_prev.pinn, q.pinn, q_star.pinn);
#endif
#ifdef WTZMU
	norm2 += second_derivative_squared(q_prev.WtTz, q.WtTz, q_star.WtTz);
	norm2 += second_derivative_squared(q_prev.WxTz, q.WxTz, q_star.WxTz);
	norm2 += second_derivative_squared(q_prev.WyTz, q.WyTz, q_star.WyTz);
	norm2 += second_derivative_squared(q_prev.WnTz, q.WnTz, q_star.WnTz);
#endif
#ifdef PI
	norm2 += second_derivative_squared(q_prev.Pi, q.Pi, q_star.Pi);
#endif

	return sqrt(norm2);
}




precision adaptive_method_norm(precision q_star_norm, precision second_derivative_norm, precision dt_prev, precision delta_0)
{
	precision tolerance = delta_0 * fmax(1.0, q_star_norm);		// shortcut formula 

	return dt_prev * sqrt(tolerance / second_derivative_norm);
}


precision adaptive_method(precision q_prev, precision q, precision q_star, precision f, precision dt_prev, precision delta_0)
{
	// compute the adaptive time step to satisfy either the absolute and relative tolerance

	precision dt_abs = dt_prev * sqrt(fabs(delta_0 / (q_prev  -  2. * q  +  q_star)));

	precision dt_rel = dt_abs * dt_abs / 2. * (sign(q) * f  +  sqrt(f * f  +  4. * fabs(q) / dt_abs / dt_abs));

	return fmax(dt_abs, dt_rel);
}


precision predict_next_time_step(hydro_variables q_prev, hydro_variables q, hydro_variables q_star, hydro_variables f, precision dt_prev, precision delta_0)
{
	// take the smallest adaptive time step of the hydro variables

	precision dt = adaptive_method(q_prev.ttt, q.ttt, q_star.ttt, f.ttt, dt_prev, delta_0);

	dt =  fmin(dt, adaptive_method(q_prev.ttx, q.ttx, q_star.ttx, f.ttx, dt_prev, delta_0));
	dt =  fmin(dt, adaptive_method(q_prev.tty, q.tty, q_star.tty, f.tty, dt_prev, delta_0));
	dt =  fmin(dt, adaptive_method(q_prev.ttn, q.ttn, q_star.ttn, f.ttn, dt_prev, delta_0));

#ifdef ANISO_HYDRO
	dt = fmin(dt, adaptive_method(q_prev.pl, q.pl, q_star.pl, f.pl, dt_prev, delta_0));
#if (PT_MATCHING == 1)
	dt = fmin(dt, adaptive_method(q_prev.pt, q.pt, q_star.pt, f.pt, dt_prev, delta_0));
#endif
#endif
#ifdef PIMUNU
	dt = fmin(dt, adaptive_method(q_prev.pitt, q.pitt, q_star.pitt, f.pitt, dt_prev, delta_0));
	dt = fmin(dt, adaptive_method(q_prev.pitx, q.pitx, q_star.pitx, f.pitx, dt_prev, delta_0));
	dt = fmin(dt, adaptive_method(q_prev.pity, q.pity, q_star.pity, f.pity, dt_prev, delta_0));
	dt = fmin(dt, adaptive_method(q_prev.pitn, q.pitn, q_star.pitn, f.pitn, dt_prev, delta_0));
	dt = fmin(dt, adaptive_method(q_prev.pixx, q.pixx, q_star.pixx, f.pixx, dt_prev, delta_0));
	dt = fmin(dt, adaptive_method(q_prev.pixy, q.pixy, q_star.pixy, f.pixy, dt_prev, delta_0));
	dt = fmin(dt, adaptive_method(q_prev.pixn, q.pixn, q_star.pixn, f.pixn, dt_prev, delta_0));
	dt = fmin(dt, adaptive_method(q_prev.piyy, q.piyy, q_star.piyy, f.piyy, dt_prev, delta_0));
	dt = fmin(dt, adaptive_method(q_prev.piyn, q.piyn, q_star.piyn, f.piyn, dt_prev, delta_0));
	dt = fmin(dt, adaptive_method(q_prev.pinn, q.pinn, q_star.pinn, f.pinn, dt_prev, delta_0));
#endif
#ifdef WTZMU
	dt = fmin(dt, adaptive_method(q_prev.WtTz, q.WtTz, q_star.WtTz, f.WtTz, dt_prev, delta_0));
	dt = fmin(dt, adaptive_method(q_prev.WxTz, q.WxTz, q_star.WxTz, f.WxTz, dt_prev, delta_0));
	dt = fmin(dt, adaptive_method(q_prev.WyTz, q.WyTz, q_star.WyTz, f.WyTz, dt_prev, delta_0));
	dt = fmin(dt, adaptive_method(q_prev.WnTz, q.WnTz, q_star.WnTz, f.WnTz, dt_prev, delta_0));
#endif
#ifdef PI
	dt = fmin(dt, adaptive_method(q_prev.Pi, q.Pi, q_star.Pi, f.Pi, dt_prev, delta_0));
#endif

	return dt;
}


precision compute_dt_source(precision t, const hydro_variables * const __restrict__ q_prev, const hydro_variables * const __restrict__ q, const hydro_variables * const __restrict__ f, precision dt_prev, lattice_parameters lattice)
{
	precision dt_source = 1./0.;

	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision alpha = lattice.alpha;
	precision delta_0 = lattice.delta_0;
	
	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				hydro_variables q_star = compute_q_star(q[s], f[s], dt_prev);

				precision q_star_norm = compute_q_star_norm(q_star);

				precision second_derivative_norm = compute_second_derivative_norm(q_prev[s], q[s], q_star);

				precision dt_next = adaptive_method_norm(q_star_norm, second_derivative_norm, dt_prev, delta_0);


				//precision dt_next = predict_next_time_step(q_prev[s], q[s], q_star, f[s], dt_prev, delta_0);

				dt_source = fmin(dt_source, dt_next);
			}
		}
	}

	FILE * dt_predict;
	dt_predict = fopen("output/dt_source.dat", "a");
	fprintf(dt_predict, "%.4f\t%.8f\n", t, dt_source);
	fclose(dt_predict);

	return fmax((1. -  alpha / 2.) * dt_prev, fmin(dt_source, (1. + alpha) * dt_prev));
}





























