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

	dt = 0.01 * dt_min * floor(100. * dt / dt_min);	// round dt to numerical precision (dt_min / 100)

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

				dt_CFL = fmin(dt_CFL, 0.125 * fmin(dx / ax, fmin(dy / ay, dn / an)));	// take the min wave propagation time
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


precision adaptive_method(precision q_prev, precision q, precision q_star, precision f, precision dt2, precision delta_0)
{
	// compute the adaptive time step to satisfy
	// either the absolute and relative tolerance

	precision second_derivative = 2. * fabs(q_prev  -  2. * q  +  q_star) / dt2;

	precision dt_abs = sqrt(2. * delta_0 / second_derivative);

	precision dt_rel = dt_abs * dt_abs / 2. * (sign(q) * f  +  sqrt(f * f  +  4. * fabs(q) / dt_abs / dt_abs));

	return fmax(dt_abs, dt_rel);
}


precision predict_next_time_step(hydro_variables q_prev, hydro_variables q, hydro_variables q_star, hydro_variables f, precision dt_prev2, precision delta_0)
{
	// take the smallest adaptive time step of the hydro variables

	precision dt = adaptive_method(q_prev.ttt, q.ttt, q_star.ttt, f.ttt, dt_prev2, delta_0);

	dt =  fmin(dt, adaptive_method(q_prev.ttx, q.ttx, q_star.ttx, f.ttx, dt_prev2, delta_0));
	dt =  fmin(dt, adaptive_method(q_prev.tty, q.tty, q_star.tty, f.tty, dt_prev2, delta_0));
	dt =  fmin(dt, adaptive_method(q_prev.ttn, q.ttn, q_star.ttn, f.ttn, dt_prev2, delta_0));

#ifdef ANISO_HYDRO
	dt = fmin(dt, adaptive_method(q_prev.pl, q.pl, q_star.pl, f.pl, dt_prev2, delta_0));
#if (PT_MATCHING == 1)
	dt = fmin(dt, adaptive_method(q_prev.pt, q.pt, q_star.pt, f.pt, dt_prev2, delta_0));
#endif
#endif
#ifdef PIMUNU
	dt = fmin(dt, adaptive_method(q_prev.pitt, q.pitt, q_star.pitt, f.pitt, dt_prev2, delta_0));
	dt = fmin(dt, adaptive_method(q_prev.pitx, q.pitx, q_star.pitx, f.pitx, dt_prev2, delta_0));
	dt = fmin(dt, adaptive_method(q_prev.pity, q.pity, q_star.pity, f.pity, dt_prev2, delta_0));
	dt = fmin(dt, adaptive_method(q_prev.pitn, q.pitn, q_star.pitn, f.pitn, dt_prev2, delta_0));
	dt = fmin(dt, adaptive_method(q_prev.pixx, q.pixx, q_star.pixx, f.pixx, dt_prev2, delta_0));
	dt = fmin(dt, adaptive_method(q_prev.pixy, q.pixy, q_star.pixy, f.pixy, dt_prev2, delta_0));
	dt = fmin(dt, adaptive_method(q_prev.pixn, q.pixn, q_star.pixn, f.pixn, dt_prev2, delta_0));
	dt = fmin(dt, adaptive_method(q_prev.piyy, q.piyy, q_star.piyy, f.piyy, dt_prev2, delta_0));
	dt = fmin(dt, adaptive_method(q_prev.piyn, q.piyn, q_star.piyn, f.piyn, dt_prev2, delta_0));
	dt = fmin(dt, adaptive_method(q_prev.pinn, q.pinn, q_star.pinn, f.pinn, dt_prev2, delta_0));
#endif
#ifdef WTZMU
	dt = fmin(dt, adaptive_method(q_prev.WtTz, q.WtTz, q_star.WtTz, f.WtTz, dt_prev2, delta_0));
	dt = fmin(dt, adaptive_method(q_prev.WxTz, q.WxTz, q_star.WxTz, f.WxTz, dt_prev2, delta_0));
	dt = fmin(dt, adaptive_method(q_prev.WyTz, q.WyTz, q_star.WyTz, f.WyTz, dt_prev2, delta_0));
	dt = fmin(dt, adaptive_method(q_prev.WnTz, q.WnTz, q_star.WnTz, f.WnTz, dt_prev2, delta_0));
#endif
#ifdef PI
	dt = fmin(dt, adaptive_method(q_prev.Pi, q.Pi, q_star.Pi, f.Pi, dt_prev2, delta_0));
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

	precision dt_prev2 = dt_prev * dt_prev;

	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				hydro_variables q_star = compute_q_star(q[s], f[s], dt_prev);

				precision dt_next = predict_next_time_step(q_prev[s], q[s], q_star, f[s], dt_prev2, delta_0);

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





























