#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_poly.h>
#include "../include/Macros.h"
#include "../include/Precision.h"
#include "../include/FluxTerms.h"
#include "../include/DynamicalVariables.h"
#include "../include/NeighborCells.h"
#include "../include/OpenMP.h"

const precision sqrt_variables = sqrt(NUMBER_CONSERVED_VARIABLES - B_FIELD_COMPONENTS);


inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


precision compute_adaptive_time_step(precision t, precision dt_CFL, precision dt_source, precision dt_min)
{
	precision dt = fmin(dt_CFL, dt_source);

	dt = 0.001 * dt_min * floor(1000. * dt / dt_min);	// round dt to numerical precision (dt_min / 100)

	dt = fmax(dt_min, dt);

#ifdef ADAPTIVE_FILE
	FILE * dt_adaptive;
	dt_adaptive = fopen("output/adaptive/dt_adaptive.dat", "a");
	fprintf(dt_adaptive, "%.8f\t%.8f\n", t, dt);
	fclose(dt_adaptive);
#endif

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

	int stride_y = nx + 4;                         // strides for neighbor cells along x, y, n (stride_x = 1)
	int stride_z = (nx + 4) * (ny + 4);            // stride formulas based from linear_column_index()

	#pragma omp parallel for collapse(3) reduction(min:dt_CFL)
	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				precision ui1[6], uj1[6], uk1[6];  // these are just filler arguments below

				precision vxi[4];                  // vx of neighbor cells along x [i-2, i-1, i+1, i+2]
				precision vyj[4];                  // vy of neighbor cells along y [j-2, j-1, j+1, j+2]
				precision vnk[4];                  // vn of neighbor cells along n [k-2, k-1, k+1, k+2]

				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				int simm = s - 2;                  // neighbor cell indices (x)
				int sim  = s - 1;
				int sip  = s + 1;
				int sipp = s + 2;

				int sjmm = s - 2*stride_y;         // neighbor cell indices (y)
				int sjm  = s - stride_y;
				int sjp  = s + stride_y;
				int sjpp = s + 2*stride_y;

				int skmm = s - 2*stride_z;         // neighbor cell indices (n)
				int skm  = s - stride_z;
				int skp  = s + stride_z;
				int skpp = s + 2*stride_z;

				precision ux = u[s].ux;            // current fluid velocity
				precision uy = u[s].uy;
			#ifndef BOOST_INVARIANT
				precision un = u[s].un;
			#else
				precision un = 0;
			#endif
				precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t2 * un * un);

				get_fluid_velocity_neighbor_cells(u[simm], u[sim], u[sip], u[sipp], u[sjmm], u[sjm], u[sjp], u[sjpp], u[skmm], u[skm], u[skp], u[skpp], ui1, uj1, uk1, vxi, vyj, vnk, t2);

				precision ax = compute_max_local_propagation_speed(vxi, ux / ut, Theta);
				precision ay = compute_max_local_propagation_speed(vyj, uy / ut, Theta);
				precision an = compute_max_local_propagation_speed(vnk, un / ut, Theta);

				dt_CFL = fmin(dt_CFL, fmin(dx / ax, fmin(dy / ay, dn / an)));   // take the minimum wave propagation time
			}
		}
	}
#ifdef ADAPTIVE_FILE
	FILE * dt_CFL_bound;
	dt_CFL_bound = fopen("output/adaptive/dt_CFL.dat", "a");
	fprintf(dt_CFL_bound, "%.8f\t%.8f\n", t, dt_CFL / 8.);
	fclose(dt_CFL_bound);
#endif

	return dt_CFL / 8.;
}


hydro_variables compute_q_star(hydro_variables q, hydro_variables f, precision dt_prev)
{
	hydro_variables q_star;                 // add euler step using old time step

	q_star.ttt = q.ttt  +  dt_prev * f.ttt;
	q_star.ttx = q.ttx  +  dt_prev * f.ttx;
	q_star.tty = q.tty  +  dt_prev * f.tty;
#ifndef BOOST_INVARIANT
	q_star.ttn = q.ttn  +  dt_prev * f.ttn;
#endif

#ifdef ANISO_HYDRO
	q_star.pl  = q.pl  +  dt_prev * f.pl;
	q_star.pt  = q.pt  +  dt_prev * f.pt;   // figure out what to do with this later
#endif

#ifdef PIMUNU
	q_star.pitt = q.pitt  +  dt_prev * f.pitt;
	q_star.pitx = q.pitx  +  dt_prev * f.pitx;
	q_star.pity = q.pity  +  dt_prev * f.pity;
	q_star.pixx = q.pixx  +  dt_prev * f.pixx;
	q_star.pixy = q.pixy  +  dt_prev * f.pixy;
	q_star.piyy = q.piyy  +  dt_prev * f.piyy;

#ifndef BOOST_INVARIANT
	q_star.pitn = q.pitn  +  dt_prev * f.pitn;
	q_star.pixn = q.pixn  +  dt_prev * f.pixn;
	q_star.piyn = q.piyn  +  dt_prev * f.piyn;
	q_star.pinn = q.pinn  +  dt_prev * f.pinn;
#else
	#ifndef ANISO_HYDRO
		q_star.pinn = q.pinn  +  dt_prev * f.pinn;
	#endif
#endif

#endif

#ifdef WTZMU
	q_star.WtTz = q.WtTz  +  dt_prev * f.WtTz;
	q_star.WxTz = q.WxTz  +  dt_prev * f.WxTz;
	q_star.WyTz = q.WyTz  +  dt_prev * f.WyTz;
	q_star.WnTz = q.WnTz  +  dt_prev * f.WnTz;
#endif

#ifdef PI
	q_star.Pi = q.Pi  +  dt_prev * f.Pi;
#endif

	return q_star;
}


precision compute_hydro_norm2(hydro_variables q)   // should probably add things with the same dimension...
{
	precision norm2 = q.ttt * q.ttt  +  q.ttx * q.ttx  +  q.tty * q.tty;

#ifndef BOOST_INVARIANT
	norm2 +=  q.ttn * q.ttn;
#endif

#ifdef ANISO_HYDRO
	norm2 += (q.pl * q.pl);
	norm2 += (q.pt * q.pt);
#endif

#ifdef PIMUNU
	norm2 += (q.pitt * q.pitt  +  q.pitx * q.pitx  +  q.pity * q.pity  +  q.pixx * q.pixx  +  q.pixy * q.pixy  +  q.piyy * q.piyy);

#ifndef BOOST_INVARIANT
	norm2 += (q.pitn * q.pitn  +  q.pixn * q.pixn  +  q.piyn * q.piyn  +  q.pinn * q.pinn);
#else
	#ifndef ANISO_HYDRO
		norm2 += (q.pinn * q.pinn);
	#endif
#endif

#endif

#ifdef WTZMU
	norm2 += (q.WtTz * q.WtTz  +  q.WxTz * q.WxTz  +  q.WyTz * q.WyTz  +  q.WnTz * q.WnTz);
#endif

#ifdef PI
	norm2 += (q.Pi * q.Pi);
#endif

	return norm2;
}


precision second_derivative_squared(precision q_prev, precision q, precision q_star)
{
	return (q_prev  -  2. * q  +  q_star) * (q_prev  -  2. * q  +  q_star);    // left out prefactor 2 / dt_prev^2
}


precision compute_second_derivative_norm(hydro_variables q_prev, hydro_variables q, hydro_variables q_star)
{
	precision norm2 = (	second_derivative_squared(q_prev.ttt, q.ttt, q_star.ttt) +
						second_derivative_squared(q_prev.ttx, q.ttx, q_star.ttx) +
						second_derivative_squared(q_prev.tty, q.tty, q_star.tty));

#ifndef BOOST_INVARIANT
	norm2 += second_derivative_squared(q_prev.ttn, q.ttn, q_star.ttn);
#endif

#ifdef ANISO_HYDRO
	norm2 += second_derivative_squared(q_prev.pl, q.pl, q_star.pl);
	norm2 += second_derivative_squared(q_prev.pt, q.pt, q_star.pt);
#endif

#ifdef PIMUNU
	norm2 += (	second_derivative_squared(q_prev.pitt, q.pitt, q_star.pitt)	+
				second_derivative_squared(q_prev.pitx, q.pitx, q_star.pitx)	+
				second_derivative_squared(q_prev.pity, q.pity, q_star.pity)	+
				second_derivative_squared(q_prev.pixx, q.pixx, q_star.pixx)	+
				second_derivative_squared(q_prev.pixy, q.pixy, q_star.pixy)	+
				second_derivative_squared(q_prev.piyy, q.piyy, q_star.piyy));
#ifndef BOOST_INVARIANT
	norm2 += (	second_derivative_squared(q_prev.pitn, q.pitn, q_star.pitn)	+
				second_derivative_squared(q_prev.pixn, q.pixn, q_star.pixn)	+
				second_derivative_squared(q_prev.piyn, q.piyn, q_star.piyn) +
				second_derivative_squared(q_prev.pinn, q.pinn, q_star.pinn));
#else
	#ifndef ANISO_HYDRO
		norm2 += second_derivative_squared(q_prev.pinn, q.pinn, q_star.pinn);
	#endif
#endif
#endif

#ifdef WTZMU
	norm2 += (	second_derivative_squared(q_prev.WtTz, q.WtTz, q_star.WtTz)	+
				second_derivative_squared(q_prev.WxTz, q.WxTz, q_star.WxTz)	+
				second_derivative_squared(q_prev.WyTz, q.WyTz, q_star.WyTz)	+
				second_derivative_squared(q_prev.WnTz, q.WnTz, q_star.WnTz));
#endif

#ifdef PI
	norm2 += second_derivative_squared(q_prev.Pi, q.Pi, q_star.Pi);
#endif

	return sqrt(norm2);
}


precision dot_product(hydro_variables q, hydro_variables f)
{
	precision dot = q.ttt * f.ttt  +  q.ttx * f.ttx  +  q.tty * f.tty;

#ifndef BOOST_INVARIANT
	dot += (q.ttn * f.ttn);
#endif

#ifdef ANISO_HYDRO
	dot += (q.pl * f.pl);
	dot += (q.pt * f.pt);

#endif

#ifdef PIMUNU
	dot += (q.pitt * f.pitt  +  q.pitx * f.pitx  +  q.pity * f.pity  +  q.pixx * f.pixx  +  q.pixy * f.pixy  +  q.piyy * f.piyy);
#ifndef BOOST_INVARIANT
	dot += (q.pitn * f.pitn  +  q.pixn * f.pixn  +  q.piyn * f.piyn  +  q.pinn * f.pinn);
#else
	#ifndef ANISO_HYDRO
		dot += (q.pinn * f.pinn);
	#endif
#endif
#endif

#ifdef WTZMU
	dot += (q.WtTz * f.WtTz  +  q.WxTz * f.WxTz  +  q.WyTz * f.WyTz  +  q.WnTz * f.WnTz);
#endif

#ifdef PI
	dot += (q.Pi * f.Pi);
#endif

	return dot;
}


precision adaptive_method_norm(precision q_norm2, precision f_norm2, precision second_derivative_norm, precision q_dot_f, precision dt_prev, precision delta_0)
{
	precision dt_abs = dt_prev * sqrt(delta_0 * sqrt_variables / second_derivative_norm);

	precision dt_abs2 = dt_abs * dt_abs / sqrt_variables;   // additionally divide out sqrt_variables
	precision dt_abs4 = dt_abs2 * dt_abs2;

	precision c = f_norm2 * dt_abs4;                        // quartic equation: x^4 - c.x^2 - b.x - a = 0  (x = dt_rel)
	precision b = 2. * q_dot_f * dt_abs4;
	precision a = q_norm2 * dt_abs4;

	precision y0, y1, y2;                                   // resolvent cubic equation: y^3 + d.y^2 + e.y + f = 0  (d,e,f)
	int cubic_roots = gsl_poly_solve_cubic(-2. * c,  c * c  +  4. * a,  - b * b, &y0, &y1, &y2);

	if(y0 > 0.)                                             // take the greatest positive solution (see gsl manual)
	{
		precision u = sqrt(y0);                             // quartic equation factored:
                                                            // x^4 - c.x^2 - b.x - a = (x^2 - u.x + t)(x^2 + u.x + v)
		precision discriminant_1 = 2. * (c  +  b / u)  -  y0;
		precision discriminant_2 = 2. * (c  -  b / u)  -  y0;

		precision dt_rel = dt_abs;

		if(discriminant_1 > 0)
		{
			dt_rel = fmax(dt_rel, (u + sqrt(discriminant_1)) / 2.);
		}
		if(discriminant_2 > 0)
		{
			dt_rel = fmax(dt_rel, (-u + sqrt(discriminant_2)) / 2.);
		}
		return dt_rel;
	}

	return dt_abs;
}


precision compute_dt_source(precision t, const hydro_variables * const __restrict__ q_prev, const hydro_variables * const __restrict__ q, const hydro_variables * const __restrict__ f, precision dt_prev, lattice_parameters lattice)
{
	precision dt_source = 1./0.;

	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision alpha = lattice.alpha;
	precision delta_0 = lattice.delta_0;

	#pragma omp parallel for collapse(3) reduction(min:dt_source)
	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				hydro_variables q_star = compute_q_star(q[s], f[s], dt_prev);

				precision q_norm = compute_hydro_norm2(q[s]);
				precision f_norm = compute_hydro_norm2(f[s]);
				precision q_dot_f = dot_product(q[s], f[s]);

				precision second_derivative_norm = compute_second_derivative_norm(q_prev[s], q[s], q_star);

				dt_source = fmin(dt_source, adaptive_method_norm(q_norm, f_norm, second_derivative_norm, q_dot_f, dt_prev, delta_0));
			}
		}
	}

	dt_source = fmax((1. -  alpha) * dt_prev, fmin(dt_source, (1. + alpha) * dt_prev));

#ifdef ADAPTIVE_FILE
	FILE * dt_predict;
	dt_predict = fopen("output/adaptive/dt_source.dat", "a");
	fprintf(dt_predict, "%.4f\t%.8f\n", t, dt_source);
	fclose(dt_predict);
#endif

	return dt_source;
}




