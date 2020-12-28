#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../include/Macros.h"
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"
#include "../include/AnisoVariables.h"
#include "../include/GhostCells.h"
#include "../include/FluxTerms.h"
#include "../include/SourceTerms.h"
#include "../include/InferredVariables.h"
#include "../include/NeighborCells.h"
#include "../include/Regulation.h"
#include "../include/OpenMP.h"

inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


void euler_step(precision t, const hydro_variables * const __restrict__ q_current, hydro_variables * const __restrict__ q_update,
const precision * const __restrict__ e_current, const precision * const __restrict__ lambda_current, const precision * const __restrict__ aT_current, const precision * const __restrict__ aL_current, const fluid_velocity * const __restrict__ u_previous, const fluid_velocity * const __restrict__ u_current, precision dt, precision dt_prev, lattice_parameters lattice, hydro_parameters hydro, int update, int RK2)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision dx = lattice.lattice_spacing_x;
	precision dy = lattice.lattice_spacing_y;
	precision dn = lattice.lattice_spacing_eta;

	precision Theta = hydro.flux_limiter;

	precision t2 = t * t;

	int stride_y = nx + 4;                                                  // strides for neighbor cells along x, y, n (stride_x = 1)
	int stride_z = (nx + 4) * (ny + 4);                                     // stride formulas based from linear_column_index()

	#pragma omp parallel for collapse(3)
	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				precision qs[NUMBER_CONSERVED_VARIABLES];       // dynamical variables at cell s
				precision E[NUMBER_CONSERVED_VARIABLES];        // total source function E at cell s
				precision  S[NUMBER_CONSERVED_VARIABLES];       // external source terms S at cell s

				precision e1[6];                                // energy density of 6 neighbor cells {i-1, i+1, j-1, j+1, k-1, k+1}

				precision qi1[2 * NUMBER_CONSERVED_VARIABLES];  // dynamical variables of 2 neighbor cells along x {i-1, i+1}
				precision qj1[2 * NUMBER_CONSERVED_VARIABLES];  // dynamical variables of 2 neighbor cells along y {j-1, j+1}
				precision qk1[2 * NUMBER_CONSERVED_VARIABLES];  // dynamical variables of 2 neighbor cells along n {k-1, k+1}

				precision qi2[2 * NUMBER_CONSERVED_VARIABLES];  // dynamical variables of 2 neighbor cells along x {i-2, i+2}
				precision qj2[2 * NUMBER_CONSERVED_VARIABLES];  // dynamical variables of 2 neighbor cells along y {j-2, j+2}
				precision qk2[2 * NUMBER_CONSERVED_VARIABLES];  // dynamical variables of 2 neighbor cells along n {k-2, k+2}

				precision ui1[6];                               // fluid velocity of 2 neighbor cells along x {i-1, i+1}
				precision uj1[6];                               // fluid velocity of 2 neighbor cells along y {j-1, j+1}
				precision uk1[6];                               // fluid velocity of 2 neighbor cells along n {k-1, k+1}

				precision vxi[4];                               // vx of 4 neighbor cells along x {i-2, i-1, i+1, i+2}
				precision vyj[4];                               // vy of 4 neighbor cells along y {j-2, j-1, j+1, j+2}
				precision vnk[4];                               // vn of 4 neighbor cells along n {k-2, k-1, k+1, k+2}

				precision  Hx_plus[NUMBER_CONSERVED_VARIABLES]; // Hx_{i + 1/2}
				precision Hx_minus[NUMBER_CONSERVED_VARIABLES]; // Hx_{i - 1/2}

				precision  Hy_plus[NUMBER_CONSERVED_VARIABLES]; // Hy_{j + 1/2}
				precision Hy_minus[NUMBER_CONSERVED_VARIABLES]; // Hy_{j - 1/2}

				precision  Hn_plus[NUMBER_CONSERVED_VARIABLES]; // Hn_{k + 1/2}
				precision Hn_minus[NUMBER_CONSERVED_VARIABLES]; // Hn_{k - 1/2}

				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				int simm = s - 2;                               // neighbor cell indices (x)
				int sim  = s - 1;
				int sip  = s + 1;
				int sipp = s + 2;

				int sjmm = s - 2*stride_y;                      // neighbor cell indices (y)
				int sjm  = s - stride_y;
				int sjp  = s + stride_y;
				int sjpp = s + 2*stride_y;

				int skmm = s - 2*stride_z;                      // neighbor cell indices (n)
				int skm  = s - stride_z;
				int skp  = s + stride_z;
				int skpp = s + 2*stride_z;

				qs[0] = q_current[s].ttt;                       // dynamical variables of cell s
				qs[1] = q_current[s].ttx;
				qs[2] = q_current[s].tty;

				int a = 3;
			#ifndef BOOST_INVARIANT
				qs[a] = q_current[s].ttn; a++;
			#endif

			#ifdef ANISO_HYDRO
				qs[a] = q_current[s].pl; a++;
				qs[a] = q_current[s].pt; a++;
			#endif

			#ifdef B_FIELD
				qs[a] = q_current[s].b; a++;
			#endif

			#ifdef PIMUNU
				qs[a] = q_current[s].pitt; a++;
				qs[a] = q_current[s].pitx; a++;
				qs[a] = q_current[s].pity; a++;
			#ifndef BOOST_INVARIANT
				qs[a] = q_current[s].pitn; a++;
			#endif
				qs[a] = q_current[s].pixx; a++;
				qs[a] = q_current[s].pixy; a++;
			#ifndef BOOST_INVARIANT
				qs[a] = q_current[s].pixn; a++;
			#endif
				qs[a] = q_current[s].piyy; a++;
			#ifndef BOOST_INVARIANT
				qs[a] = q_current[s].piyn; a++;
				qs[a] = q_current[s].pinn; a++;
			#else
			#ifndef ANISO_HYDRO
				qs[a] = q_current[s].pinn; a++;
			#endif
			#endif
			#endif

			#ifdef WTZMU
				qs[a] = q_current[s].WtTz; a++;
				qs[a] = q_current[s].WxTz; a++;
				qs[a] = q_current[s].WyTz; a++;
				qs[a] = q_current[s].WnTz; a++;
			#endif

			#ifdef PI
				qs[a] = q_current[s].Pi; a++;
			#endif

				precision e_s = e_current[s];                   // energy density

			#ifdef ANISO_HYDRO
			#ifdef LATTICE_QCD
				precision lambda_s = lambda_current[s];         // anisotropic variables
				precision aT_s = aT_current[s];
				precision aL_s = aL_current[s];
			#else
				precision lambda_s = 0;
				precision aT_s = 0;
				precision aL_s = 0;
			#endif
			#endif

				precision ux = u_current[s].ux;                 // current fluid velocity
				precision uy = u_current[s].uy;
			#ifndef BOOST_INVARIANT
				precision un = u_current[s].un;
			#else
				precision un = 0;
			#endif

				precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t2 * un * un);

				precision ux_p = u_previous[s].ux;              // previous fluid velocity
				precision uy_p = u_previous[s].uy;
			#ifndef BOOST_INVARIANT
				precision un_p = u_previous[s].un;
			#else
				precision un_p = 0;
			#endif

				// for numerical spatial derivatives
				get_energy_density_neighbor_cells(e_current, e1, sim, sip, sjm, sjp, skm, skp);

				get_fluid_velocity_neighbor_cells(u_current[simm], u_current[sim], u_current[sip], u_current[sipp], u_current[sjmm], u_current[sjm], u_current[sjp], u_current[sjpp], u_current[skmm], u_current[skm], u_current[skp], u_current[skpp], ui1, uj1, uk1, vxi, vyj, vnk, t2);

				get_hydro_variables_neighbor_cells(q_current[sim],  q_current[sip],  qi1);
				get_hydro_variables_neighbor_cells(q_current[simm], q_current[sipp], qi2);

				get_hydro_variables_neighbor_cells(q_current[sjm],  q_current[sjp],  qj1);
				get_hydro_variables_neighbor_cells(q_current[sjmm], q_current[sjpp], qj2);

				get_hydro_variables_neighbor_cells(q_current[skm],  q_current[skp],  qk1);
				get_hydro_variables_neighbor_cells(q_current[skmm], q_current[skpp], qk2);

				// compute source term S
			#ifdef ANISO_HYDRO
				source_terms_aniso_hydro(S, qs, e_s, lambda_s, aT_s, aL_s, t, qi1, qj1, qk1, e1, ui1, uj1, uk1, ux, uy, un, ux_p, uy_p, un_p, dt_prev, dx, dy, dn, hydro);
			#else
				source_terms_viscous_hydro(S, qs, e_s, t, qi1, qj1, qk1, e1, ui1, uj1, uk1, ux, uy, un, ux_p, uy_p, un_p, dt_prev, dx, dy, dn, hydro);
			#endif

				// compute flux terms (H)
				flux_terms(Hx_plus, Hx_minus, qs, qi1, qi2, vxi, ux / ut, Theta);
				flux_terms(Hy_plus, Hy_minus, qs, qj1, qj2, vyj, uy / ut, Theta);
				flux_terms(Hn_plus, Hn_minus, qs, qk1, qk2, vnk, un / ut, Theta);

				// store the results
				if(!update)                                     // store total source function qI <== E
				{
					for(int n = 0; n < NUMBER_CONSERVED_VARIABLES; n++)
					{
						E[n] = S[n]  +  (Hx_minus[n] - Hx_plus[n]) / dx  +  (Hy_minus[n] - Hy_plus[n]) / dy  +  (Hn_minus[n] - Hn_plus[n]) / dn;
					}

					q_update[s].ttt = E[0];
					q_update[s].ttx = E[1];
					q_update[s].tty = E[2];

					a = 3;
				#ifndef BOOST_INVARIANT
					q_update[s].ttn = E[a]; a++;
				#endif

				#ifdef ANISO_HYDRO
					q_update[s].pl = E[a]; a++;
					q_update[s].pt = E[a]; a++;
				#endif

				#ifdef B_FIELD
					q_update[s].b = E[a]; a++;
				#endif

				#ifdef PIMUNU
					q_update[s].pitt = E[a]; a++;
					q_update[s].pitx = E[a]; a++;
					q_update[s].pity = E[a]; a++;
				#ifndef BOOST_INVARIANT
					q_update[s].pitn = E[a]; a++;
				#endif
					q_update[s].pixx = E[a]; a++;
					q_update[s].pixy = E[a]; a++;
				#ifndef BOOST_INVARIANT
					q_update[s].pixn = E[a]; a++;
				#endif
					q_update[s].piyy = E[a]; a++;
				#ifndef BOOST_INVARIANT
					q_update[s].piyn = E[a]; a++;
					q_update[s].pinn = E[a]; a++;
				#else
				#ifndef ANISO_HYDRO
					q_update[s].pinn = E[a]; a++;
				#endif
				#endif
				#endif

				#ifdef WTZMU
					q_update[s].WtTz = E[a]; a++;
					q_update[s].WxTz = E[a]; a++;
					q_update[s].WyTz = E[a]; a++;
					q_update[s].WnTz = E[a]; a++;
				#endif

				#ifdef PI
					q_update[s].Pi = E[a]; a++;
				#endif
				}
				else if(!RK2)                                   // store first intermediate euler step qI <== q + E.dt
				{
					for(int n = 0; n < NUMBER_CONSERVED_VARIABLES; n++)
					{
						qs[n] += dt * (S[n]  +  (Hx_minus[n] - Hx_plus[n]) / dx  +  (Hy_minus[n] - Hy_plus[n]) / dy  +  (Hn_minus[n] - Hn_plus[n]) / dn);
					}

					q_update[s].ttt = qs[0];
					q_update[s].ttx = qs[1];
					q_update[s].tty = qs[2];

					a = 3;
				#ifndef BOOST_INVARIANT
					q_update[s].ttn = qs[a]; a++;
				#endif

				#ifdef ANISO_HYDRO
					q_update[s].pl = qs[a]; a++;
					q_update[s].pt = qs[a]; a++;
				#endif

				#ifdef B_FIELD
					q_update[s].b = qs[a]; a++;
				#endif

				#ifdef PIMUNU
					q_update[s].pitt = qs[a]; a++;
					q_update[s].pitx = qs[a]; a++;
					q_update[s].pity = qs[a]; a++;
				#ifndef BOOST_INVARIANT
					q_update[s].pitn = qs[a]; a++;
				#endif
					q_update[s].pixx = qs[a]; a++;
					q_update[s].pixy = qs[a]; a++;
				#ifndef BOOST_INVARIANT
					q_update[s].pixn = qs[a]; a++;
				#endif
					q_update[s].piyy = qs[a]; a++;
				#ifndef BOOST_INVARIANT
					q_update[s].piyn = qs[a]; a++;
					q_update[s].pinn = qs[a]; a++;
				#else
				#ifndef ANISO_HYDRO
					q_update[s].pinn = qs[a]; a++;
				#endif
				#endif
				#endif

				#ifdef WTZMU
					q_update[s].WtTz = qs[a]; a++;
					q_update[s].WxTz = qs[a]; a++;
					q_update[s].WyTz = qs[a]; a++;
					q_update[s].WnTz = qs[a]; a++;
				#endif

				#ifdef PI
					q_update[s].Pi = qs[a]; a++;
				#endif
				}
				else                                            // store RK2 update Q <== (q + (qI + EI.dt))/2
				{
					for(int n = 0; n < NUMBER_CONSERVED_VARIABLES; n++)
					{
						qs[n] += dt * (S[n]  +  (Hx_minus[n] - Hx_plus[n]) / dx  +  (Hy_minus[n] - Hy_plus[n]) / dy  +  (Hn_minus[n] - Hn_plus[n]) / dn);
					}

					q_update[s].ttt = (q[s].ttt  +  qs[0]) / 2.;    // here q is the extern variable, not
					q_update[s].ttx = (q[s].ttx  +  qs[1]) / 2.;    // one of the euler_step(..) arguments
					q_update[s].tty = (q[s].tty  +  qs[2]) / 2.;

					a = 3;
				#ifndef BOOST_INVARIANT
					q_update[s].ttn = (q[s].ttn  +  qs[a]) / 2.; a++;
				#endif

				#ifdef ANISO_HYDRO
					q_update[s].pl = (q[s].pl  +  qs[a]) / 2.; a++;
					q_update[s].pt = (q[s].pt  +  qs[a]) / 2.; a++;
				#endif

				#ifdef B_FIELD
					q_update[s].b = (q[s].b  +  qs[a]) / 2.; a++;
				#endif

				#ifdef PIMUNU
					q_update[s].pitt = (q[s].pitt  +  qs[a]) / 2.; a++;
					q_update[s].pitx = (q[s].pitx  +  qs[a]) / 2.; a++;
					q_update[s].pity = (q[s].pity  +  qs[a]) / 2.; a++;
				#ifndef BOOST_INVARIANT
					q_update[s].pitn = (q[s].pitn  +  qs[a]) / 2.; a++;
				#endif
					q_update[s].pixx = (q[s].pixx  +  qs[a]) / 2.; a++;
					q_update[s].pixy = (q[s].pixy  +  qs[a]) / 2.; a++;
				#ifndef BOOST_INVARIANT
					q_update[s].pixn = (q[s].pixn  +  qs[a]) / 2.; a++;
				#endif
					q_update[s].piyy = (q[s].piyy  +  qs[a]) / 2.; a++;
				#ifndef BOOST_INVARIANT
					q_update[s].piyn = (q[s].piyn  +  qs[a]) / 2.; a++;
					q_update[s].pinn = (q[s].pinn  +  qs[a]) / 2.; a++;
				#else
				#ifndef ANISO_HYDRO
					q_update[s].pinn = (q[s].pinn  +  qs[a]) / 2.; a++;
				#endif
				#endif
				#endif

				#ifdef WTZMU
					q_update[s].WtTz = (q[s].WtTz  +  qs[a]) / 2.; a++;
					q_update[s].WxTz = (q[s].WxTz  +  qs[a]) / 2.; a++;
					q_update[s].WyTz = (q[s].WyTz  +  qs[a]) / 2.; a++;
					q_update[s].WnTz = (q[s].WnTz  +  qs[a]) / 2.; a++;
				#endif

				#ifdef PI
					q_update[s].Pi = (q[s].Pi  +  qs[a]) / 2.; a++;
				#endif
				}
			}
		}
	}
}


void recompute_euler_step(const hydro_variables * const __restrict__ q_current, hydro_variables * const __restrict__ q_update, precision dt, lattice_parameters lattice)
 {
 	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	// recompute qI after computing adaptive step

	#pragma omp parallel for collapse(3)
	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				q_update[s].ttt = q_current[s].ttt  +  dt * q_update[s].ttt;    // qI <== (q + E.dt_new)
				q_update[s].ttx = q_current[s].ttx  +  dt * q_update[s].ttx;
				q_update[s].tty = q_current[s].tty  +  dt * q_update[s].tty;
			#ifndef BOOST_INVARIANT
				q_update[s].ttn = q_current[s].ttn  +  dt * q_update[s].ttn;
			#endif

			#ifdef ANISO_HYDRO
				q_update[s].pl  = q_current[s].pl  +  dt * q_update[s].pl;
				q_update[s].pt  = q_current[s].pt  +  dt * q_update[s].pt;
			#endif

			#ifdef B_FIELD
				q_update[s].b  = q_current[s].b  +  dt * q_update[s].b;
			#endif

			#ifdef PIMUNU
				q_update[s].pitt = q_current[s].pitt  +  dt * q_update[s].pitt;
				q_update[s].pitx = q_current[s].pitx  +  dt * q_update[s].pitx;
				q_update[s].pity = q_current[s].pity  +  dt * q_update[s].pity;
			#ifndef BOOST_INVARIANT
				q_update[s].pitn = q_current[s].pitn  +  dt * q_update[s].pitn;
			#endif
				q_update[s].pixx = q_current[s].pixx  +  dt * q_update[s].pixx;
				q_update[s].pixy = q_current[s].pixy  +  dt * q_update[s].pixy;
			#ifndef BOOST_INVARIANT
				q_update[s].pixn = q_current[s].pixn  +  dt * q_update[s].pixn;
			#endif
				q_update[s].piyy = q_current[s].piyy  +  dt * q_update[s].piyy;
			#ifndef BOOST_INVARIANT
				q_update[s].piyn = q_current[s].piyn  +  dt * q_update[s].piyn;
				q_update[s].pinn = q_current[s].pinn  +  dt * q_update[s].pinn;
			#else
			#ifndef ANISO_HYDRO
				q_update[s].pinn = q_current[s].pinn  +  dt * q_update[s].pinn;
			#endif
			#endif
			#endif

			#ifdef WTZMU
				q_update[s].WtTz = q_current[s].WtTz  +  dt * q_update[s].WtTz;
				q_update[s].WxTz = q_current[s].WxTz  +  dt * q_update[s].WxTz;
				q_update[s].WyTz = q_current[s].WyTz  +  dt * q_update[s].WyTz;
				q_update[s].WnTz = q_current[s].WnTz  +  dt * q_update[s].WnTz;
			#endif

			#ifdef PI
				q_update[s].Pi = q_current[s].Pi  +  dt * q_update[s].Pi;
			#endif
			}
		}
	}
}


void evolve_hydro_one_time_step(int n, precision t, precision dt, precision dt_prev, lattice_parameters lattice, hydro_parameters hydro, int update, bool hit_CFL)
{
	// first intermediate euler step
	int RK2 = 0;

	if(lattice.adaptive_time_step && !hit_CFL)                                      // qI <== q + E.dt
	{
		if(n == 0)
		{
			euler_step(t, q, qI, e, lambda, aT, aL, up, u, dt, dt_prev, lattice, hydro, update, RK2);
		}
		else
		{
			recompute_euler_step(q, qI, dt, lattice);
		}
	}
	else
	{
		euler_step(t, q, qI, e, lambda, aT, aL, up, u, dt, dt_prev, lattice, hydro, update, RK2);
	}

	t += dt;

	swap_fluid_velocity(&u, &up);                                                   // swap u <==> up

#ifdef ANISO_HYDRO
	set_inferred_variables_aniso_hydro(qI, e, u, t, lattice, hydro);                // compute aniso (e, u)
	// todo: try regulating before solving for aniso variables (haven't done yet)
#ifdef LATTICE_QCD
	set_anisotropic_variables(qI, e, lambda, aT, aL, lattice, hydro);               // compute aniso X = (lambda, aT, aL)
#endif
	regulate_residual_currents(t, qI, e, u, lattice, hydro, RK2);                   // regulate aniso qI
#else
	set_inferred_variables_viscous_hydro(qI, e, u, t, lattice, hydro);              // compute viscous (e, u)
	regulate_viscous_currents(t, qI, e, u, lattice, hydro, RK2);                    // regulate viscous qI
#endif
	set_ghost_cells(qI, e, u, lattice);                                             // for (qI, e, u)

	// second intermediate Euler step
	dt_prev = dt;                                                                   // update previous time step
	RK2 = 1;                                                                        // so that Q <== RK2 update

	euler_step(t, qI, Q, e, lambda, aT, aL, up, u, dt, dt_prev, lattice, hydro, update, RK2);

#ifdef ANISO_HYDRO
	set_inferred_variables_aniso_hydro(Q, e, u, t, lattice, hydro);                 // compute aniso (e, u)
#ifdef LATTICE_QCD
	set_anisotropic_variables(Q, e, lambda, aT, aL, lattice, hydro);                // compute aniso X = (lambda, aT, aL)
#endif
	regulate_residual_currents(t, Q, e, u, lattice, hydro, RK2);                    // regulate aniso Q
#else
	set_inferred_variables_viscous_hydro(Q, e, u, t, lattice, hydro);               // compute viscous (e, u)
	regulate_viscous_currents(t, Q, e, u, lattice, hydro, RK2);                     // regulate viscous Q
#endif
	swap_hydro_variables(&q, &Q);                                                   // swap q <==> Q
	set_ghost_cells(q, e, u, lattice);                                              // for (q, e, u)
}







