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
#include "../include/RungeKutta.h"
#include "../include/OpenMP.h"


inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


void euler_step(precision t, const hydro_variables * const __restrict__ q_current, hydro_variables * const __restrict__ q_update,
const precision * const __restrict__ e_current, const precision * const __restrict__ lambda_current, const precision * const __restrict__ aT_current, const precision * const __restrict__ aL_current, const fluid_velocity * const __restrict__ u_previous, const fluid_velocity * const __restrict__ u_current, precision dt, precision dt_prev, lattice_parameters lattice, hydro_parameters hydro, int stage)
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

				if(!stage)							// get source function q^(1) <-- E(q)
				{
					for(int n = 0; n < NUMBER_CONSERVED_VARIABLES; n++)
					{
						qs[n] = S[n]  +  (Hx_minus[n] - Hx_plus[n]) / dx  +  (Hy_minus[n] - Hy_plus[n]) / dy  +  (Hn_minus[n] - Hn_plus[n]) / dn;
					}
				}
				else								// get intermediate Euler step q^(l) <-- dt.E(q^(l-1))
				{
					for(int n = 0; n < NUMBER_CONSERVED_VARIABLES; n++)
					{
						qs[n] = dt * (S[n]  +  (Hx_minus[n] - Hx_plus[n]) / dx  +  (Hy_minus[n] - Hy_plus[n]) / dy  +  (Hn_minus[n] - Hn_plus[n]) / dn);
					}
				}

				// store results

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
		}
	}
}


void recompute_euler_step(hydro_variables * const __restrict__ q_update, precision dt, lattice_parameters lattice)
 {
 	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	// recompute first intermediate Euler step after computing adaptive step dt
	// q_update is passed E(q)

	#pragma omp parallel for collapse(3)
	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				q_update[s].ttt *= dt;    // q1 <== (q + E.dt_new)
				q_update[s].ttx *= dt;
				q_update[s].tty *= dt;
			#ifndef BOOST_INVARIANT
				q_update[s].ttn *= dt;
			#endif

			#ifdef ANISO_HYDRO
				q_update[s].pl *= dt;
				q_update[s].pt *= dt;
			#endif

			#ifdef B_FIELD
				q_update[s].b *= dt;
			#endif

			#ifdef PIMUNU
				q_update[s].pitt *= dt;
				q_update[s].pitx *= dt;
				q_update[s].pity *= dt;
			#ifndef BOOST_INVARIANT
				q_update[s].pitn *= dt;
			#endif
				q_update[s].pixx *= dt;
				q_update[s].pixy *= dt;
			#ifndef BOOST_INVARIANT
				q_update[s].pixn *= dt;
			#endif
				q_update[s].piyy *= dt;
			#ifndef BOOST_INVARIANT
				q_update[s].piyn *= dt;
				q_update[s].pinn *= dt;
			#else
			#ifndef ANISO_HYDRO
				q_update[s].pinn *= dt;
			#endif
			#endif
			#endif

			#ifdef WTZMU
				q_update[s].WtTz *= dt;
				q_update[s].WxTz *= dt;
				q_update[s].WyTz *= dt;
				q_update[s].WnTz *= dt;
			#endif

			#ifdef PI
				q_update[s].Pi *= dt;
			#endif
			}
		}
	}
}

void evolve_hydro_one_time_step_2_stages(int n, precision t, precision dt, precision dt_prev, lattice_parameters lattice, hydro_parameters hydro, int stage, bool hit_CFL)
{
	int end_RK = 0;
	precision t0 = t;																		// time at start of update

	if(lattice.adaptive_time_step && !hit_CFL)                                      		// q1 <-- dt.E(t,q)
	{
		if(n == 0)
		{
			euler_step(t, q, q1, e, lambda, aT, aL, up, u, dt, dt_prev, lattice, hydro, stage);
		}
		else
		{
			recompute_euler_step(q1, dt, lattice);
		}
	}
	else
	{
		euler_step(t, q, q1, e, lambda, aT, aL, up, u, dt, dt_prev, lattice, hydro, stage);
	}

	get_RK_stage_1(q, q1, A11, lattice);													// q1 <-- q + a11.q1

	t = t0 + C1*dt;																			// set time, previous time step for second stage
	dt_prev = t - t0;

	swap_fluid_velocity(&u, &up);                                                   		// swap u <--> up before second stage

#ifdef ANISO_HYDRO
	set_inferred_variables_aniso_hydro(t, q1, e, u, lattice, hydro);                		// compute aniso (e, u)
	// todo: try regulating before solving for aniso variables (haven't done yet)
#ifdef LATTICE_QCD
	set_anisotropic_variables(q1, e, lambda, aT, aL, lattice, hydro);               		// compute aniso X = (lambda, aT, aL)
#endif
	regulate_residual_currents(t, q1, e, u, lattice, hydro, end_RK);                   		// regulate aniso q1
#else
	set_inferred_variables_viscous_hydro(t, q1, e, u, lattice, hydro);              		// compute viscous (e, u)
	regulate_viscous_currents(t, q1, e, u, lattice, hydro, end_RK);                    		// regulate viscous q1
#endif

	set_ghost_cells(q1, e, u, lattice);                                             		// for (q1, e, u)

	euler_step(t, q1, q2, e, lambda, aT, aL, up, u, dt, dt_prev, lattice, hydro, stage);	// q2 <-- dt.E(t + c1.dt, q1)

	get_RK_stage_2(q, q1, q2, B0, B1, B2, lattice);											// q2 <-- b0.q + b1.q1 + b2.q2

	end_RK = 1;																				// finished RK update for q2
	t = t0 + dt;																			// time of updated step

#ifdef ANISO_HYDRO
	set_inferred_variables_aniso_hydro(t, q2, e, u, lattice, hydro);                		// compute aniso (e, u)
#ifdef LATTICE_QCD
	set_anisotropic_variables(q2, e, lambda, aT, aL, lattice, hydro);               		// compute aniso X = (lambda, aT, aL)
#endif
	regulate_residual_currents(t, q2, e, u, lattice, hydro, end_RK);                   		// regulate aniso q2
#else
	set_inferred_variables_viscous_hydro(t, q2, e, u, lattice, hydro);              		// compute viscous (e, u)
	regulate_viscous_currents(t, q2, e, u, lattice, hydro, end_RK);                    		// regulate viscous q2
#endif

	swap_hydro_variables(&q, &q2);                                                  		// swap q <--> q2
	set_ghost_cells(q, e, u, lattice);                                              		// for (q, e, u)
}


void evolve_hydro_one_time_step_3_stages(int n, precision t, precision dt, precision dt_prev, lattice_parameters lattice, hydro_parameters hydro, int stage, bool hit_CFL)
{
#if (STAGES == 3)
	int end_RK = 0;
	precision t0 = t;

	if(lattice.adaptive_time_step && !hit_CFL)                                      		// q1 <-- dt.E(t,q)
	{
		if(n == 0)
		{
			euler_step(t, q, q1, e, lambda, aT, aL, up, u, dt, dt_prev, lattice, hydro, stage);
		}
		else
		{
			recompute_euler_step(q1, dt, lattice);
		}
	}
	else
	{
		euler_step(t, q, q1, e, lambda, aT, aL, up, u, dt, dt_prev, lattice, hydro, stage);
	}

	get_RK_stage_1(q, q1, A11, lattice);													// q1 <-- q + a11.q1

	t = t0 + C1*dt;																			// set time, previous time step for second stage
	dt_prev = t - t0;                                                                   		

#ifdef ANISO_HYDRO
	set_inferred_variables_aniso_hydro(t, q1, e, uI, lattice, hydro);                		// compute aniso (e, uI)
#ifdef LATTICE_QCD
	set_anisotropic_variables(q1, e, lambda, aT, aL, lattice, hydro);               		// compute aniso X = (lambda, aT, aL)
#endif
	regulate_residual_currents(t, q1, e, uI, lattice, hydro, end_RK);                   	// regulate aniso q1
#else
	set_inferred_variables_viscous_hydro(t, q1, e, uI, lattice, hydro);             		// compute viscous (e, uI)
	regulate_viscous_currents(t, q1, e, uI, lattice, hydro, end_RK);                    	// regulate viscous q1
#endif

	set_ghost_cells(q1, e, uI, lattice);                                             		// for (q1, e, uI)

	euler_step(t, q1, q2, e, lambda, aT, aL, u, uI, dt, dt_prev, lattice, hydro, stage);	// q2 <-- dt.E(t + c1.dt, q1)

	get_RK_stage_2(q, q1, q2, A20, A21, A22, lattice);										// q2 <-- a20.q + a21.q1 + a22.q2

	t = t0 + C2*dt;																			// set time, previous time step for third stage
	dt_prev = t - t0;																		

#ifdef ANISO_HYDRO
	set_inferred_variables_aniso_hydro(t, q2, e, uI, lattice, hydro);                		// compute aniso (e, uI)
#ifdef LATTICE_QCD
	set_anisotropic_variables(q2, e, lambda, aT, aL, lattice, hydro);               		// compute aniso X = (lambda, aT, aL)
#endif
	regulate_residual_currents(t, q2, e, uI, lattice, hydro, end_RK);                   	// regulate aniso q2
#else
	set_inferred_variables_viscous_hydro(t, q2, e, uI, lattice, hydro);              		// compute viscous (e, uI)
	regulate_viscous_currents(t, q2, e, uI, lattice, hydro, end_RK);                    	// regulate viscous q2
#endif

	set_ghost_cells(q2, e, uI, lattice);                                            		// for (q2, e, uI)

	euler_step(t, q2, q3, e, lambda, aT, aL, u, uI, dt, dt_prev, lattice, hydro, stage);	// q3 <-- dt.E(t + c2.dt, q2)							

	get_RK_stage_3(q, q1, q2, q3, B0, B1, B2, B3, lattice);									// q3 <-- b0.q + b1.q1 + b2.q2 + b3.q3

	end_RK = 1;																				// finished RK update for q3
	t = t0 + dt;																			// time of updated step

	swap_fluid_velocity(&u, &up);                                                   		// swap u <--> up after RK update

#ifdef ANISO_HYDRO
	set_inferred_variables_aniso_hydro(t, q3, e, u, lattice, hydro);                		// compute aniso (e, u)
#ifdef LATTICE_QCD
	set_anisotropic_variables(q3, e, lambda, aT, aL, lattice, hydro);               		// compute aniso X = (lambda, aT, aL)
#endif
	regulate_residual_currents(t, q3, e, u, lattice, hydro, end_RK);                   		// regulate aniso q3
#else
	set_inferred_variables_viscous_hydro(t, q3, e, u, lattice, hydro);              		// compute viscous (e, u)
	regulate_viscous_currents(t, q3, e, u, lattice, hydro, end_RK);                    		// regulate viscous q3
#endif

	swap_hydro_variables(&q, &q3);                                                  		// swap q <--> q3
	set_ghost_cells(q, e, u, lattice);                                              		// for (q, e, u)

#endif
}


void evolve_hydro_one_time_step_4_stages(int n, precision t, precision dt, precision dt_prev, lattice_parameters lattice, hydro_parameters hydro, int stage, bool hit_CFL)
{
#if (STAGES == 4)
	int end_RK = 0;
	precision t0 = t;

	if(lattice.adaptive_time_step && !hit_CFL)                                      		// q1 <-- dt.E(t,q)
	{
		if(n == 0)
		{
			euler_step(t, q, q1, e, lambda, aT, aL, up, u, dt, dt_prev, lattice, hydro, stage);
		}
		else
		{
			recompute_euler_step(q1, dt, lattice);
		}
	}
	else
	{
		euler_step(t, q, q1, e, lambda, aT, aL, up, u, dt, dt_prev, lattice, hydro, stage);
	}

	get_RK_stage_1(q, q1, A11, lattice);													// q1 <-- q + a11.q1

	t = t0 + C1*dt;																			// set time, previous time step for second stage
	dt_prev = t - t0;                                                                   		

#ifdef ANISO_HYDRO
	set_inferred_variables_aniso_hydro(t, q1, e, uI, lattice, hydro);                		// compute aniso (e, uI)
#ifdef LATTICE_QCD
	set_anisotropic_variables(q1, e, lambda, aT, aL, lattice, hydro);               		// compute aniso X = (lambda, aT, aL)
#endif
	regulate_residual_currents(t, q1, e, uI, lattice, hydro, end_RK);                   	// regulate aniso q1
#else
	set_inferred_variables_viscous_hydro(t, q1, e, uI, lattice, hydro);             		// compute viscous (e, uI)
	regulate_viscous_currents(t, q1, e, uI, lattice, hydro, end_RK);                    	// regulate viscous q1
#endif

	set_ghost_cells(q1, e, uI, lattice);                                             		// for (q1, e, uI)

	euler_step(t, q1, q2, e, lambda, aT, aL, u, uI, dt, dt_prev, lattice, hydro, stage);	// q2 <-- dt.E(t + c1.dt, q1)

	get_RK_stage_2(q, q1, q2, A20, A21, A22, lattice);										// q2 <-- a20.q + a21.q1 + a22.q2

	t = t0 + C2*dt;																			// set time, previous time step for third stage
	dt_prev = t - t0;																		

#ifdef ANISO_HYDRO
	set_inferred_variables_aniso_hydro(t, q2, e, uI, lattice, hydro);                		// compute aniso (e, uI)
#ifdef LATTICE_QCD
	set_anisotropic_variables(q2, e, lambda, aT, aL, lattice, hydro);               		// compute aniso X = (lambda, aT, aL)
#endif
	regulate_residual_currents(t, q2, e, uI, lattice, hydro, end_RK);                   	// regulate aniso q2
#else
	set_inferred_variables_viscous_hydro(t, q2, e, uI, lattice, hydro);              		// compute viscous (e, uI)
	regulate_viscous_currents(t, q2, e, uI, lattice, hydro, end_RK);                    	// regulate viscous q2
#endif

	set_ghost_cells(q2, e, uI, lattice);                                            		// for (q2, e, uI)

	euler_step(t, q2, q3, e, lambda, aT, aL, u, uI, dt, dt_prev, lattice, hydro, stage);	// q3 <-- dt.E(t + c2.dt, q2)												
	get_RK_stage_3(q, q1, q2, q3, A30, A31, A32, A33, lattice);								// q3 <-- a30.q + a31.q1 + a32.q2 + a33.q3

	t = t0 + C3*dt;																			// set time, previous time step for fourth stage
	dt_prev = t - t0;	

#ifdef ANISO_HYDRO
	set_inferred_variables_aniso_hydro(t, q3, e, uI, lattice, hydro);                		// compute aniso (e, uI)
#ifdef LATTICE_QCD
	set_anisotropic_variables(q3, e, lambda, aT, aL, lattice, hydro);               		// compute aniso X = (lambda, aT, aL)
#endif
	regulate_residual_currents(t, q3, e, uI, lattice, hydro, end_RK);                   	// regulate aniso q3
#else
	set_inferred_variables_viscous_hydro(t, q3, e, uI, lattice, hydro);              		// compute viscous (e, uI)
	regulate_viscous_currents(t, q3, e, uI, lattice, hydro, end_RK);                    	// regulate viscous q3
#endif

	set_ghost_cells(q3, e, uI, lattice);                                            		// for (q3, e, uI)

	euler_step(t, q3, q4, e, lambda, aT, aL, u, uI, dt, dt_prev, lattice, hydro, stage);	// q4 <-- dt.E(t + c3.dt, q3)												
	get_RK_stage_4(q, q1, q2, q3, q4, B0, B1, B2, B3, B4, lattice);							// q4 <-- b0.q + b1.q1 + b2.q2 + b3.q3 + b4.q4

	end_RK = 1;																				// finished RK update for q4
	t = t0 + dt;																			// time of updated step

	swap_fluid_velocity(&u, &up);                                                   		// swap u <--> up after RK update

#ifdef ANISO_HYDRO
	set_inferred_variables_aniso_hydro(t, q4, e, u, lattice, hydro);                		// compute aniso (e, u)
#ifdef LATTICE_QCD
	set_anisotropic_variables(q4, e, lambda, aT, aL, lattice, hydro);               		// compute aniso X = (lambda, aT, aL)
#endif
	regulate_residual_currents(t, q4, e, u, lattice, hydro, end_RK);                   		// regulate aniso q4
#else
	set_inferred_variables_viscous_hydro(t, q4, e, u, lattice, hydro);              		// compute viscous (e, u)
	regulate_viscous_currents(t, q4, e, u, lattice, hydro, end_RK);                    		// regulate viscous q4
#endif

	swap_hydro_variables(&q, &q4);                                                  		// swap q <--> q4
	set_ghost_cells(q, e, u, lattice);                                              		// for (q, e, u)
#endif
}


