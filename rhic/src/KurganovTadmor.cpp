#include <stdlib.h>
#include <math.h>
#include "../include/Macros.h"
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"
#include "../include/GhostCells.h"
#include "../include/FluxTerms.h"
#include "../include/SourceTerms.h"
#include "../include/InferredVariables.h"
#include "../include/NeighborCells.h"
#include "../include/Regulation.h"

inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


void euler_step(precision t, const hydro_variables * const __restrict__ q_current, hydro_variables * const __restrict__ q_update,
const precision * const __restrict__ e_current, precision * const __restrict__ e_update, const fluid_velocity * const __restrict__ u_previous, const fluid_velocity * const __restrict__ u_current, fluid_velocity * const __restrict__ u_update, precision dt, precision dt_prev, lattice_parameters lattice, hydro_parameters hydro, int update, int RK2)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision dx = lattice.lattice_spacing_x;
	precision dy = lattice.lattice_spacing_y;
	precision dn = lattice.lattice_spacing_eta;

	precision Theta = hydro.flux_limiter;

	precision t2 = t * t;							// current time
	precision t4 = t2 * t2;

	precision t_update = t + dt;					// time after first intermediate Euler step
	precision t2_update = t_update * t_update;
	precision t4_update = t2_update * t2_update;

	int stride_y = nx + 4;							// strides for neighbor cells along x, y, n (stride_x = 1)
	int stride_z = (nx + 4) * (ny + 4);				// stride formulas based from linear_column_index()

	precision qs[NUMBER_CONSERVED_VARIABLES];		// current variables at cell s
	precision f[NUMBER_CONSERVED_VARIABLES];		// source function at cell s
	precision  S[NUMBER_CONSERVED_VARIABLES];		// external source terms in hydrodynamic equations

	precision e1[6];								// primary variables of neighbors [i-1, i+1, j-1, j+1, k-1, k+1]

	precision qi1[2 * NUMBER_CONSERVED_VARIABLES];	// conserved variables of neighbor cells along x [i-1, i+1]
	precision qj1[2 * NUMBER_CONSERVED_VARIABLES];	// conserved variables of neighbor cells along y [j-1, j+1]
	precision qk1[2 * NUMBER_CONSERVED_VARIABLES];	// conserved variables of neighbor cells along n [k-1, k+1]

	precision qi2[2 * NUMBER_CONSERVED_VARIABLES];	// conserved variables of neighbor cells along x [i-2, i+2]
	precision qj2[2 * NUMBER_CONSERVED_VARIABLES];	// conserved variables of neighbor cells along y [j-2, j+2]
	precision qk2[2 * NUMBER_CONSERVED_VARIABLES];	// conserved variables of neighbor cells along n [k-2, k+2]

	precision ui1[6];								// fluid velocity of neighbor cells along x [i-1, i+1]
	precision uj1[6];								// fluid velocity of neighbor cells along y [j-1, j+1]
	precision uk1[6];								// fluid velocity of neighbor cells along n [k-1, k+1]

	precision vxi[4];								// vx of neighbor cells along x [i-2, i-1, i+1, i+2]
	precision vyj[4];								// vy of neighbor cells along y [j-2, j-1, j+1, j+2]
	precision vnk[4];								// vn of neighbor cells along n [k-2, k-1, k+1, k+2]

	precision  Hx_plus[NUMBER_CONSERVED_VARIABLES];	// Hx_{i + 1/2}
	precision Hx_minus[NUMBER_CONSERVED_VARIABLES];	// Hx_{i - 1/2}

	precision  Hy_plus[NUMBER_CONSERVED_VARIABLES];	// Hy_{j + 1/2}
	precision Hy_minus[NUMBER_CONSERVED_VARIABLES];	// Hy_{j - 1/2}

	precision  Hn_plus[NUMBER_CONSERVED_VARIABLES];	// Hn_{k + 1/2}
	precision Hn_minus[NUMBER_CONSERVED_VARIABLES];	// Hn_{k - 1/2}

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

				qs[0] = q_current[s].ttt;	// conserved variables of cell s
				qs[1] = q_current[s].ttx;
				qs[2] = q_current[s].tty;
				qs[3] = q_current[s].ttn;

				int a = 4;					// counter

			#ifdef ANISO_HYDRO
				qs[a] = q_current[s].pl; a++;
			#if (PT_MATCHING == 1)
				qs[a] = q_current[s].pt; a++;
			#endif
			#endif
			#ifdef PIMUNU
				qs[a] = q_current[s].pitt; a++;
				qs[a] = q_current[s].pitx; a++;
				qs[a] = q_current[s].pity; a++;
				qs[a] = q_current[s].pitn; a++;
				qs[a] = q_current[s].pixx; a++;
				qs[a] = q_current[s].pixy; a++;
				qs[a] = q_current[s].pixn; a++;
				qs[a] = q_current[s].piyy; a++;
				qs[a] = q_current[s].piyn; a++;
				qs[a] = q_current[s].pinn; a++;
			#endif
			#ifdef WTZMU
				qs[a] = q_current[s].WtTz; a++;
				qs[a] = q_current[s].WxTz; a++;
				qs[a] = q_current[s].WyTz; a++;
				qs[a] = q_current[s].WnTz;
			#endif
			#ifdef PI
				qs[a] = q_current[s].Pi;
			#endif

				precision e_s = e_current[s];

				precision ux = u_current[s].ux;		// current fluid velocity
				precision uy = u_current[s].uy;
				precision un = u_current[s].un;
				precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t2 * un * un);

				precision ux_p = u_previous[s].ux;	// previous fluid velocity
				precision uy_p = u_previous[s].uy;
				precision un_p = u_previous[s].un;

			#if (PT_MATCHING == 0)
				get_primary_neighbor_cells(e_current, e1, sim, sip, sjm, sjp, skm, skp);
			#endif

				get_fluid_velocity_neighbor_cells(u_current[simm], u_current[sim], u_current[sip], u_current[sipp], u_current[sjmm], u_current[sjm], u_current[sjp], u_current[sjpp], u_current[skmm], u_current[skm], u_current[skp], u_current[skpp], ui1, uj1, uk1, vxi, vyj, vnk, t2);

				get_hydro_neighbor_cells(q_current[sim],  q_current[sip],  qi1);
				get_hydro_neighbor_cells(q_current[simm], q_current[sipp], qi2);

				get_hydro_neighbor_cells(q_current[sjm],  q_current[sjp],  qj1);
				get_hydro_neighbor_cells(q_current[sjmm], q_current[sjpp], qj2);

				get_hydro_neighbor_cells(q_current[skm],  q_current[skp],  qk1);
				get_hydro_neighbor_cells(q_current[skmm], q_current[skpp], qk2);


				// compute external source term (S)
			#ifdef ANISO_HYDRO
				source_terms_aniso_hydro(S, qs, e_s, t, qi1, qj1, qk1, e1, ui1, uj1, uk1, ux, uy, un, ux_p, uy_p, un_p, dt_prev, dx, dy, dn, hydro);
			#else
				source_terms_viscous_hydro(S, qs, e_s, t, qi1, qj1, qk1, e1, ui1, uj1, uk1, ux, uy, un, ux_p, uy_p, un_p, dt_prev, dx, dy, dn, hydro);
			#endif

				// compute flux terms (H)
				flux_terms(Hx_plus, Hx_minus, qs, qi1, qi2, vxi, ux / ut, Theta);
				flux_terms(Hy_plus, Hy_minus, qs, qj1, qj2, vyj, uy / ut, Theta);
				flux_terms(Hn_plus, Hn_minus, qs, qk1, qk2, vnk, un / ut, Theta);


				if(!update)
				{
					for(int n = 0; n < NUMBER_CONSERVED_VARIABLES; n++)
					{
						f[n] = S[n]  +  (Hx_minus[n] - Hx_plus[n]) / dx  +  (Hy_minus[n] - Hy_plus[n]) / dy  +  (Hn_minus[n] - Hn_plus[n]) / dn;
					}

					q_update[s].ttt = f[0];
					q_update[s].ttx = f[1];
					q_update[s].tty = f[2];
					q_update[s].ttn = f[3];

					a = 4;
				#ifdef ANISO_HYDRO
					q_update[s].pl = f[a]; a++;
				#if (PT_MATCHING == 1)
					q_update[s].pt = f[a]; a++;
				#endif
				#endif
				#ifdef PIMUNU
					q_update[s].pitt = f[a]; a++;
					q_update[s].pitx = f[a]; a++;
					q_update[s].pity = f[a]; a++;
					q_update[s].pitn = f[a]; a++;
					q_update[s].pixx = f[a]; a++;
					q_update[s].pixy = f[a]; a++;
					q_update[s].pixn = f[a]; a++;
					q_update[s].piyy = f[a]; a++;
					q_update[s].piyn = f[a]; a++;
					q_update[s].pinn = f[a]; a++;
				#endif
				#ifdef WTZMU
					q_update[s].WtTz = f[a]; a++;
					q_update[s].WxTz = f[a]; a++;
					q_update[s].WyTz = f[a]; a++;
					q_update[s].WnTz = f[a];
				#endif
				#ifdef PI
					q_update[s].Pi = f[a];
				#endif
				}
				else if(!RK2)
				{
					for(int n = 0; n < NUMBER_CONSERVED_VARIABLES; n++)
					{
						qs[n] += dt * (S[n]  +  (Hx_minus[n] - Hx_plus[n]) / dx  +  (Hy_minus[n] - Hy_plus[n]) / dy  +  (Hn_minus[n] - Hn_plus[n]) / dn);
					}

					hydro_variables q_intermediate;

					q_intermediate.ttt = qs[0];
					q_intermediate.ttx = qs[1];
					q_intermediate.tty = qs[2];
					q_intermediate.ttn = qs[3];

					a = 4;
				#ifdef ANISO_HYDRO
					q_intermediate.pl = qs[a]; a++;
				#if (PT_MATCHING == 1)
					q_intermediate.pt = qs[a]; a++;
				#endif
				#endif
				#ifdef PIMUNU
					q_intermediate.pitt = qs[a]; a++;
					q_intermediate.pitx = qs[a]; a++;
					q_intermediate.pity = qs[a]; a++;
					q_intermediate.pitn = qs[a]; a++;
					q_intermediate.pixx = qs[a]; a++;
					q_intermediate.pixy = qs[a]; a++;
					q_intermediate.pixn = qs[a]; a++;
					q_intermediate.piyy = qs[a]; a++;
					q_intermediate.piyn = qs[a]; a++;
					q_intermediate.pinn = qs[a]; a++;
				#endif
				#ifdef WTZMU
					q_intermediate.WtTz = qs[a]; a++;
					q_intermediate.WxTz = qs[a]; a++;
					q_intermediate.WyTz = qs[a]; a++;
					q_intermediate.WnTz = qs[a];
				#endif
				#ifdef PI
					q_intermediate.Pi = qs[a];
				#endif

					inferred_variables root = solve_inferred_variables_aniso_hydro(q_intermediate, t_update, t2_update, t4_update, hydro);

					e_update[s] = root.energy;
					u_update[s].ux = root.ux;
					u_update[s].uy = root.uy;
					u_update[s].un = root.un;

					q_intermediate = regulate_viscous_currents_aniso(q_intermediate, root, t_update, t2_update, t4_update, hydro);

					q_update[s] = q_intermediate;
				}
				else
				{
					for(int n = 0; n < NUMBER_CONSERVED_VARIABLES; n++)
					{
						qs[n] += dt * (S[n]  +  (Hx_minus[n] - Hx_plus[n]) / dx  +  (Hy_minus[n] - Hy_plus[n]) / dy  +  (Hn_minus[n] - Hn_plus[n]) / dn);
					}

					hydro_variables q_start = q[s];
					hydro_variables q_final;

					q_final.ttt = (q_start.ttt  +  qs[0]) / 2.;
					q_final.ttx = (q_start.ttx  +  qs[1]) / 2.;
					q_final.tty = (q_start.tty  +  qs[2]) / 2.;
					q_final.ttn = (q_start.ttn  +  qs[3]) / 2.;

					a = 4;
				#ifdef ANISO_HYDRO
					q_final.pl = (q_start.pl  +  qs[a]) / 2.; a++;
				#if (PT_MATCHING == 1)
					q_final.pt = (q_start.pt  +  qs[a]) / 2.; a++;
				#endif
				#endif
				#ifdef PIMUNU
					q_final.pitt = (q_start.pitt  +  qs[a]) / 2.; a++;
					q_final.pitx = (q_start.pitx  +  qs[a]) / 2.; a++;
					q_final.pity = (q_start.pity  +  qs[a]) / 2.; a++;
					q_final.pitn = (q_start.pitn  +  qs[a]) / 2.; a++;
					q_final.pixx = (q_start.pixx  +  qs[a]) / 2.; a++;
					q_final.pixy = (q_start.pixy  +  qs[a]) / 2.; a++;
					q_final.pixn = (q_start.pixn  +  qs[a]) / 2.; a++;
					q_final.piyy = (q_start.piyy  +  qs[a]) / 2.; a++;
					q_final.piyn = (q_start.piyn  +  qs[a]) / 2.; a++;
					q_final.pinn = (q_start.pinn  +  qs[a]) / 2.; a++;
				#endif
				#ifdef WTZMU
					q_final.WtTz = (q_start.WtTz  +  qs[a]) / 2.; a++;
					q_final.WxTz = (q_start.WxTz  +  qs[a]) / 2.; a++;
					q_final.WyTz = (q_start.WyTz  +  qs[a]) / 2.; a++;
					q_final.WnTz = (q_start.WnTz  +  qs[a]) / 2.;
				#endif
				#ifdef PI
					q_final.Pi = (q_start.Pi  +  qs[a]) / 2.;
				#endif

					inferred_variables root = solve_inferred_variables_aniso_hydro(q_final, t, t2, t4, hydro);

					e_update[s] = root.energy;
					u_update[s].ux = root.ux;
					u_update[s].uy = root.uy;
					u_update[s].un = root.un;

					q_final = regulate_viscous_currents_aniso(q_final, root, t, t2, t4, hydro);

					q_update[s] = q_final;
				}
			}
		}
	}
}


void add_euler_step(precision t, const hydro_variables * const __restrict__ q_current, hydro_variables * const __restrict__ q_update, precision * const __restrict__ e_update, fluid_velocity * const __restrict__ u_update, precision dt, lattice_parameters lattice, hydro_parameters hydro)
 {
 	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision t2 = t * t;
	precision t4 = t2 * t2;

	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				hydro_variables f = q_update[s];

				hydro_variables q_intermediate;

				q_intermediate.ttt = q_current[s].ttt  +  dt * f.ttt;		// before this update: Q stores the source function f(q,t)
				q_intermediate.ttx = q_current[s].ttx  +  dt * f.ttx;		// i.e. Q = q + dt.f(q,t)
				q_intermediate.tty = q_current[s].tty  +  dt * f.tty;
				q_intermediate.ttn = q_current[s].ttn  +  dt * f.ttn;

			#ifdef ANISO_HYDRO
				q_intermediate.pl  = q_current[s].pl  +  dt * f.pl;
			#if (PT_MATCHING == 1)
				q_intermediate.pt  = q_current[s].pt  +  dt * f.pt;
			#endif
			#endif
			#ifdef PIMUNU
				q_intermediate.pitt = q_current[s].pitt  +  dt * f.pitt;
				q_intermediate.pitx = q_current[s].pitx  +  dt * f.pitx;
				q_intermediate.pity = q_current[s].pity  +  dt * f.pity;
				q_intermediate.pitn = q_current[s].pitn  +  dt * f.pitn;
				q_intermediate.pixx = q_current[s].pixx  +  dt * f.pixx;
				q_intermediate.pixy = q_current[s].pixy  +  dt * f.pixy;
				q_intermediate.pixn = q_current[s].pixn  +  dt * f.pixn;
				q_intermediate.piyy = q_current[s].piyy  +  dt * f.piyy;
				q_intermediate.piyn = q_current[s].piyn  +  dt * f.piyn;
				q_intermediate.pinn = q_current[s].pinn  +  dt * f.pinn;
			#endif
			#ifdef WTZMU
				q_intermediate.WtTz = q_current[s].WtTz  +  dt * f.WtTz;
				q_intermediate.WxTz = q_current[s].WxTz  +  dt * f.WxTz;
				q_intermediate.WyTz = q_current[s].WyTz  +  dt * f.WyTz;
				q_intermediate.WnTz = q_current[s].WnTz  +  dt * f.WnTz;
			#endif
			#ifdef PI
				q_intermediate.Pi = q_current[s].Pi  +  dt * f.Pi;
			#endif

				inferred_variables root = solve_inferred_variables_aniso_hydro(q_intermediate, t, t2, t4, hydro);

				e_update[s] = root.energy;
				u_update[s].ux = root.ux;
				u_update[s].uy = root.uy;
				u_update[s].un = root.un;

				q_intermediate = regulate_viscous_currents_aniso(q_intermediate, root, t, t2, t4, hydro);

				q_update[s] = q_intermediate;
			}
		}
	}
}


void evolve_hydro_one_time_step(int n, precision t, precision dt, precision dt_prev, lattice_parameters lattice, hydro_parameters hydro, int update)
{
	int RK2 = 0;

	if(lattice.adaptive_time_step)											// first intermediate euler step
	{																		// qI = q + dt.(S - dHx/dx - dHy/dy - dHn/dn)
		if(n == 0)
		{
			euler_step(t, q, qI, e, E, up, u, uI, dt, dt_prev, lattice, hydro, update, RK2);
		}
		else
		{
			add_euler_step(t + dt, q, qI, E, uI, dt, lattice, hydro);
		}
	}
	else
	{
		euler_step(t, q, qI, e, E, up, u, uI, dt, dt_prev, lattice, hydro, update, RK2);
	}

	t += dt;																// next intermediate time step

	// next I want to move inferred variables


#ifdef ANISO_HYDRO
	//set_inferred_variables_aniso_hydro(qI, E, uI, t, lattice, hydro);		// compute (uI, e)
	//regulate_residual_currents(t, qI, E, uI, lattice, hydro);				// regulate qI
#else
	//set_inferred_variables_viscous_hydro(qI, E, uI, t, lattice, hydro);
	//regulate_viscous_currents(t, qI, E, uI, lattice, hydro);
#endif

	set_ghost_cells(qI, E, uI, lattice);									// set (qI, uI, e) ghost cells

	dt_prev = dt;															// update previous time step

	RK2 = 1;

	euler_step(t, qI, Q, E, e, u, uI, up, dt, dt_prev, lattice, hydro, update, RK2);	// second intermediate euler step
																			// Q = qI + dt.(S - dHx/dx - dHy/dy - dHn/dn)


#ifdef ANISO_HYDRO
	//set_inferred_variables_aniso_hydro(Q, e, up, t, lattice, hydro);			// compute (u, e)
	//regulate_residual_currents(t, Q, e, up, lattice, hydro);					// regulate Q
#else
	//set_inferred_variables_viscous_hydro(Q, e, up, t, lattice, hydro);
	//regulate_viscous_currents(t, Q, e, up, lattice, hydro);
#endif

	swap_hydro_variables(&q, &Q);											// swap q and Q
	swap_fluid_velocity(&u, &up);											// swap u and up

	set_ghost_cells(q, e, u, lattice);										// set (Q, u, e) ghost cells


}







