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


void euler_step(precision t, const hydro_variables * const __restrict__ Q_current, hydro_variables * const __restrict__ Q_update,
const precision * const __restrict__ e, const fluid_velocity * const __restrict__ u, const fluid_velocity * const __restrict__ up, precision dt, precision dt_prev, lattice_parameters lattice, hydro_parameters hydro, int update)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	precision dx = lattice.lattice_spacing_x;
	precision dy = lattice.lattice_spacing_y;
	precision dn = lattice.lattice_spacing_eta;

	precision Theta = hydro.flux_limiter;

	precision t2 = t * t;
	int stride_y = nx + 4;							// strides for neighbor cells along x, y, n (stride_x = 1)
	int stride_z = (nx + 4) * (ny + 4);				// stride formulas based from linear_column_index()

	precision Qs[NUMBER_CONSERVED_VARIABLES];		// current variables at cell s
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

				Qs[0] = Q_current[s].ttt;	// conserved variables of cell s
				Qs[1] = Q_current[s].ttx;
				Qs[2] = Q_current[s].tty;
				Qs[3] = Q_current[s].ttn;

				int a = 4;					// counter

			#ifdef ANISO_HYDRO
				Qs[a] = Q_current[s].pl; a++;
			#if (PT_MATCHING == 1)
				Qs[a] = Q_current[s].pt; a++;
			#endif
			#endif
			#ifdef PIMUNU
				Qs[a] = Q_current[s].pitt; a++;
				Qs[a] = Q_current[s].pitx; a++;
				Qs[a] = Q_current[s].pity; a++;
				Qs[a] = Q_current[s].pitn; a++;
				Qs[a] = Q_current[s].pixx; a++;
				Qs[a] = Q_current[s].pixy; a++;
				Qs[a] = Q_current[s].pixn; a++;
				Qs[a] = Q_current[s].piyy; a++;
				Qs[a] = Q_current[s].piyn; a++;
				Qs[a] = Q_current[s].pinn; a++;
			#endif
			#ifdef WTZMU
				Qs[a] = Q_current[s].WtTz; a++;
				Qs[a] = Q_current[s].WxTz; a++;
				Qs[a] = Q_current[s].WyTz; a++;
				Qs[a] = Q_current[s].WnTz;
			#endif
			#ifdef PI
				Qs[a] = Q_current[s].Pi;
			#endif

				precision e_s = e[s];

				precision ux = u[s].ux;		// current fluid velocity 
				precision uy = u[s].uy;
				precision un = u[s].un;
				precision ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t2 * un * un);

				precision ux_p = up[s].ux;	// previous fluid velocity
				precision uy_p = up[s].uy;
				precision un_p = up[s].un;

			#if (PT_MATCHING == 0)
				get_primary_neighbor_cells(e, e1, sim, sip, sjm, sjp, skm, skp);
			#endif

				get_fluid_velocity_neighbor_cells(u[simm], u[sim], u[sip], u[sipp], u[sjmm], u[sjm], u[sjp], u[sjpp], u[skmm], u[skm], u[skp], u[skpp], ui1, uj1, uk1, vxi, vyj, vnk, t2);

				get_hydro_neighbor_cells(Q_current[sim],  Q_current[sip],  qi1);
				get_hydro_neighbor_cells(Q_current[simm], Q_current[sipp], qi2);

				get_hydro_neighbor_cells(Q_current[sjm],  Q_current[sjp],  qj1);
				get_hydro_neighbor_cells(Q_current[sjmm], Q_current[sjpp], qj2);

				get_hydro_neighbor_cells(Q_current[skm],  Q_current[skp],  qk1);
				get_hydro_neighbor_cells(Q_current[skmm], Q_current[skpp], qk2);


				// compute external source term (S)
			#ifdef ANISO_HYDRO
				source_terms_aniso_hydro(S, Qs, e_s, t, qi1, qj1, qk1, e1, ui1, uj1, uk1, ux, uy, un, ux_p, uy_p, un_p, dt_prev, dx, dy, dn, hydro);
			#else
				source_terms_viscous_hydro(S, Qs, e_s, t, qi1, qj1, qk1, e1, ui1, uj1, uk1, ux, uy, un, ux_p, uy_p, un_p, dt_prev, dx, dy, dn, hydro);
			#endif

				// compute flux terms (H)
				flux_terms(Hx_plus, Hx_minus, Qs, qi1, qi2, vxi, ux / ut, Theta);
				flux_terms(Hy_plus, Hy_minus, Qs, qj1, qj2, vyj, uy / ut, Theta);
				flux_terms(Hn_plus, Hn_minus, Qs, qk1, qk2, vnk, un / ut, Theta);

				for(int n = 0; n < NUMBER_CONSERVED_VARIABLES; n++)
				{	
					if(update) 	// add euler step to Qs (usual time evolution)
					{
						Qs[n] += dt * (S[n]  +  (Hx_minus[n] - Hx_plus[n]) / dx  +  (Hy_minus[n] - Hy_plus[n]) / dy  +  (Hn_minus[n] - Hn_plus[n]) / dn);
					}
					else 		// only store total source function in Qs (for computing the adaptive time step)
					{
						Qs[n] = S[n]  +  (Hx_minus[n] - Hx_plus[n]) / dx  +  (Hy_minus[n] - Hy_plus[n]) / dy  +  (Hn_minus[n] - Hn_plus[n]) / dn;
					}
				}

				Q_update[s].ttt = Qs[0];		// store the final results
				Q_update[s].ttx = Qs[1];
				Q_update[s].tty = Qs[2];
				Q_update[s].ttn = Qs[3];

				a = 4;							// reset counter

			#ifdef ANISO_HYDRO
				Q_update[s].pl = Qs[a]; a++;
			#if (PT_MATCHING == 1)
				Q_update[s].pt = Qs[a]; a++;
			#endif
			#endif
			#ifdef PIMUNU
				Q_update[s].pitt = Qs[a]; a++;
				Q_update[s].pitx = Qs[a]; a++;
				Q_update[s].pity = Qs[a]; a++;
				Q_update[s].pitn = Qs[a]; a++;
				Q_update[s].pixx = Qs[a]; a++;
				Q_update[s].pixy = Qs[a]; a++;
				Q_update[s].pixn = Qs[a]; a++;
				Q_update[s].piyy = Qs[a]; a++;
				Q_update[s].piyn = Qs[a]; a++;
				Q_update[s].pinn = Qs[a]; a++;
			#endif
			#ifdef WTZMU
				Q_update[s].WtTz = Qs[a]; a++;
				Q_update[s].WxTz = Qs[a]; a++;
				Q_update[s].WyTz = Qs[a]; a++;
				Q_update[s].WnTz = Qs[a];
			#endif
			#ifdef PI
				Q_update[s].Pi = Qs[a];
			#endif
			}
		}
	}
}


void add_euler_step(const hydro_variables * const __restrict__ q, hydro_variables * const __restrict__ Q, precision dt, lattice_parameters lattice)
 {
 	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;
 	
	for(int k = 2; k < nz + 2; k++)			
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				Q[s].ttt = q[s].ttt  +  dt * Q[s].ttt;		// before this update: Q stores the source function f(q,t)
				Q[s].ttx = q[s].ttx  +  dt * Q[s].ttx;		// i.e. Q = q + dt.f(q,t)
				Q[s].tty = q[s].tty  +  dt * Q[s].tty;
				Q[s].ttn = q[s].ttn  +  dt * Q[s].ttn;

			#ifdef ANISO_HYDRO
				Q[s].pl  = q[s].pl  +  dt * Q[s].pl;
			#if (PT_MATCHING == 1)
				Q[s].pt  = q[s].pt  +  dt * Q[s].pt;
			#endif
			#endif
			#ifdef PIMUNU
				Q[s].pitt = q[s].pitt  +  dt * Q[s].pitt;
				Q[s].pitx = q[s].pitx  +  dt * Q[s].pitx;
				Q[s].pity = q[s].pity  +  dt * Q[s].pity;
				Q[s].pitn = q[s].pitn  +  dt * Q[s].pitn;
				Q[s].pixx = q[s].pixx  +  dt * Q[s].pixx;
				Q[s].pixy = q[s].pixy  +  dt * Q[s].pixy;
				Q[s].pixn = q[s].pixn  +  dt * Q[s].pixn;
				Q[s].piyy = q[s].piyy  +  dt * Q[s].piyy;
				Q[s].piyn = q[s].piyn  +  dt * Q[s].piyn;
				Q[s].pinn = q[s].pinn  +  dt * Q[s].pinn;
			#endif
			#ifdef WTZMU
				Q[s].WtTz = q[s].WtTz  +  dt * Q[s].WtTz;
				Q[s].WxTz = q[s].WxTz  +  dt * Q[s].WxTz;
				Q[s].WyTz = q[s].WyTz  +  dt * Q[s].WyTz;
				Q[s].WnTz = q[s].WnTz  +  dt * Q[s].WnTz;
			#endif
			#ifdef PI
				Q[s].Pi = q[s].Pi  +  dt * Q[s].Pi;
			#endif
			}
		}
	}
}


void runge_kutta2(const hydro_variables * const __restrict__ q, hydro_variables * const __restrict__ Q, lattice_parameters lattice)
 {
 	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;
 	
	for(int k = 2; k < nz + 2; k++)			// 2nd order Runge-Kutta update
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				Q[s].ttt = (q[s].ttt  +  Q[s].ttt) / 2.;
				Q[s].ttx = (q[s].ttx  +  Q[s].ttx) / 2.;
				Q[s].tty = (q[s].tty  +  Q[s].tty) / 2.;
				Q[s].ttn = (q[s].ttn  +  Q[s].ttn) / 2.;

			#ifdef ANISO_HYDRO
				Q[s].pl  = (q[s].pl  +  Q[s].pl)  / 2.;
			#if (PT_MATCHING == 1)
				Q[s].pt  = (q[s].pt  +  Q[s].pt)  / 2.;
			#endif
			#endif
			#ifdef PIMUNU
				Q[s].pitt = (q[s].pitt  +  Q[s].pitt) / 2.;
				Q[s].pitx = (q[s].pitx  +  Q[s].pitx) / 2.;
				Q[s].pity = (q[s].pity  +  Q[s].pity) / 2.;
				Q[s].pitn = (q[s].pitn  +  Q[s].pitn) / 2.;
				Q[s].pixx = (q[s].pixx  +  Q[s].pixx) / 2.;
				Q[s].pixy = (q[s].pixy  +  Q[s].pixy) / 2.;
				Q[s].pixn = (q[s].pixn  +  Q[s].pixn) / 2.;
				Q[s].piyy = (q[s].piyy  +  Q[s].piyy) / 2.;
				Q[s].piyn = (q[s].piyn  +  Q[s].piyn) / 2.;
				Q[s].pinn = (q[s].pinn  +  Q[s].pinn) / 2.;
			#endif
			#ifdef WTZMU
				Q[s].WtTz = (q[s].WtTz  +  Q[s].WtTz) / 2.;
				Q[s].WxTz = (q[s].WxTz  +  Q[s].WxTz) / 2.;
				Q[s].WyTz = (q[s].WyTz  +  Q[s].WyTz) / 2.;
				Q[s].WnTz = (q[s].WnTz  +  Q[s].WnTz) / 2.;
			#endif
			#ifdef PI
				Q[s].Pi = (q[s].Pi  +  Q[s].Pi) / 2.;
			#endif
			}
		}
	}
}


void evolve_hydro_one_time_step(int n, precision t, precision dt, precision dt_prev, lattice_parameters lattice, hydro_parameters hydro, int update)
{
	if(lattice.adaptive_time_step)											// first intermediate euler step 
	{																		// qI = q + dt.(S - dHx/dx - dHy/dy - dHn/dn)
		if(n == 0)	
		{
			euler_step(t, q, qI, e, u, up, dt, dt_prev, lattice, hydro, update);
		}
		else
		{
			add_euler_step(q, qI, dt, lattice);																		
		}
	}
	else
	{
		euler_step(t, q, qI, e, u, up, dt, dt_prev, lattice, hydro, update);
	}

	t += dt;																// next intermediate time step

#ifdef ANISO_HYDRO
	set_inferred_variables_aniso_hydro(qI, e, uI, t, lattice, hydro);		// compute (uI, e)
	regulate_residual_currents(t, qI, e, uI, lattice, hydro);				// regulate qI 
#else
	set_inferred_variables_viscous_hydro(qI, e, uI, t, lattice, hydro);
	regulate_viscous_currents(t, qI, e, uI, lattice, hydro);
#endif

	set_ghost_cells(qI, e, uI, lattice);									// set (qI, uI, e) ghost cells

	dt_prev = dt;															// update previous time step 

	euler_step(t, qI, Q, e, uI, u, dt, dt_prev, lattice, hydro, update);	// second intermediate euler step
																			// Q = qI + dt.(S - dHx/dx - dHy/dy - dHn/dn)
	runge_kutta2(q, Q, lattice);											// 2nd order runge kutta: Q = (q + Q) / 2

	swap_fluid_velocity(&up, &u);											// swap up and u

#ifdef ANISO_HYDRO
	set_inferred_variables_aniso_hydro(Q, e, u, t, lattice, hydro);			// compute (u, e)
	regulate_residual_currents(t, Q, e, u, lattice, hydro);					// regulate Q 
#else
	set_inferred_variables_viscous_hydro(Q, e, u, t, lattice, hydro);
	regulate_viscous_currents(t, Q, e, u, lattice, hydro);
#endif

	set_ghost_cells(Q, e, u, lattice);										// set (Q, u, e) ghost cells
	
	swap_hydro_variables(&q, &Q);											// swap q and Q
}







