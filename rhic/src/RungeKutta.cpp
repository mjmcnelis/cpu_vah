
#include <stdlib.h>
#include "../include/Macros.h"
#include "../include/DynamicalVariables.h"
#include "../include/OpenMP.h"


inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


void get_RK_stage_1(const hydro_variables * const __restrict__ q, hydro_variables * const __restrict__ q1, precision a1, lattice_parameters lattice)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	// compute q1, which is the q-argument of the second stage
	// here q1 is passed the first intermediate euler step dt.E(t,q)

	#pragma omp parallel for collapse(3)
	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				q1[s].ttt = q[s].ttt  +  a1 * q1[s].ttt;		// q1 <-- q + a1.q1
				q1[s].ttx = q[s].ttx  +  a1 * q1[s].ttx;
				q1[s].tty = q[s].tty  +  a1 * q1[s].tty;
			#ifndef BOOST_INVARIANT
				q1[s].ttn = q[s].ttn  +  a1 * q1[s].ttn;
			#endif

			#ifdef ANISO_HYDRO
				q1[s].pl = q[s].pl  +  a1 * q1[s].pl;
				q1[s].pt = q[s].pt  +  a1 * q1[s].pt;
			#endif

			#ifdef B_FIELD
				q1[s].b = q[s].b  +  a1 * q1[s].b;
			#endif

			#ifdef PIMUNU
				q1[s].pitt = q[s].pitt  +  a1 * q1[s].pitt;
				q1[s].pitx = q[s].pitx  +  a1 * q1[s].pitx;
				q1[s].pity = q[s].pity  +  a1 * q1[s].pity;
			#ifndef BOOST_INVARIANT
				q1[s].pitn = q[s].pitn  +  a1 * q1[s].pitn;
			#endif
				q1[s].pixx = q[s].pixx  +  a1 * q1[s].pixx;
				q1[s].pixy = q[s].pixy  +  a1 * q1[s].pixy;
			#ifndef BOOST_INVARIANT
				q1[s].pixn = q[s].pixn  +  a1 * q1[s].pixn;
			#endif
				q1[s].piyy = q[s].piyy  +  a1 * q1[s].piyy;
			#ifndef BOOST_INVARIANT
				q1[s].piyn = q[s].piyn  +  a1 * q1[s].piyn;
				q1[s].pinn = q[s].pinn  +  a1 * q1[s].pinn;
			#else
			#ifndef ANISO_HYDRO
				q1[s].pinn = q[s].pinn  +  a1 * q1[s].pinn;
			#endif
			#endif
			#endif

			#ifdef WTZMU
				q1[s].WtTz = q[s].WtTz  +  a1 * q1[s].WtTz;
				q1[s].WxTz = q[s].WxTz  +  a1 * q1[s].WxTz;
				q1[s].WyTz = q[s].WyTz  +  a1 * q1[s].WyTz;
				q1[s].WnTz = q[s].WnTz  +  a1 * q1[s].WnTz;
			#endif

			#ifdef PI
				q1[s].Pi = q[s].Pi  +  a1 * q1[s].Pi;
			#endif
			}
		}
	}
}


void get_RK_stage_2(const hydro_variables * const __restrict__ q, const hydro_variables * const __restrict__ q1, hydro_variables * const __restrict__ q2, precision a0, precision a1, precision a2, lattice_parameters lattice)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	// compute q2, which is either a RK update (s = 2 stages) or the q argument of a third stage
	// here q2 is passed the second intermediate Euler step dt.E(t + c1.dt, q1)

	#pragma omp parallel for collapse(3)
	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);

				// q2 <-- a0.q + a1.q1 + a2.q2

				q2[s].ttt = a0 * q[s].ttt  +  a1 * q1[s].ttt  +  a2 * q2[s].ttt;		
				q2[s].ttx = a0 * q[s].ttx  +  a1 * q1[s].ttx  +  a2 * q2[s].ttx;
				q2[s].tty = a0 * q[s].tty  +  a1 * q1[s].tty  +  a2 * q2[s].tty;
			#ifndef BOOST_INVARIANT
				q2[s].ttn = a0 * q[s].ttn  +  a1 * q1[s].ttn  +  a2 * q2[s].ttn;
			#endif

			#ifdef ANISO_HYDRO
				q2[s].pl = a0 * q[s].pl  +  a1 * q1[s].pl  +  a2 * q2[s].pl;
				q2[s].pt = a0 * q[s].pt  +  a1 * q1[s].pt  +  a2 * q2[s].pt;
			#endif

			#ifdef B_FIELD
				q2[s].b = a0 * q[s].b  +  a1 * q1[s].b  +  a2 * q2[s].b;
			#endif

			#ifdef PIMUNU
				q2[s].pitt = a0 * q[s].pitt  +  a1 * q1[s].pitt  +  a2 * q2[s].pitt;
				q2[s].pitx = a0 * q[s].pitx  +  a1 * q1[s].pitx  +  a2 * q2[s].pitx;
				q2[s].pity = a0 * q[s].pity  +  a1 * q1[s].pity  +  a2 * q2[s].pity;
			#ifndef BOOST_INVARIANT
				q2[s].pitn = a0 * q[s].pitn  +  a1 * q1[s].pitn  +  a2 * q2[s].pitn;
			#endif
				q2[s].pixx = a0 * q[s].pixx  +  a1 * q1[s].pixx  +  a2 * q2[s].pixx;
				q2[s].pixy = a0 * q[s].pixy  +  a1 * q1[s].pixy  +  a2 * q2[s].pixy;
			#ifndef BOOST_INVARIANT
				q2[s].pixn = a0 * q[s].pixn  +  a1 * q1[s].pixn  +  a2 * q2[s].pixn;
			#endif
				q2[s].piyy = a0 * q[s].piyy  +  a1 * q1[s].piyy  +  a2 * q2[s].piyy;
			#ifndef BOOST_INVARIANT
				q2[s].piyn = a0 * q[s].piyn  +  a1 * q1[s].piyn  +  a2 * q2[s].piyn;
				q2[s].pinn = a0 * q[s].pinn  +  a1 * q1[s].pinn  +  a2 * q2[s].pinn;
			#else
			#ifndef ANISO_HYDRO
				q2[s].pinn = a0 * q[s].pinn  +  a1 * q1[s].pinn  +  a2 * q2[s].pinn;
			#endif
			#endif
			#endif

			#ifdef WTZMU
				q2[s].WtTz = a0 * q[s].WtTz  +  a1 * q1[s].WtTz  +  a2 * q2[s].WtTz;
				q2[s].WxTz = a0 * q[s].WxTz  +  a1 * q1[s].WxTz  +  a2 * q2[s].WxTz;
				q2[s].WyTz = a0 * q[s].WyTz  +  a1 * q1[s].WyTz  +  a2 * q2[s].WyTz;
				q2[s].WnTz = a0 * q[s].WnTz  +  a1 * q1[s].WnTz  +  a2 * q2[s].WnTz;
			#endif

			#ifdef PI
				q2[s].Pi = a0 * q[s].Pi  +  a1 * q1[s].Pi  +  a2 * q2[s].Pi;
			#endif
			}
		}
	}
}


void get_RK_stage_3(const hydro_variables * const __restrict__ q, const hydro_variables * const __restrict__ q1, const hydro_variables * const __restrict__ q2, hydro_variables * const __restrict__ q3, precision a0, precision a1, precision a2, precision a3, lattice_parameters lattice)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	// compute q3, which is either a RK update (s = 3 stages) or the q argument of a fourth stage
	// here q3 is passed the second intermediate Euler step dt.E(t + c2.dt, q2)

	#pragma omp parallel for collapse(3)
	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);		

				// q3 <-- a0.q + a1.q1 + a2.q2 + a3.q3

				q3[s].ttt = a0 * q[s].ttt  +  a1 * q1[s].ttt  +  a2 * q2[s].ttt  +  a3 * q3[s].ttt;	
				q3[s].ttx = a0 * q[s].ttx  +  a1 * q1[s].ttx  +  a2 * q2[s].ttx  +  a3 * q3[s].ttx;
				q3[s].tty = a0 * q[s].tty  +  a1 * q1[s].tty  +  a2 * q2[s].tty  +  a3 * q3[s].tty;
			#ifndef BOOST_INVARIANT
				q3[s].ttn = a0 * q[s].ttn  +  a1 * q1[s].ttn  +  a2 * q2[s].ttn  +  a3 * q3[s].ttn;
			#endif

			#ifdef ANISO_HYDRO
				q3[s].pl = a0 * q[s].pl  +  a1 * q1[s].pl  +  a2 * q2[s].pl  +  a3 * q3[s].pl;
				q3[s].pt = a0 * q[s].pt  +  a1 * q1[s].pt  +  a2 * q2[s].pt  +  a3 * q3[s].pt;
			#endif

			#ifdef B_FIELD
				q3[s].b = a0 * q[s].b  +  a1 * q1[s].b  +  a2 * q2[s].b  +  a3 * q3[s].b;
			#endif

			#ifdef PIMUNU
				q3[s].pitt = a0 * q[s].pitt  +  a1 * q1[s].pitt  +  a2 * q2[s].pitt  +  a3 * q3[s].pitt;
				q3[s].pitx = a0 * q[s].pitx  +  a1 * q1[s].pitx  +  a2 * q2[s].pitx  +  a3 * q3[s].pitx;
				q3[s].pity = a0 * q[s].pity  +  a1 * q1[s].pity  +  a2 * q2[s].pity  +  a3 * q3[s].pity;
			#ifndef BOOST_INVARIANT
				q3[s].pitn = a0 * q[s].pitn  +  a1 * q1[s].pitn  +  a2 * q2[s].pitn  +  a3 * q3[s].pitn;
			#endif
				q3[s].pixx = a0 * q[s].pixx  +  a1 * q1[s].pixx  +  a2 * q2[s].pixx  +  a3 * q3[s].pixx;
				q3[s].pixy = a0 * q[s].pixy  +  a1 * q1[s].pixy  +  a2 * q2[s].pixy  +  a3 * q3[s].pixy;
			#ifndef BOOST_INVARIANT
				q3[s].pixn = a0 * q[s].pixn  +  a1 * q1[s].pixn  +  a2 * q2[s].pixn  +  a3 * q3[s].pixn;
			#endif
				q3[s].piyy = a0 * q[s].piyy  +  a1 * q1[s].piyy  +  a2 * q2[s].piyy  +  a3 * q3[s].piyy;
			#ifndef BOOST_INVARIANT
				q3[s].piyn = a0 * q[s].piyn  +  a1 * q1[s].piyn  +  a2 * q2[s].piyn  +  a3 * q3[s].piyn;
				q3[s].pinn = a0 * q[s].pinn  +  a1 * q1[s].pinn  +  a2 * q2[s].pinn  +  a3 * q3[s].pinn;
			#else
			#ifndef ANISO_HYDRO
				q3[s].pinn = a0 * q[s].pinn  +  a1 * q1[s].pinn  +  a2 * q2[s].pinn  +  a3 * q3[s].pinn;
			#endif
			#endif
			#endif

			#ifdef WTZMU
				q3[s].WtTz = a0 * q[s].WtTz  +  a1 * q1[s].WtTz  +  a2 * q2[s].WtTz  +  a3 * q3[s].WtTz;
				q3[s].WxTz = a0 * q[s].WxTz  +  a1 * q1[s].WxTz  +  a2 * q2[s].WxTz  +  a3 * q3[s].WxTz;
				q3[s].WyTz = a0 * q[s].WyTz  +  a1 * q1[s].WyTz  +  a2 * q2[s].WyTz  +  a3 * q3[s].WyTz;
				q3[s].WnTz = a0 * q[s].WnTz  +  a1 * q1[s].WnTz  +  a2 * q2[s].WnTz  +  a3 * q3[s].WnTz;
			#endif

			#ifdef PI
				q3[s].Pi = a0 * q[s].Pi  +  a1 * q1[s].Pi  +  a2 * q2[s].Pi  +  a3 * q3[s].Pi;
			#endif
			}
		}
	}
}


void get_RK_stage_4(const hydro_variables * const __restrict__ q, const hydro_variables * const __restrict__ q1, const hydro_variables * const __restrict__ q2, const hydro_variables * const __restrict__ q3, hydro_variables * const __restrict__ q4, precision a0, precision a1, precision a2, precision a3, precision a4, lattice_parameters lattice)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	// compute q4, which is either a RK update (s = 4 stages) or the q argument of a fifth stage
	// here q4 is passed the second intermediate Euler step dt.E(t + c3.dt, q3)

	#pragma omp parallel for collapse(3)
	for(int k = 2; k < nz + 2; k++)
	{
		for(int j = 2; j < ny + 2; j++)
		{
			for(int i = 2; i < nx + 2; i++)
			{
				int s = linear_column_index(i, j, k, nx + 4, ny + 4);		

				// q4 <-- a0.q + a1.q1 + a2.q2 + a3.q3 + a4.q4

				q4[s].ttt = a0 * q[s].ttt  +  a1 * q1[s].ttt  +  a2 * q2[s].ttt  +  a3 * q3[s].ttt  +  a4 * q4[s].ttt;	
				q4[s].ttx = a0 * q[s].ttx  +  a1 * q1[s].ttx  +  a2 * q2[s].ttx  +  a3 * q3[s].ttx  +  a4 * q4[s].ttx;
				q4[s].tty = a0 * q[s].tty  +  a1 * q1[s].tty  +  a2 * q2[s].tty  +  a3 * q3[s].tty  +  a4 * q4[s].tty;
			#ifndef BOOST_INVARIANT
				q4[s].ttn = a0 * q[s].ttn  +  a1 * q1[s].ttn  +  a2 * q2[s].ttn  +  a3 * q3[s].ttn  +  a4 * q4[s].ttn;
			#endif

			#ifdef ANISO_HYDRO
				q4[s].pl = a0 * q[s].pl  +  a1 * q1[s].pl  +  a2 * q2[s].pl  +  a3 * q3[s].pl  +  a4 * q4[s].pl;
				q4[s].pt = a0 * q[s].pt  +  a1 * q1[s].pt  +  a2 * q2[s].pt  +  a3 * q3[s].pt  +  a4 * q4[s].pt;
			#endif

			#ifdef B_FIELD
				q4[s].b = a0 * q[s].b  +  a1 * q1[s].b  +  a2 * q2[s].b  +  a3 * q3[s].b  +  a4 * q4[s].b;
			#endif

			#ifdef PIMUNU
				q4[s].pitt = a0 * q[s].pitt  +  a1 * q1[s].pitt  +  a2 * q2[s].pitt  +  a3 * q3[s].pitt  +  a4 * q4[s].pitt;
				q4[s].pitx = a0 * q[s].pitx  +  a1 * q1[s].pitx  +  a2 * q2[s].pitx  +  a3 * q3[s].pitx  +  a4 * q4[s].pitx;
				q4[s].pity = a0 * q[s].pity  +  a1 * q1[s].pity  +  a2 * q2[s].pity  +  a3 * q3[s].pity  +  a4 * q4[s].pity;
			#ifndef BOOST_INVARIANT
				q4[s].pitn = a0 * q[s].pitn  +  a1 * q1[s].pitn  +  a2 * q2[s].pitn  +  a3 * q3[s].pitn  +  a4 * q4[s].pitn;
			#endif
				q4[s].pixx = a0 * q[s].pixx  +  a1 * q1[s].pixx  +  a2 * q2[s].pixx  +  a3 * q3[s].pixx  +  a4 * q4[s].pixx;
				q4[s].pixy = a0 * q[s].pixy  +  a1 * q1[s].pixy  +  a2 * q2[s].pixy  +  a3 * q3[s].pixy  +  a4 * q4[s].pixy;
			#ifndef BOOST_INVARIANT
				q4[s].pixn = a0 * q[s].pixn  +  a1 * q1[s].pixn  +  a2 * q2[s].pixn  +  a3 * q3[s].pixn  +  a4 * q4[s].pixn;
			#endif
				q4[s].piyy = a0 * q[s].piyy  +  a1 * q1[s].piyy  +  a2 * q2[s].piyy  +  a3 * q3[s].piyy  +  a4 * q4[s].piyy;
			#ifndef BOOST_INVARIANT
				q4[s].piyn = a0 * q[s].piyn  +  a1 * q1[s].piyn  +  a2 * q2[s].piyn  +  a3 * q3[s].piyn  +  a4 * q4[s].piyn;
				q4[s].pinn = a0 * q[s].pinn  +  a1 * q1[s].pinn  +  a2 * q2[s].pinn  +  a3 * q3[s].pinn  +  a4 * q4[s].pinn;
			#else
			#ifndef ANISO_HYDRO
				q4[s].pinn = a0 * q[s].pinn  +  a1 * q1[s].pinn  +  a2 * q2[s].pinn  +  a3 * q3[s].pinn  +  a4 * q4[s].pinn;
			#endif
			#endif
			#endif

			#ifdef WTZMU
				q4[s].WtTz = a0 * q[s].WtTz  +  a1 * q1[s].WtTz  +  a2 * q2[s].WtTz  +  a3 * q3[s].WtTz  +  a4 * q4[s].WtTz;
				q4[s].WxTz = a0 * q[s].WxTz  +  a1 * q1[s].WxTz  +  a2 * q2[s].WxTz  +  a3 * q3[s].WxTz  +  a4 * q4[s].WxTz;
				q4[s].WyTz = a0 * q[s].WyTz  +  a1 * q1[s].WyTz  +  a2 * q2[s].WyTz  +  a3 * q3[s].WyTz  +  a4 * q4[s].WyTz;
				q4[s].WnTz = a0 * q[s].WnTz  +  a1 * q1[s].WnTz  +  a2 * q2[s].WnTz  +  a3 * q3[s].WnTz  +  a4 * q4[s].WnTz;
			#endif

			#ifdef PI
				q4[s].Pi = a0 * q[s].Pi  +  a1 * q1[s].Pi  +  a2 * q2[s].Pi  +  a3 * q3[s].Pi  +  a4 * q4[s].Pi;
			#endif
			}
		}
	}
}