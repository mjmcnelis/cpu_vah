#include <iostream>
#include "../include/GhostCells.h"
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"
using namespace std;


inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


void ghost_cell_BC(hydro_variables * const __restrict__ q, precision * const __restrict__ e, fluid_velocity * const __restrict__ u, int s, int sBC)
{
	e[s] = e[sBC];	// set the ghost cell boundary conditions
	u[s] = u[sBC];
	q[s] = q[sBC];	// s = ghost cell index	 |	sBC = physical cell at boundary index
}


void set_ghost_cells(hydro_variables * const __restrict__ q, precision * const __restrict__ e, fluid_velocity * const __restrict__ u, lattice_parameters lattice)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	// loop over the physical (y,z) lattice points
	for(int j = 2; j < ny + 2; j++)
	{
		for(int k = 2; k < nz + 2; k++)
		{
			for(int i = 0; i <= 1; i++)				// set left ghost cells (x)
			{
				int s   = linear_column_index(i, j, k, nx + 4, ny + 4);
				int sBC = linear_column_index(2, j, k, nx + 4, ny + 4);
				ghost_cell_BC(q, e, u, s, sBC);
			}

			for(int i = nx + 2; i <= nx + 3; i++)	// set right ghost cells (x)
			{
				int s   = linear_column_index(     i, j, k, nx + 4, ny + 4);
				int sBC = linear_column_index(nx + 1, j, k, nx + 4, ny + 4);
				ghost_cell_BC(q, e, u, s, sBC);
			}
		}
	}

	// loop over the physical (x,y) and (x,z) lattice points
	for(int i = 2; i < nx + 2; i++)
	{
		for(int k = 2; k < nz + 2; k++)
		{
			for(int j = 0; j <= 1; j++)				// set left ghost cells (y)
			{
				int s   = linear_column_index(i, j, k, nx + 4, ny + 4);
				int sBC = linear_column_index(i, 2, k, nx + 4, ny + 4);
				ghost_cell_BC(q, e, u, s, sBC);
			}

			for(int j = ny + 2; j <= ny + 3; j++)	// set right ghost cells (y)
			{
				int s   = linear_column_index(i,      j, k, nx + 4, ny + 4);
				int sBC = linear_column_index(i, ny + 1, k, nx + 4, ny + 4);
				ghost_cell_BC(q, e, u, s, sBC);
			}
		}

		for(int j = 2; j < ny + 2; j++)
		{
			for(int k = 0; k <= 1; k++)				// set left ghost cells (z)
			{
				int s   = linear_column_index(i, j, k, nx + 4, ny + 4);
				int sBC = linear_column_index(i, j, 2, nx + 4, ny + 4);
				ghost_cell_BC(q, e, u, s, sBC);
			}

			for(int k = nz + 2; k <= nz + 3; k++)	// set right ghost cells (z)
			{
				int s   = linear_column_index(i, j,      k, nx + 4, ny + 4);
				int sBC = linear_column_index(i, j, nz + 1, nx + 4, ny + 4);
				ghost_cell_BC(q, e, u, s, sBC);
			}
		}
	}
}




