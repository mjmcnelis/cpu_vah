#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"
#include "../include/OpenMP.h"

inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


void set_ghost_cell_boundary_conditions(hydro_variables * const __restrict__ q, precision * const __restrict__ e, fluid_velocity * const __restrict__ u, int s, int sBC)
{
	e[s] = e[sBC];  // s   = index of ghost cell
	u[s] = u[sBC];  // sBC = index of physical cell at boundary
	q[s] = q[sBC];
}


void set_ghost_cells(hydro_variables * const __restrict__ q, precision * const __restrict__ e, fluid_velocity * const __restrict__ u, lattice_parameters lattice)
{
	int nx = lattice.lattice_points_x;
	int ny = lattice.lattice_points_y;
	int nz = lattice.lattice_points_eta;

	// loop over physical (y,z) lattice points
	#pragma omp parallel for collapse(2)
	for(int j = 2; j < ny + 2; j++)
	{
		for(int k = 2; k < nz + 2; k++)
		{
			for(int i = 0; i <= 1; i++)             // set left ghost cells (x)
			{
				int s   = linear_column_index(i, j, k, nx + 4, ny + 4);
				int sBC = linear_column_index(2, j, k, nx + 4, ny + 4);
				set_ghost_cell_boundary_conditions(q, e, u, s, sBC);
			}

			for(int i = nx + 2; i <= nx + 3; i++)   // set right ghost cells (x)
			{
				int s   = linear_column_index(     i, j, k, nx + 4, ny + 4);
				int sBC = linear_column_index(nx + 1, j, k, nx + 4, ny + 4);
				set_ghost_cell_boundary_conditions(q, e, u, s, sBC);
			}
		}
	}

	// loop over physical (x,z) lattice points
	#pragma omp parallel for collapse(2)
	for(int i = 2; i < nx + 2; i++)
	{
		for(int k = 2; k < nz + 2; k++)
		{
			for(int j = 0; j <= 1; j++)             // set left ghost cells (y)
			{
				int s   = linear_column_index(i, j, k, nx + 4, ny + 4);
				int sBC = linear_column_index(i, 2, k, nx + 4, ny + 4);
				set_ghost_cell_boundary_conditions(q, e, u, s, sBC);
			}

			for(int j = ny + 2; j <= ny + 3; j++)   // set right ghost cells (y)
			{
				int s   = linear_column_index(i,      j, k, nx + 4, ny + 4);
				int sBC = linear_column_index(i, ny + 1, k, nx + 4, ny + 4);
				set_ghost_cell_boundary_conditions(q, e, u, s, sBC);
			}
		}
    }
    
    // loop over physical (x,y) lattice points
    #pragma omp parallel for collapse(2)
    for(int i = 2; i < nx + 2; i++)
    {
		for(int j = 2; j < ny + 2; j++)
		{
			for(int k = 0; k <= 1; k++)             // set left ghost cells (z)
			{
				int s   = linear_column_index(i, j, k, nx + 4, ny + 4);
				int sBC = linear_column_index(i, j, 2, nx + 4, ny + 4);
				set_ghost_cell_boundary_conditions(q, e, u, s, sBC);
			}

			for(int k = nz + 2; k <= nz + 3; k++)   // set right ghost cells (z)
			{
				int s   = linear_column_index(i, j,      k, nx + 4, ny + 4);
				int sBC = linear_column_index(i, j, nz + 1, nx + 4, ny + 4);
				set_ghost_cell_boundary_conditions(q, e, u, s, sBC);
			}
		}
	}
}




