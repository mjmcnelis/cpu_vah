
#include "../include/GhostCells.h"
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"


inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


void ghost_cell_BC(CONSERVED_VARIABLES * const __restrict__ q, precision * const __restrict__ e, FLUID_VELOCITY * const __restrict__ u, int s, int sBC)
{
	// s = ghost cell index	 |	sBC = physical cell at boundary index

	// set the ghost cells
	e[s] = e[sBC];			// this is no longer needed

	u->ut[s] = u->ut[sBC];
	u->ux[s] = u->ux[sBC];
	u->uy[s] = u->uy[sBC];
	u->un[s] = u->un[sBC];

	q->ttt[s] = q->ttt[sBC];
	q->ttx[s] = q->ttx[sBC];
	q->tty[s] = q->tty[sBC];
	q->ttn[s] = q->ttn[sBC];

	q->pl[s] = q->pl[sBC];
	q->pt[s] = q->pt[sBC];

#ifdef PIMUNU
	q->pitt[s] = q->pitt[sBC];
	q->pitx[s] = q->pitx[sBC];
	q->pity[s] = q->pity[sBC];
	q->pitn[s] = q->pitn[sBC];
	q->pixx[s] = q->pixx[sBC];
	q->pixy[s] = q->pixy[sBC];
	q->pixn[s] = q->pixn[sBC];
	q->piyy[s] = q->piyy[sBC];
	q->piyn[s] = q->piyn[sBC];
	q->pinn[s] = q->pinn[sBC];
#endif
#ifdef WTZMU
	q->WtTz[s] = q->WtTz[sBC];
	q->WxTz[s] = q->WxTz[sBC];
	q->WyTz[s] = q->WyTz[sBC];
	q->WnTz[s] = q->WnTz[sBC];
#endif
}



void set_ghost_cells(CONSERVED_VARIABLES * const __restrict__ q, precision * const __restrict__ e, FLUID_VELOCITY * const __restrict__ u, int nx, int ny, int nz)
{
	int s;		// ghost cell index
	int sBC;	// physical boundary index

	// loop over the physical (y,z) lattice points
	for(int j = 2; j < ny + 2; j++)
	{
		for(int k = 2; k < nz + 2; k++)
		{
			for(int i = 0; i <= 1; i++)				// set left ghost cells (x)
			{
				s = linear_column_index(i, j, k, nx + 4, ny + 4);
				sBC = linear_column_index(2, j, k, nx + 4, ny + 4);
				ghost_cell_BC(q, e, u, s, sBC);
			}

			for(int i = nx + 2; i <= nx + 3; i++)	// set right ghost cells (x)
			{
				s = linear_column_index(i, j, k, nx + 4, ny + 4);
				sBC = linear_column_index(nx + 1, j, k, nx + 4, ny + 4);
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
				s = linear_column_index(i, j, k, nx + 4, ny + 4);
				sBC = linear_column_index(i, 2, k, nx + 4, ny + 4);
				ghost_cell_BC(q, e, u, s, sBC);
			}

			for(int j = ny + 2; j <= ny + 3; j++)	// set right ghost cells (y)
			{
				s = linear_column_index(i, j, k, nx + 4, ny + 4);
				sBC = linear_column_index(i, ny + 1, k, nx + 4, ny + 4);
				ghost_cell_BC(q, e, u, s, sBC);
			}
		}

		for(int j = 2; j < ny + 2; j++)
		{
			for(int k = 0; k <= 1; k++)				// set left ghost cells (z)
			{
				s = linear_column_index(i, j, k, nx + 4, ny + 4);
				sBC = linear_column_index(i, j, 2, nx + 4, ny + 4);
				ghost_cell_BC(q, e, u, s, sBC);
			}

			for(int k = nz + 2; k <= nz + 3; k++)	// set right ghost cells (z)
			{
				s = linear_column_index(i, j, k, nx + 4, ny + 4);
				sBC = linear_column_index(i, j, nz + 1, nx + 4, ny + 4);
				ghost_cell_BC(q, e, u, s, sBC);
			}
		}
	}
}




