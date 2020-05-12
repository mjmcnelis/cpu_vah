
#ifndef FREEZEOUTFINDER_H_
#define FREEZEOUTFINDER_H_

#include <stdlib.h>
#include <stdio.h>
#include <fstream>
//#include <vector>
#include "Precision.h"
#include "Parameters.h"
#include "DynamicalVariables.h"

class freezeout_finder
{
	private:
		double t_prev;										// lower time bracket

		int nx;												// spatial grid parameters
		int ny;
		int nz;

		double dx;
		double dy;
		double dz;
		double tau_coarse_factor;							// coarse graining factor in tau

		int dimension;										// dimension of freezeout surface
		double *lattice_spacing;							// spacetime lattice spacing is uniform

		double e_switch;									// freezeout energy density
		int independent_hydro_variables;					// default value is 10 (independent components of Tmunu)

		double ***cube;										// 2+1d hypercube
		double ****hypercube;								// 3+1d hypercube

		double *****hydro_evolution;						// for storing hydro information on 2.nx.ny.nz hypergrid

		std::ofstream freezeout_surface_file;
		// FILE * freezeout_surface_file;						// freezeout surface file

	public:

		freezeout_finder(lattice_parameters lattice, hydro_parameters hydro);
		~freezeout_finder();

		// not sure what this is for yet
		//std::vector<FO_Element> fo_surf;					// for holding freezeout cell info (for JETSCAPE?)
		void set_hydro_evolution(double t_set, hydro_variables * const __restrict__ q, precision * const __restrict__ e, fluid_velocity * const __restrict__ u);
		void swap_and_set_hydro_evolution(hydro_variables * const __restrict__ q, precision * const __restrict__ e, fluid_velocity * const __restrict__ u);

		void construct_energy_density_cube(double ****energy_density, int ix, int iy);
		void construct_energy_density_hypercube(double ****energy_density, int ix, int iy, int iz);

		void find_2d_freezeout_cells(double t_current, hydro_parameters hydro);
		void find_3d_freezeout_cells(double t_current, hydro_parameters hydro);

		void close_file_and_free_memory();
};


#endif







