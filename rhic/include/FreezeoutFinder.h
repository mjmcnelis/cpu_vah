
#ifndef FREEZEOUTFINDER_H_
#define FREEZEOUTFINDER_H_

#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <vector>
#include "Parameters.h"
#include "DynamicalVariables.h"

class freezeout_finder
{
	private:
		int nx;												// spatial grid parameters
		int ny;
		int nz;
		double dt;
		double dx;
		double dy; 
		double dz;
		double tau_coarse_factor;							// coarse graining factor in tau

		int dimension;										// dimension of freezeout surface 
		double *lattice_spacing;							// spacetime lattice spacing is uniform

		
		int independent_hydro_variables;					// default value is 10 (independent components of Tmunu)

		

	public:
		
		freezeout_finder(lattice_parameters lattice);
		~freezeout_finder();

		double ****hyperCube4D;								// 3+1d hypercube 
		double ***hyperCube3D;								// 2+1d hypercube

		// not sure what this is for yet 
		//std::vector<FO_Element> fo_surf;					// for holding freezeout cell info (for JETSCAPE?)


		// why need energy density evolution if already in hydrodynamic_evolution?


		double ****energy_density_evolution;				// for storing hydro information on Nx.Ny.Nz.tau_coarse_factor hypergrid
		double *****hydrodynamic_evolution;					// ux, uy, un, e, pl, pt, piTxx, piTxy, WTzx, WTzy

		std::ofstream freezeout_surface_file;

		void swap_and_set_energy_density_hydrodynamic_evolution(hydro_variables * const __restrict__ q, precision * const __restrict__ e, fluid_velocity * const __restrict__ u);

		void write_to_file();
		void close_file();
		void free_memory();
};


#endif
