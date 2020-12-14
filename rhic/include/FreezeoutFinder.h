
#ifndef FREEZEOUTFINDER_H_
#define FREEZEOUTFINDER_H_

#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include "Precision.h"
#include "Macros.h"
#include "Parameters.h"
#include "DynamicalVariables.h"
#include "FreezeoutSurface.h"
using namespace std;


class freezeout_finder
{
	private:
		double max_radius;									// maximum radius of freezeout surface

		double tau_coord;									// spacetime coordinates of the max radius
		double x_coord;
		double y_coord;
		double eta_coord;

		double t_prev;										// lower time bracket

		int nx;												// spatial grid parameters
		int ny;
		int nz;

		double dx;
		double dy;
		double dz;

		int dimension;										// dimension of freezeout surface
		double *lattice_spacing;							// spacetime lattice spacing is uniform

		double e_switch;									// freezeout energy density
		int independent_hydro_variables;					// default value is 10 (independent components of Tmunu)

		double ***cube;										// 2+1d hypercube
		double ****hypercube;								// 3+1d hypercube

		float *****hydro_info;								// for storing hydro variables on 2.nx.ny.nz hypergrid


	#ifndef JETSCAPE
		std::ofstream freezeout_surface_file;				// file for output/surface.dat
	#endif

	#ifdef FREEZEOUT_SLICE
		std::vector<float> tau_slice_x;
		std::vector<float> x_slice_x;
	#ifndef BOOST_INVARIANT
		std::vector<float> tau_slice_z;
		std::vector<float> eta_slice_z;
	#endif
	#endif

	public:

		freezeout_finder(lattice_parameters lattice, hydro_parameters hydro);
		~freezeout_finder();

		freezeout_surface surface;							// particlization hypersurface

		void load_initial_grid(double t_set, hydro_variables * const __restrict__ q, precision * const __restrict__ e, fluid_velocity * const __restrict__ u);
		void load_current_grid(hydro_variables * const __restrict__ q, precision * const __restrict__ e, fluid_velocity * const __restrict__ u);

		void construct_energy_density_cube(float ****energy_density, int ix, int iy);
		void construct_energy_density_hypercube(float ****energy_density, int ix, int iy, int iz);

		void construct_energy_density_cube_test(double ***cube_test, float ****energy_density, int ix, int iy);
		void construct_energy_density_hypercube_test(double ****hypercube_test, float ****energy_density, int ix, int iy, int iz);

		void find_2d_freezeout_cells(double t_current, hydro_parameters hydro);
		void find_3d_freezeout_cells(double t_current, hydro_parameters hydro);

		void free_finder_memory(int sample);
};


#endif







