
#ifndef FREEZEOUTFINDER_H_
#define FREEZEOUTFINDER_H_

#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <vector>
#include "Precision.h"
#include "Macros.h"
#include "Parameters.h"
#include "DynamicalVariables.h"

//#define freezeout_precision float


class freezeout_element
{
 private:		// standard vh format

 public:
   float tau;	// freezeout cell position
   float x;
   float y;
   float eta;

   float dat;	// covariant surface element vector
   float dax;
   float day;
   float dan;

   float ux;	// contravariant fluid velocity
   float uy;
   float un;

   float E;		// energy density [fm^-4]
   float T;		// temperature [fm^-1]
   float P;		// equilibrium pressure [fm^-4]

   				// contravariant shear stress
   float pixx;	// [fm^-4]
   float pixy;	// [fm^-4]
   float pixn;	// [fm^-5]
   float piyy;	// [fm^-4]
   float piyn;	// [fm^-5]

   float Pi;	// bulk pressure [fm^-4]

   freezeout_element(double tau_in, double x_in, double y_in, double eta_in, double dat_in, double dax_in, double day_in, double dan_in, double ux_in, double uy_in, double un_in, double E_in, double T_in, double P_in, double pixx_in, double pixy_in, double pixn_in, double piyy_in, double piyn_in, double Pi_in);
   ~freezeout_element();
};




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

		float *****hydro_evolution;							// for storing hydro information on 2.nx.ny.nz hypergrid

	#ifndef JETSCAPE
		std::ofstream freezeout_surface_file;				// file for output/surface.dat
	#endif

	public:

		freezeout_finder(lattice_parameters lattice, hydro_parameters hydro);
		~freezeout_finder();

	#ifdef JETSCAPE
		std::vector<freezeout_element> freezeout_surface;	// freezeout surface elements stored in vector for JETSCAPE
	#endif


		void set_hydro_evolution(double t_set, hydro_variables * const __restrict__ q, precision * const __restrict__ e, fluid_velocity * const __restrict__ u);
		void swap_and_set_hydro_evolution(hydro_variables * const __restrict__ q, precision * const __restrict__ e, fluid_velocity * const __restrict__ u);

		void construct_energy_density_cube(float ****energy_density, int ix, int iy);
		void construct_energy_density_hypercube(float ****energy_density, int ix, int iy, int iz);

		void find_2d_freezeout_cells(double t_current, hydro_parameters hydro);
		void find_3d_freezeout_cells(double t_current, hydro_parameters hydro);

		void files_and_free_memory(int sample);
};


#endif







