
#ifndef FREEZEOUTFINDER_H_
#define FREEZEOUTFINDER_H_

#include "Precision.h"
#include "Parameters.h"

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
	public:
		
		freezeout_finder(lattice_parameters lattice_in);
		~freezeout_finder();
};


#endif
