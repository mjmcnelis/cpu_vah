
#include "../include/Macros.h"
#include "../include/FreezeoutFinder.h"
#include "../include/Memory.h"
using namespace std;

inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


freezeout_finder::freezeout_finder(lattice_parameters lattice)
{
	freezeout_surface_file.open("output/surface.dat");

	independent_hydro_variables = 10;

	nx = lattice.lattice_points_x;
	ny = lattice.lattice_points_y;
	nz = lattice.lattice_points_eta;

	dt = lattice.fixed_time_step;
	dx = lattice.lattice_spacing_x;
	dy = lattice.lattice_spacing_y;
	dz = lattice.lattice_spacing_eta;
	tau_coarse_factor = lattice.tau_coarse_factor;

	// initialize cornelius variables for freezeout surface finding (see example_4d() in example_cornelius)
#ifdef BOOST_INVARIANT
	dimension = 3;
	lattice_spacing = new double[dimension];
	lattice_spacing[0] = dt * tau_coarse_factor;  
    lattice_spacing[1] = dx;                              
    lattice_spacing[2] = dy;

	if(!(nx > 1 && ny > 1))
	{
		printf("freezeout_finder::freezeout_finder error: 2d spatial grid needs to have finite size\n");
		exit(-1);
 	}
 	//printf("%lf\t%lf\t%lf\n", lattice_spacing[0], lattice_spacing[1], lattice_spacing[2]);
#else
	dimension = 4;
	lattice_spacing = new double[dimension];
	lattice_spacing[0] = dt * tau_coarse_factor;  
    lattice_spacing[1] = dx;                              
    lattice_spacing[2] = dy;
    lattice_spacing[3] = dz;

    if(!(nx > 1 && ny > 1 && nz > 1))
	{
		printf("freezeout_finder::freezeout_finder error: 3d spatial grid needs to have finite size\n");
		exit(-1);
 	}
#endif


 	// allocate memory
 	hyperCube4D = calloc_4d_array(hyperCube4D, 2, 2, 2, 2);
 	hyperCube3D = calloc_3d_array(hyperCube3D, 2, 2, 2);

 	energy_density_evolution = calloc_4d_array(energy_density_evolution, 2, nx, ny, nz);
 	hydrodynamic_evolution = calloc_5d_array(hydrodynamic_evolution, independent_hydro_variables, 2, nx, ny, nz);
}

freezeout_finder::~freezeout_finder()
{
	
}



void freezeout_finder::swap_and_set_energy_density_hydrodynamic_evolution(hydro_variables * const __restrict__ q, precision * const __restrict__ e, fluid_velocity * const __restrict__ u)
{
  for(int i = 2; i < nx + 2; i++)
  {
    for(int j = 2; j < ny + 2; j++)
    {
      for(int k = 2; k < nz + 2; k++)
      {
        // hydro variables from previous freezeout finder call written to zeroth index 0:

        energy_density_evolution[0][i-2][j-2][k-2] = energy_density_evolution[1][i-2][j-2][k-2];

        for(int n = 0; n < independent_hydro_variables; n++)
        {
          hydrodynamic_evolution[n][0][i-2][j-2][k-2] = hydrodynamic_evolution[n][1][i-2][j-2][k-2];
        }


        // current hydro variables written to first index 1:
        int s = linear_column_index(i, j, k, nx + 4, ny + 4);
        
        energy_density_evolution[1][i-2][j-2][k-2] = e[s];					// technically don't need this (just a matter of convenience)
       

        hydrodynamic_evolution[0][1][i-2][j-2][k-2] = u[s].ux;				// can't forget this is the order I'm doing it in 
        hydrodynamic_evolution[1][1][i-2][j-2][k-2] = u[s].uy;

    #ifndef BOOST_INVARIANT
        hydrodynamic_evolution[2][1][i-2][j-2][k-2] = u[s].un;
    #endif

        hydrodynamic_evolution[3][1][i-2][j-2][k-2] = e[s];

    #ifdef ANISO_HYDRO														// independent anisotropic hydro variables 
        hydrodynamic_evolution[4][1][i-2][j-2][k-2] = q[s].pl;
        hydrodynamic_evolution[5][1][i-2][j-2][k-2] = q[s].pt;

	#ifdef PIMUNU
        hydrodynamic_evolution[6][1][i-2][j-2][k-2] = q[s].pixx;
        hydrodynamic_evolution[7][1][i-2][j-2][k-2] = q[s].pixy;
	#endif

	#ifdef WTZMU
        hydrodynamic_evolution[8][1][i-2][j-2][k-2] = q[s].WTzx;
        hydrodynamic_evolution[9][1][i-2][j-2][k-2] = q[s].WTzy;
	#endif

    #else																	// independent viscous hydro variables

	#ifdef PIMUNU
        hydrodynamic_evolution[4][1][i-2][j-2][k-2] = q[s].pixx;
        hydrodynamic_evolution[5][1][i-2][j-2][k-2] = q[s].pixy;
        hydrodynamic_evolution[6][1][i-2][j-2][k-2] = q[s].pixn;
        hydrodynamic_evolution[7][1][i-2][j-2][k-2] = q[s].piyy;
        hydrodynamic_evolution[8][1][i-2][j-2][k-2] = q[s].piyn;
	#endif

	#ifdef PI
        hydrodynamic_evolution[9][1][i-2][j-2][k-2] = q[s].Pi;
	#endif

    #endif
      }
    } 
  } 
}



void freezeout_finder::call_boost_invariant_freezeout_finder(int dim, int nx, int ny, int nz, int n, double t0, double dt, double t, double dx, double dy, double dz, double *lattice_spacing, double freezeoutEnergyDensity,
  double ****hyperCube4D, double ***hyperCube3D, double ****energy_density_evoution, double *****hydrodynamic_evoution,
  std::ofstream& freezeoutSurfaceFile, std::vector<FO_Element>& fo_surf, EOS eqnOfState)
  {
    //besides writing centroid and normal to file, write all the hydro variables
    //#pragma omp parallel for collapse(2)
    for (int ix = 0; ix < nx-1; ix++)
    {
      for (int iy = 0; iy < ny-1; iy++)
      {
        Cornelius cor;
        cor.init(dim, freezeoutEnergyDensity, lattice_spacing);
        //write the values of energy density to all corners of the hyperCube
        writeEnergyDensityToHypercube3D(hyperCube3D, energy_density_evoution, 0, ix, iy);
        //use cornelius to find the centroid and normal vector of each hyperCube
        cor.find_surface_3d(hyperCube3D);
        //write centroid and normal of each surface element to file
        for (int i = 0; i < cor.get_Nelements(); i++)
        {
          //declare a new fo cell to hold info, later push back to vector
          FO_Element fo_cell;

          //first write the position of the centroid of surface element
          double cell_tau = t0 + ((double)n) * dt; //check if this is the correct time!
          double cell_x = (double)ix * dx  - (((double)(nx-1)) / 2.0 * dx);
          double cell_y = (double)iy * dy  - (((double)(ny-1)) / 2.0 * dy);
          double cell_z = 0.0;

          double tau_frac = cor.get_centroid_elem(i,0) / lattice_spacing[0];
          double x_frac = cor.get_centroid_elem(i,1) / lattice_spacing[1];
          double y_frac = cor.get_centroid_elem(i,2) / lattice_spacing[2];

          //cell position
          double tau = cor.get_centroid_elem(i,0) + cell_tau;
          double x = cor.get_centroid_elem(i,1) + cell_x;
          double y = cor.get_centroid_elem(i,2) + cell_y;
          double eta = 0.0;

          //covariant surface normal vector
          double ds0 = t * cor.get_normal_elem(i,0);
          double ds1 = t * cor.get_normal_elem(i,1);
          double ds2 = t * cor.get_normal_elem(i,2);
          double ds3 = 0.0;

          //contravariant flow velocity
          double ux = interpolateVariable3D(hydrodynamic_evoution, 0, 0, ix, iy, tau_frac, x_frac, y_frac);
          double uy = interpolateVariable3D(hydrodynamic_evoution, 1, 0, ix, iy, tau_frac, x_frac, y_frac);
          double un = interpolateVariable3D(hydrodynamic_evoution, 2, 0, ix, iy, tau_frac, x_frac, y_frac);

          //energy density, Temperature, Pressure
          double eps = interpolateVariable3D(hydrodynamic_evoution, 3, 0, ix, iy, tau_frac, x_frac, y_frac);
          double T = eqnOfState.effectiveTemperature(eps);
          double P = eqnOfState.equilibriumPressure(eps);

          //contravariant components of shear stress
          double pixx = interpolateVariable3D(hydrodynamic_evoution, 4, 0, ix, iy, tau_frac, x_frac, y_frac);
          double pixy = interpolateVariable3D(hydrodynamic_evoution, 5, 0, ix, iy, tau_frac, x_frac, y_frac);
          double pixn = interpolateVariable3D(hydrodynamic_evoution, 6, 0, ix, iy, tau_frac, x_frac, y_frac);
          double piyy = interpolateVariable3D(hydrodynamic_evoution, 7, 0, ix, iy, tau_frac, x_frac, y_frac);
          double piyn = interpolateVariable3D(hydrodynamic_evoution, 8, 0, ix, iy, tau_frac, x_frac, y_frac);
          //bulk pressure
          double Pi = interpolateVariable3D(hydrodynamic_evoution, 9, 0, ix, iy, tau_frac, x_frac, y_frac);

          #pragma omp critical
          {
            freezeoutSurfaceFile << tau  << " " <<  x   << " " <<  y   << " " << eta  << " "
                                 << ds0  << " " << ds1  << " " << ds2  << " " << ds3  << " "
                                 <<  ux  << " " <<  uy  << " " <<  un  << " "
                                 << eps  << " " <<  T   << " " <<  P   << " "
                                 << pixx << " " << pixy << " " << pixn << " " << piyy << " " << piyn << " "
                                 <<  Pi  << endl;
          }
          //add the fo cell to fo surface
          //#pragma omp critical
          //if (SAVE_FO_SURF_VECTOR) fo_surf.push_back(fo_cell);

        } //for (int i = 0; i < cor.get_Nelements(); i++)
      } // for (int iy = 0; iy < ny-1; iy++)
    } // for (int ix = 0; ix < nx-1; ix++)
  }


void freezeout_finder::write_to_file()
{
	freezeout_surface_file << 1.234 << " " << 5.6 << endl;
}


void freezeout_finder::close_file()
{
	freezeout_surface_file.close();
}


void freezeout_finder::free_memory()
{
	printf("\nFreeing freezeout_finder memory...\n");

	delete [] lattice_spacing; 
	free_4d_array(energy_density_evolution, 2, nx, ny);
	free_5d_array(hydrodynamic_evolution, independent_hydro_variables, 2, nx, ny);
	free_4d_array(hyperCube4D, 2, 2, 2);
	free_3d_array(hyperCube3D, 2, 2);
}















