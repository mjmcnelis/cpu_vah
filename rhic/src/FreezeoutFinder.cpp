
#include "../include/Macros.h"
#include "../include/FreezeoutFinder.h"
#include "../include/Memory.h"
#include "../include/Hydrodynamics.h"
using namespace std;

inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


freezeout_finder::freezeout_finder(lattice_parameters lattice, hydro_parameters hydro)
{
	freezeout_surface_file.open("output/surface.dat");

	independent_hydro_variables = 10;

  precision T_switch = hydro.freezeout_temperature_GeV;
  //precision e_switch = equilibrium_energy_density(T_switch / hbarc, hydro.conformal_eos_prefactor);
  e_switch = equilibrium_energy_density_new(T_switch / hbarc, hydro.conformal_eos_prefactor);     

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
 	hydro_evolution = calloc_5d_array(hydro_evolution, independent_hydro_variables, 2, nx, ny, nz);
}

freezeout_finder::~freezeout_finder()
{
	
}



void freezeout_finder::swap_and_set_energy_density_hydro_evolution(hydro_variables * const __restrict__ q, precision * const __restrict__ e, fluid_velocity * const __restrict__ u)
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
          hydro_evolution[n][0][i-2][j-2][k-2] = hydro_evolution[n][1][i-2][j-2][k-2];
        }

        int s = linear_column_index(i, j, k, nx + 4, ny + 4); 
        

        // current hydro variables written to first index 1:

        energy_density_evolution[1][i-2][j-2][k-2] = e[s];		// technically don't need this (just a matter of convenience)
       

        hydro_evolution[0][1][i-2][j-2][k-2] = u[s].ux;				// can't forget this is the order I'm doing it in 
        hydro_evolution[1][1][i-2][j-2][k-2] = u[s].uy;

      #ifndef BOOST_INVARIANT
        hydro_evolution[2][1][i-2][j-2][k-2] = u[s].un;
      #endif

        hydro_evolution[3][1][i-2][j-2][k-2] = e[s];

      #ifdef ANISO_HYDRO														          // independent anisotropic hydro variables 
        hydro_evolution[4][1][i-2][j-2][k-2] = q[s].pl;
        hydro_evolution[5][1][i-2][j-2][k-2] = q[s].pt;

      #ifdef PIMUNU
        hydro_evolution[6][1][i-2][j-2][k-2] = q[s].pixx;
        hydro_evolution[7][1][i-2][j-2][k-2] = q[s].pixy;
      #endif

      #ifdef WTZMU
        hydro_evolution[8][1][i-2][j-2][k-2] = q[s].WTzx;
        hydro_evolution[9][1][i-2][j-2][k-2] = q[s].WTzy;
      #endif

      #else																	                 // independent viscous hydro variables

      #ifdef PIMUNU
        hydro_evolution[4][1][i-2][j-2][k-2] = q[s].pixx;
        hydro_evolution[5][1][i-2][j-2][k-2] = q[s].pixy;
        hydro_evolution[6][1][i-2][j-2][k-2] = q[s].pixn;
        hydro_evolution[7][1][i-2][j-2][k-2] = q[s].piyy;
        hydro_evolution[8][1][i-2][j-2][k-2] = q[s].piyn;
      #endif

      #ifdef PI
        hydro_evolution[9][1][i-2][j-2][k-2] = q[s].Pi;
      #endif

      #endif
      }
    } 
  } 
}



void freezeout_finder::call_2d_freezeout_finder(double t_current)
{
  // write centroid and normal to file, write all the hydro variables to file

  for(int ix = 0; ix < nx-1; ix++)                
  {
    for(int iy = 0; iy < ny-1; iy++)
    {
      double cell_t = t_current;                                                // position of current fluid cell        
      double cell_x = (double)ix * dx  - (((double)(nx-1)) / 2.0 * dx);         
      double cell_y = (double)iy * dy  - (((double)(ny-1)) / 2.0 * dy);
      double cell_z = 0;


      Cornelius cornelius;                                                      // initialize cornelius
      cornelius.init(dimension, e_switch, lattice_spacing);
      
      writeEnergyDensityToHypercube3D(hyperCube3D, energy_density_evoution, 0, ix, iy);   // construct hyper cube

      
      cornelius.find_surface_3d(hyperCube3D);                                 // use cornelius to find the centroid and normal vector of each hyper cube

      for(int s = 0; s < cornelius.get_Nelements(); s++)
      {
        // is this used at all?
        //FO_Element fo_cell;                                                 // declare a new fo cell to hold info, later push back to vector

      
        double t_frac = cornelius.get_centroid_elem(s, 0) / lattice_spacing[0];
        double x_frac = cornelius.get_centroid_elem(s, 1) / lattice_spacing[1];
        double y_frac = cornelius.get_centroid_elem(s, 2) / lattice_spacing[2];

        
        double tau = cornelius.get_centroid_elem(s, 0) + cell_t;                // centroid position of freezeout cell
        double x = cornelius.get_centroid_elem(s, 1) + cell_x;
        double y = cornelius.get_centroid_elem(s, 2) + cell_y;
        double eta = 0;

        
        double ds0 = t_current * cornelius.get_normal_elem(s, 0);               // covariant surface normal vector
        double ds1 = t_current * cornelius.get_normal_elem(s, 1);
        double ds2 = t_current * cornelius.get_normal_elem(s, 2);               // don't I want to use tau instead of t_current 
        double ds3 = 0;


        //contravariant flow velocity
        double ux = interpolateVariable3D(hydro_evolution, 0, 0, ix, iy, t_frac, x_frac, y_frac);
        double uy = interpolateVariable3D(hydro_evolution, 1, 0, ix, iy, t_frac, x_frac, y_frac);
        double un = interpolateVariable3D(hydro_evolution, 2, 0, ix, iy, t_frac, x_frac, y_frac);


        // energy density, Temperature, Pressure
        double e = interpolateVariable3D(hydro_evolution, 3, 0, ix, iy, t_frac, x_frac, y_frac);

        equation_of_state_new eos(e);
        double T = eos.effective_temperature();
        double p = eos.equilibrium_pressure();


        //contravariant components of shear stress
        double pixx = interpolateVariable3D(hydro_evolution, 4, 0, ix, iy, t_frac, x_frac, y_frac);
        double pixy = interpolateVariable3D(hydro_evolution, 5, 0, ix, iy, t_frac, x_frac, y_frac);
        double pixn = interpolateVariable3D(hydro_evolution, 6, 0, ix, iy, t_frac, x_frac, y_frac);
        double piyy = interpolateVariable3D(hydro_evolution, 7, 0, ix, iy, t_frac, x_frac, y_frac);
        double piyn = interpolateVariable3D(hydro_evolution, 8, 0, ix, iy, t_frac, x_frac, y_frac);
        //bulk pressure
        double Pi = interpolateVariable3D(hydro_evolution, 9, 0, ix, iy, t_frac, x_frac, y_frac);

        
        freezeout_surface_file  << tau  << " " <<  x   << " " <<  y   << " " << eta  << " "
                                << ds0  << " " << ds1  << " " << ds2  << " " << ds3  << " "
                                <<  ux  << " " <<  uy  << " " <<  un  << " "
                                << eps  << " " <<  T   << " " <<  p   << " "
                                << pixx << " " << pixy << " " << pixn << " " << piyy << " " << piyn << " "
                                <<  Pi  << endl;

      } //for (int i = 0; i < cor.get_Nelements(); i++)

    } // for (int iy = 0; iy < ny-1; iy++)

  } // for (int ix = 0; ix < nx-1; ix++)
}


void freezeout_finder::close_file_and_free_memory()
{
	printf("\nFreeing freezeout_finder memory...\n");

  freezeout_surface_file.close();

	delete [] lattice_spacing; 
	free_4d_array(energy_density_evolution, 2, nx, ny);
	free_5d_array(hydro_evolution, independent_hydro_variables, 2, nx, ny);
	free_4d_array(hyperCube4D, 2, 2, 2);
	free_3d_array(hyperCube3D, 2, 2);
}















