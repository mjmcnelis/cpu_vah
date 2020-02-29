
#include "../include/Macros.h"
#include "../include/FreezeoutFinder.h"
#include "../include/Memory.h"
#include "../include/Hydrodynamics.h"
#include "../include/EquationOfState.h"
#include "../../cornelius-c++-1.3/cornelius.cpp"
using namespace std;

inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


freezeout_finder::freezeout_finder(lattice_parameters lattice, hydro_parameters hydro)
{
	freezeout_surface_file.open("output/surface.dat");
  //freezeout_surface_file = fopen("output/surface.dat", "w");

	independent_hydro_variables = 10;

  //precision e_switch = equilibrium_energy_density(hydro.freezeout_temperature_GeV / hbarc, hydro.conformal_eos_prefactor);
  e_switch = (double)equilibrium_energy_density_new(hydro.freezeout_temperature_GeV / hbarc, hydro.conformal_eos_prefactor);

	nx = lattice.lattice_points_x;
	ny = lattice.lattice_points_y;
	nz = lattice.lattice_points_eta;

	dt = (double)lattice.fixed_time_step;
	dx = (double)lattice.lattice_spacing_x;
	dy = (double)lattice.lattice_spacing_y;
	dz = (double)lattice.lattice_spacing_eta;
	tau_coarse_factor = lattice.tau_coarse_factor;

	// initialize cornelius variables for freezeout surface finding (see example_4d() in example_cornelius)
#ifdef BOOST_INVARIANT
	dimension = 3;
	lattice_spacing = new double[dimension];
	lattice_spacing[0] = dt * (double)tau_coarse_factor;
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
 	hypercube = calloc_4d_array(hypercube, 2, 2, 2, 2);
 	cube = calloc_3d_array(cube, 2, 2, 2);
 	hydro_evolution = calloc_5d_array(hydro_evolution, independent_hydro_variables, 2, nx, ny, nz);
}

freezeout_finder::~freezeout_finder()
{

}


void freezeout_finder::set_hydro_evolution(hydro_variables * const __restrict__ q, precision * const __restrict__ e, fluid_velocity * const __restrict__ u)
{
  for(int i = 2; i < nx + 2; i++)
  {
    for(int j = 2; j < ny + 2; j++)
    {
      for(int k = 2; k < nz + 2; k++)
      {
        int s = linear_column_index(i, j, k, nx + 4, ny + 4);

        // set current hydro variables written to first index 1:

        hydro_evolution[0][1][i-2][j-2][k-2] = (double)u[s].ux;
        hydro_evolution[1][1][i-2][j-2][k-2] = (double)u[s].uy;

      #ifndef BOOST_INVARIANT
        hydro_evolution[2][1][i-2][j-2][k-2] = (double)u[s].un;
      #endif

        hydro_evolution[3][1][i-2][j-2][k-2] = (double)e[s];

      #ifdef ANISO_HYDRO                                              // independent anisotropic hydro variables
        hydro_evolution[4][1][i-2][j-2][k-2] = (double)q[s].pl;
        hydro_evolution[5][1][i-2][j-2][k-2] = (double)q[s].pt;

      #ifdef PIMUNU
        hydro_evolution[6][1][i-2][j-2][k-2] = (double)q[s].pixx;
        hydro_evolution[7][1][i-2][j-2][k-2] = (double)q[s].pixy;
      #endif

      #ifdef WTZMU
        hydro_evolution[8][1][i-2][j-2][k-2] = (double)q[s].WTzx;
        hydro_evolution[9][1][i-2][j-2][k-2] = (double)q[s].WTzy;
      #endif

      #else                                                           // independent viscous hydro variables

      #ifdef PIMUNU
        hydro_evolution[4][1][i-2][j-2][k-2] = (double)q[s].pixx;
        hydro_evolution[5][1][i-2][j-2][k-2] = (double)q[s].pixy;
      #ifndef BOOST_INVARIANT
        hydro_evolution[6][1][i-2][j-2][k-2] = (double)q[s].pixn;
      #endif
        hydro_evolution[7][1][i-2][j-2][k-2] = (double)q[s].piyy;
      #ifndef BOOST_INVARIANT
        hydro_evolution[8][1][i-2][j-2][k-2] = (double)q[s].piyn;
      #endif
      #endif

      #ifdef PI
        hydro_evolution[9][1][i-2][j-2][k-2] = (double)q[s].Pi;
      #endif

      #endif
      }
    }
  }
}


void freezeout_finder::swap_and_set_hydro_evolution(hydro_variables * const __restrict__ q, precision * const __restrict__ e, fluid_velocity * const __restrict__ u)
{
  for(int i = 2; i < nx + 2; i++)
  {
    for(int j = 2; j < ny + 2; j++)
    {
      for(int k = 2; k < nz + 2; k++)
      {
        int s = linear_column_index(i, j, k, nx + 4, ny + 4);

        // swap hydro variables from previous freezeout finder call written to zeroth index 0:

        for(int n = 0; n < independent_hydro_variables; n++)
        {
          hydro_evolution[n][0][i-2][j-2][k-2] = hydro_evolution[n][1][i-2][j-2][k-2];
        }

        // set current hydro variables written to first index 1:

        hydro_evolution[0][1][i-2][j-2][k-2] = (double)u[s].ux;
        hydro_evolution[1][1][i-2][j-2][k-2] = (double)u[s].uy;

      #ifndef BOOST_INVARIANT
        hydro_evolution[2][1][i-2][j-2][k-2] = (double)u[s].un;
      #endif

        hydro_evolution[3][1][i-2][j-2][k-2] = (double)e[s];

      #ifdef ANISO_HYDRO														                  // independent anisotropic hydro variables
        hydro_evolution[4][1][i-2][j-2][k-2] = (double)q[s].pl;
        hydro_evolution[5][1][i-2][j-2][k-2] = (double)q[s].pt;

      #ifdef PIMUNU
        hydro_evolution[6][1][i-2][j-2][k-2] = (double)q[s].pixx;
        hydro_evolution[7][1][i-2][j-2][k-2] = (double)q[s].pixy;
      #endif

      #ifdef WTZMU
        hydro_evolution[8][1][i-2][j-2][k-2] = (double)q[s].WTzx;
        hydro_evolution[9][1][i-2][j-2][k-2] = (double)q[s].WTzy;
      #endif

      #else																	                          // independent viscous hydro variables

      #ifdef PIMUNU
        hydro_evolution[4][1][i-2][j-2][k-2] = (double)q[s].pixx;
        hydro_evolution[5][1][i-2][j-2][k-2] = (double)q[s].pixy;
      #ifndef BOOST_INVARIANT
        hydro_evolution[6][1][i-2][j-2][k-2] = (double)q[s].pixn;
      #endif
        hydro_evolution[7][1][i-2][j-2][k-2] = (double)q[s].piyy;
      #ifndef BOOST_INVARIANT
        hydro_evolution[8][1][i-2][j-2][k-2] = (double)q[s].piyn;
      #endif
      #endif

      #ifdef PI
        hydro_evolution[9][1][i-2][j-2][k-2] = (double)q[s].Pi;
      #endif

      #endif
      }
    }
  }
}


void freezeout_finder::construct_energy_density_cube(double ****energy_density, int ix, int iy)
{
  // cube[it][ix][iy]
  // energy_density[it][ix][iy][iz]

  cube[0][0][0] = energy_density[0][ix  ][iy  ][0];
  cube[1][0][0] = energy_density[1][ix  ][iy  ][0];
  cube[0][1][0] = energy_density[0][ix+1][iy  ][0];
  cube[0][0][1] = energy_density[0][ix  ][iy+1][0];
  cube[1][1][0] = energy_density[1][ix+1][iy  ][0];
  cube[1][0][1] = energy_density[1][ix  ][iy+1][0];
  cube[0][1][1] = energy_density[0][ix+1][iy+1][0];
  cube[1][1][1] = energy_density[1][ix+1][iy+1][0];
}


double linear_interpolate_3d(double ****f, int ix, int iy, double x0, double x1, double x2)
{
  // 3d linear interpolation

  return    (1. - x0)  *  (1. - x1)  *  (1. - x2)  *  f[0][ix  ][iy  ][0]
          + (x0)       *  (1. - x1)  *  (1. - x2)  *  f[1][ix  ][iy  ][0]
          + (1. - x0)  *  (x1)       *  (1. - x2)  *  f[0][ix+1][iy  ][0]
          + (1. - x0)  *  (1. - x1)  *  (x2)       *  f[0][ix  ][iy+1][0]
          + (x0)       *  (x1)       *  (1. - x2)  *  f[1][ix+1][iy  ][0]
          + (x0)       *  (1. - x1)  *  (x2)       *  f[1][ix  ][iy+1][0]
          + (1. - x0)  *  (x1)       *  (x2)       *  f[0][ix+1][iy+1][0]
          + (x0)       *  (x1)       *  (x2)       *  f[1][ix+1][iy+1][0];
}



void freezeout_finder::find_2d_freezeout_cells(double t_current, hydro_parameters hydro)
{
  // write freezeout cell's centroid / normal vector and hydro variable to file

  Cornelius cornelius;                                                          // initialize cornelius (can move to class later)
  cornelius.init(dimension, e_switch, lattice_spacing);

  double conformal_prefactor = hydro.conformal_eos_prefactor;

  for(int i = 0; i < nx - 1; i++)
  {
    for(int j = 0; j < ny - 1; j++)
    {
      double cell_t = t_current - lattice_spacing[0];                           // lower left corner of freezeout grid
      double cell_x = dx * ((double)i  -  (double)(nx - 1) / 2.);
      double cell_y = dy * ((double)j  -  (double)(ny - 1) / 2.);
      double cell_z = 0;

      construct_energy_density_cube(hydro_evolution[3], i, j);                  // energy density cube
      cornelius.find_surface_3d(cube);                                          // find centroid and normal vector of each cube

      int freezeout_cells = cornelius.get_Nelements();

      for(int n = 0; n < freezeout_cells; n++)
      {
        // is this used at all?
        //FO_Element fo_cell;                                                 // declare a new fo cell to hold info, later push back to vector

        double t_frac = cornelius.get_centroid_elem(n, 0) / lattice_spacing[0];
        double x_frac = cornelius.get_centroid_elem(n, 1) / lattice_spacing[1];
        double y_frac = cornelius.get_centroid_elem(n, 2) / lattice_spacing[2];

        double t = cornelius.get_centroid_elem(n, 0) + cell_t;                // centroid position of freezeout cell
        double x = cornelius.get_centroid_elem(n, 1) + cell_x;
        double y = cornelius.get_centroid_elem(n, 2) + cell_y;
        double eta = 0;

        // if(!(t_frac >= 0 && t_frac <= 1) || !(x_frac >= 0 && x_frac <= 1) || !(y_frac >= 0 && y_frac <= 1))
        // {
        //   printf("freezeout_finder::find_2d_freezeout_cells error: freezeout cell outside cube\n");
        //   exit(-1);
        // }

        double ds0 = t_current * cornelius.get_normal_elem(n, 0);               // covariant surface normal vector
        double ds1 = t_current * cornelius.get_normal_elem(n, 1);
        double ds2 = t_current * cornelius.get_normal_elem(n, 2);               // don't I want to use tau instead of t_current
        double ds3 = 0;

        // interpolate contravariant flow velocity
        double ux = linear_interpolate_3d(hydro_evolution[0], i, j, t_frac, x_frac, y_frac);
        double uy = linear_interpolate_3d(hydro_evolution[1], i, j, t_frac, x_frac, y_frac);
      #ifndef BOOST_INVARIANT
        double un = linear_interpolate_3d(hydro_evolution[2], i, j, t_frac, x_frac, y_frac);
      #else
        double un = 0;
      #endif

        // interpolate thermodynamic variables
        double e = linear_interpolate_3d(hydro_evolution[3], i, j, t_frac, x_frac, y_frac);

        equation_of_state_new eos(e, conformal_prefactor);
        double T = eos.T;
        double p = eos.equilibrium_pressure();

        // interpolate anisotropic or viscous hydrodynamic variables
      #ifdef ANISO_HYDRO
        double pl = linear_interpolate_3d(hydro_evolution[4], i, j, t_frac, x_frac, y_frac);
        double pt = linear_interpolate_3d(hydro_evolution[5], i, j, t_frac, x_frac, y_frac);

      #ifdef PIMUNU
        double pixx = linear_interpolate_3d(hydro_evolution[6], i, j, t_frac, x_frac, y_frac);
        double pixy = linear_interpolate_3d(hydro_evolution[7], i, j, t_frac, x_frac, y_frac);
      #else
        double pixx = 0;
        double pixy = 0;
      #endif

      #ifdef WTZMU
        double WTzx = linear_interpolate_3d(hydro_evolution[8], i, j, t_frac, x_frac, y_frac);
        double WTzy = linear_interpolate_3d(hydro_evolution[9], i, j, t_frac, x_frac, y_frac);
      #else
        double WTzx = 0;
        double WTzy = 0;
      #endif

        //fprintf(freezeout_surface_file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", t, x, y, eta, ds0, ds1, ds2, ds3, ux);

        // freezeout_surface_file  << t    << " " << x    << " " << y    << " " << eta  << " "
        //                         << ds0  << " " << ds1  << " " << ds2  << " " << ds3  << " "
        //                         << ux   << endl;

        freezeout_surface_file  << t    << " " << x    << " " << y    << " " << eta  << " "
                                << ds0  << " " << ds1  << " " << ds2  << " " << ds3  << " "
                                << ux   << " " << uy   << " " << un   << " "
                                << e    << " " << T    << " " << p    << " "
                                << pl   << " " << pt   << " "
                                << pixx << " " << pixy << " "
                                << WTzx << " " << WTzy << endl;
      #else

      #ifdef PIMUNU
        double pixx = linear_interpolate_3d(hydro_evolution[4], i, j, t_frac, x_frac, y_frac);
        double pixy = linear_interpolate_3d(hydro_evolution[5], i, j, t_frac, x_frac, y_frac);
      #ifndef BOOST_INVARIANT
        double pixn = linear_interpolate_3d(hydro_evolution[6], i, j, t_frac, x_frac, y_frac);
      #else
        double pixn = 0;
      #endif
        double piyy = linear_interpolate_3d(hydro_evolution[7], i, j, t_frac, x_frac, y_frac);
      #ifndef BOOST_INVARIANT
        double piyn = linear_interpolate_3d(hydro_evolution[8], i, j, t_frac, x_frac, y_frac);
      #else
        double piyn = 0;
      #endif

      #else
        double pixx = 0;
        double pixy = 0;
        double pixn = 0;
        double piyy = 0;
        double piyn = 0;
      #endif

      #ifdef PI
        double Pi = linear_interpolate_3d(hydro_evolution[9], i, j, t_frac, x_frac, y_frac);
      #else
        double Pi = 0;
      #endif

        freezeout_surface_file  << t    << " " << x    << " " << y    << " " << eta  << " "
                                << ds0  << " " << ds1  << " " << ds2  << " " << ds3  << " "
                                << ux   << " " << uy   << " " << un   << " "
                                << e    << " " << T    << " " << p    << " "
                                << pixx << " " << pixy << " " << pixn << " " << piyy << " " << piyn << " "
                                << Pi   << endl;
      #endif

      }

    }


  }
}



void freezeout_finder::construct_energy_density_hypercube(double ****energy_density, int ix, int iy, int iz)
{
  // hypercube[it][ix][iy][iz]
  // energy_density[it][ix][iy][iz]

  hypercube[0][0][0][0] = energy_density[0][ix  ][iy  ][iz  ];
  hypercube[1][0][0][0] = energy_density[1][ix  ][iy  ][iz  ];
  hypercube[0][1][0][0] = energy_density[0][ix+1][iy  ][iz  ];
  hypercube[0][0][1][0] = energy_density[0][ix  ][iy+1][iz  ];
  hypercube[0][0][0][1] = energy_density[0][ix  ][iy  ][iz+1];
  hypercube[1][1][0][0] = energy_density[1][ix+1][iy  ][iz  ];
  hypercube[1][0][1][0] = energy_density[1][ix  ][iy+1][iz  ];
  hypercube[1][0][0][1] = energy_density[1][ix  ][iy  ][iz+1];
  hypercube[0][1][1][0] = energy_density[0][ix+1][iy+1][iz  ];
  hypercube[0][1][0][1] = energy_density[0][ix+1][iy  ][iz+1];
  hypercube[0][0][1][1] = energy_density[0][ix  ][iy+1][iz+1];
  hypercube[1][1][1][0] = energy_density[1][ix+1][iy+1][iz  ];
  hypercube[1][1][0][1] = energy_density[1][ix+1][iy  ][iz+1];
  hypercube[1][0][1][1] = energy_density[1][ix  ][iy+1][iz+1];
  hypercube[0][1][1][1] = energy_density[0][ix+1][iy+1][iz+1];
  hypercube[1][1][1][1] = energy_density[1][ix+1][iy+1][iz+1];
}


double linear_interpolate_4d(double ****f, int ix, int iy, int iz, double x0, double x1, double x2, double x3)
{
  // 4d linear interpolation

  return  (1. - x0)  *  (1. - x1)  *  (1. - x2)  *  (1. - x3)  *  f[0][ix  ][iy  ][iz  ]
        + (x0)       *  (1. - x1)  *  (1. - x2)  *  (1. - x3)  *  f[1][ix  ][iy  ][iz  ]
        + (1. - x0)  *  (x1)       *  (1. - x2)  *  (1. - x3)  *  f[0][ix+1][iy  ][iz  ]
        + (1. - x0)  *  (1. - x1)  *  (x2)       *  (1. - x3)  *  f[0][ix  ][iy+1][iz  ]
        + (1. - x0)  *  (1. - x1)  *  (1. - x2)  *  (x3)       *  f[0][ix  ][iy  ][iz+1]
        + (x0)       *  (x1)       *  (1. - x2)  *  (1. - x3)  *  f[1][ix+1][iy  ][iz  ]
        + (x0)       *  (1. - x1)  *  (x2)       *  (1. - x3)  *  f[1][ix  ][iy+1][iz  ]
        + (x0)       *  (1. - x1)  *  (1. - x2)  *  (x3)       *  f[1][ix  ][iy  ][iz+1]
        + (1. - x0)  *  (x1)       *  (x2)       *  (1. - x3)  *  f[0][ix+1][iy+1][iz  ]
        + (1. - x0)  *  (x1)       *  (1. - x2)  *  (x3)       *  f[0][ix+1][iy  ][iz+1]
        + (1. - x0)  *  (1. - x1)  *  (x2)       *  (x3)       *  f[0][ix  ][iy+1][iz+1]
        + (x0)       *  (x1)       *  (x2)       *  (1. - x3)  *  f[1][ix+1][iy+1][iz  ]
        + (x0)       *  (x1)       *  (1. - x2)  *  (x3)       *  f[1][ix+1][iy  ][iz+1]
        + (x0)       *  (1. - x1)  *  (x2)       *  (x3)       *  f[1][ix  ][iy+1][iz+1]
        + (1. - x0)  *  (x1)       *  (x2)       *  (x3)       *  f[0][ix+1][iy+1][iz+1]
        + (x0)       *  (x1)       *  (x2)       *  (x3)       *  f[1][ix+1][iy+1][iz+1];
}


void freezeout_finder::find_3d_freezeout_cells(double t_current, hydro_parameters hydro)
{
  // write freezeout cell's centroid / normal vector and hydro variable to file

  Cornelius cornelius;                                                          // initialize cornelius (can move to class later)
  cornelius.init(dimension, e_switch, lattice_spacing);

  double conformal_prefactor = hydro.conformal_eos_prefactor;

  for(int i = 0; i < nx - 1; i++)
  {
    for(int j = 0; j < ny - 1; j++)
    {
      for(int k = 0; k < nz - 1; k++)
      {
        double cell_t = t_current - lattice_spacing[0];                           // lower left corner of freezeout grid
        double cell_x = dx * ((double)i  -  (double)(nx - 1) / 2.);
        double cell_y = dy * ((double)j  -  (double)(ny - 1) / 2.);
        double cell_z = dz * ((double)k  -  (double)(nz - 1) / 2.);

        construct_energy_density_hypercube(hydro_evolution[3], i, j, k);          // energy density hyper cube
        cornelius.find_surface_4d(hypercube);                                     // find centroid and normal vector of each hypercube

        int freezeout_cells = cornelius.get_Nelements();

        for(int n = 0; n < freezeout_cells; n++)
        {
          // is this used at all?
          //FO_Element fo_cell;                                                 // declare a new fo cell to hold info, later push back to vector

          double t_frac = cornelius.get_centroid_elem(n, 0) / lattice_spacing[0];
          double x_frac = cornelius.get_centroid_elem(n, 1) / lattice_spacing[1];
          double y_frac = cornelius.get_centroid_elem(n, 2) / lattice_spacing[2];
          double z_frac = cornelius.get_centroid_elem(n, 3) / lattice_spacing[3];

          double t   = cornelius.get_centroid_elem(n, 0) + cell_t;              // centroid position of freezeout cell
          double x   = cornelius.get_centroid_elem(n, 1) + cell_x;
          double y   = cornelius.get_centroid_elem(n, 2) + cell_y;
          double eta = cornelius.get_centroid_elem(n, 3) + cell_z;

          // if(!(t_frac >= 0 && t_frac <= 1) || !(x_frac >= 0 && x_frac <= 1) || !(y_frac >= 0 && y_frac <= 1) || !(z_frac >= 0 && z_frac <= 1))
          // {
          //   printf("freezeout_finder::find_2d_freezeout_cells error: freezeout cell outside cube\n");
          //   exit(-1);
          // }

          // I thought needed to scale by t_current for 2+1d case only?
          double ds0 = t_current * cornelius.get_normal_elem(n, 0);               // covariant surface normal vector
          double ds1 = t_current * cornelius.get_normal_elem(n, 1);
          double ds2 = t_current * cornelius.get_normal_elem(n, 2);               // don't I want to use tau instead of t_current
          double ds3 = t_current * cornelius.get_normal_elem(n, 3);

          // interpolate contravariant flow velocity
          double ux = linear_interpolate_4d(hydro_evolution[0], i, j, k, t_frac, x_frac, y_frac, z_frac);
          double uy = linear_interpolate_4d(hydro_evolution[1], i, j, k, t_frac, x_frac, y_frac, z_frac);
        #ifndef BOOST_INVARIANT
          double un = linear_interpolate_4d(hydro_evolution[2], i, j, k, t_frac, x_frac, y_frac, z_frac);
        #else
          double un = 0;
        #endif

          // interpolate thermodynamic variables
          double e = linear_interpolate_4d(hydro_evolution[3], i, j, k, t_frac, x_frac, y_frac, z_frac);

          equation_of_state_new eos(e, conformal_prefactor);
          double T = eos.T;
          double p = eos.equilibrium_pressure();

          // interpolate anisotropic or viscous hydrodynamic variables
        #ifdef ANISO_HYDRO
          double pl = linear_interpolate_4d(hydro_evolution[4], i, j, k, t_frac, x_frac, y_frac, z_frac);
          double pt = linear_interpolate_4d(hydro_evolution[5], i, j, k, t_frac, x_frac, y_frac, z_frac);

        #ifdef PIMUNU
          double pixx = linear_interpolate_4d(hydro_evolution[6], i, j, k, t_frac, x_frac, y_frac, z_frac);
          double pixy = linear_interpolate_4d(hydro_evolution[7], i, j, k, t_frac, x_frac, y_frac, z_frac);
        #else
          double pixx = 0;
          double pixy = 0;
        #endif

        #ifdef WTZMU
          double WTzx = linear_interpolate_4d(hydro_evolution[8], i, j, k, t_frac, x_frac, y_frac, z_frac);
          double WTzy = linear_interpolate_4d(hydro_evolution[9], i, j, k, t_frac, x_frac, y_frac, z_frac);
        #else
          double WTzx = 0;
          double WTzy = 0;
        #endif

          //fprintf(freezeout_surface_file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", t, x, y, eta, ds0, ds1, ds2, ds3, ux);

          // freezeout_surface_file  << t    << " " << x    << " " << y    << " " << eta  << " "
          //                         << ds0  << " " << ds1  << " " << ds2  << " " << ds3  << " "
          //                         << ux   << endl;

          freezeout_surface_file  << t    << " " << x    << " " << y    << " " << eta  << " "
                                  << ds0  << " " << ds1  << " " << ds2  << " " << ds3  << " "
                                  << ux   << " " << uy   << " " << un   << " "
                                  << e    << " " << T    << " " << p    << " "
                                  << pl   << " " << pt   << " "
                                  << pixx << " " << pixy << " "
                                  << WTzx << " " << WTzy << endl;
        #else

        #ifdef PIMUNU
          double pixx = linear_interpolate_4d(hydro_evolution[4], i, j, k, t_frac, x_frac, y_frac, z_frac);
          double pixy = linear_interpolate_4d(hydro_evolution[5], i, j, k, t_frac, x_frac, y_frac, z_frac);
        #ifndef BOOST_INVARIANT
          double pixn = linear_interpolate_4d(hydro_evolution[6], i, j, k, t_frac, x_frac, y_frac, z_frac);
        #else
          double pixn = 0;
        #endif
          double piyy = linear_interpolate_4d(hydro_evolution[7], i, j, k, t_frac, x_frac, y_frac, z_frac);
        #ifndef BOOST_INVARIANT
          double piyn = linear_interpolate_4d(hydro_evolution[8], i, j, k, t_frac, x_frac, y_frac, z_frac);
        #else
          double piyn = 0;
        #endif

        #else
          double pixx = 0;
          double pixy = 0;
          double pixn = 0;
          double piyy = 0;
          double piyn = 0;
        #endif

        #ifdef PI
          double Pi = linear_interpolate_4d(hydro_evolution[9], i, j, k, t_frac, x_frac, y_frac, z_frac);
        #else
          double Pi = 0;
        #endif

          freezeout_surface_file  << t    << " " << x    << " " << y    << " " << eta  << " "
                                  << ds0  << " " << ds1  << " " << ds2  << " " << ds3  << " "
                                  << ux   << " " << uy   << " " << un   << " "
                                  << e    << " " << T    << " " << p    << " "
                                  << pixx << " " << pixy << " " << pixn << " " << piyy << " " << piyn << " "
                                  << Pi   << endl;
        #endif

        }
      }
    }
  }
}


void freezeout_finder::close_file_and_free_memory()
{
	printf("\nFreeing freezeout_finder memory...\n");

  freezeout_surface_file.close();
  //fclose(freezeout_surface_file);

	delete [] lattice_spacing;
	free_5d_array(hydro_evolution, independent_hydro_variables, 2, nx, ny);
	free_4d_array(hypercube, 2, 2, 2);
	free_3d_array(cube, 2, 2);
}















