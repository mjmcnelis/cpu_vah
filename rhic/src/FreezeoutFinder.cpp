#include "../include/Macros.h"
#include "../include/FreezeoutFinder.h"
#include "../include/Memory.h"
#include "../include/Hydrodynamics.h"
#include "../include/EquationOfState.h"
#include "../../cornelius-c++-1.3/cornelius.cpp"
#include "../include/OpenMP.h"
using namespace std;

inline int linear_column_index(int i, int j, int k, int nx, int ny)
{
	return i  +  nx * (j  +  ny * k);
}


freezeout_finder::freezeout_finder(lattice_parameters lattice, hydro_parameters hydro)
{
  if(hydro.run_hydro != 2)
  {
    return;
  }

  max_radius = 0;                             // default values
  tau_coord = 0;
  x_coord = 0;
  y_coord = 0;
  eta_coord = 0;

	independent_hydro_variables = 10;

  e_switch = (double)equilibrium_energy_density_new(hydro.freezeout_temperature_GeV / hbarc, hydro.conformal_eos_prefactor);

	nx = lattice.lattice_points_x;
	ny = lattice.lattice_points_y;
	nz = lattice.lattice_points_eta;

	dx = (double)lattice.lattice_spacing_x;
	dy = (double)lattice.lattice_spacing_y;
	dz = (double)lattice.lattice_spacing_eta;

	// initialize cornelius variables for freezeout surface finder (see example_4d() in example_cornelius)
#ifdef BOOST_INVARIANT
	dimension = 3;
	lattice_spacing = new double[dimension];     // lattice_spacing[0] is continuously updated in find_2d_freezeout_cells()
  lattice_spacing[1] = dx;
  lattice_spacing[2] = dy;

	if(!(nx > 1 && ny > 1))
	{
		printf("freezeout_finder::freezeout_finder error: 2d spatial grid needs to have finite size\n");
		exit(-1);
 	}
 	//printf("%lf\t%lf\n", lattice_spacing[1], lattice_spacing[2]);
#else
	dimension = 4;
	lattice_spacing = new double[dimension];     // lattice_spacing[0] is continuously updated in find_3d_freezeout_cells()
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
 	hydro_info = calloc_5d_array(hydro_info, independent_hydro_variables, 2, nx, ny, nz);
}

freezeout_finder::~freezeout_finder()
{

}


void freezeout_finder::load_initial_grid(double t_set, hydro_variables * const __restrict__ q, precision * const __restrict__ e, fluid_velocity * const __restrict__ u)
{
  t_prev = t_set;

  #pragma omp parallel for collapse(3)
  for(int i = 2; i < nx + 2; i++)
  {
    for(int j = 2; j < ny + 2; j++)
    {
      for(int k = 2; k < nz + 2; k++)
      {
        int s = linear_column_index(i, j, k, nx + 4, ny + 4);

        // set current hydro variables written to first index 1:

        hydro_info[0][1][i-2][j-2][k-2] = (float)u[s].ux;
        hydro_info[1][1][i-2][j-2][k-2] = (float)u[s].uy;

      #ifndef BOOST_INVARIANT
        hydro_info[2][1][i-2][j-2][k-2] = (float)u[s].un;
      #endif

        hydro_info[3][1][i-2][j-2][k-2] = (float)e[s];

      #ifdef ANISO_HYDRO                                              // independent anisotropic hydro variables
        hydro_info[4][1][i-2][j-2][k-2] = (float)q[s].pl;
        hydro_info[5][1][i-2][j-2][k-2] = (float)q[s].pt;

      #ifdef PIMUNU
        hydro_info[6][1][i-2][j-2][k-2] = (float)q[s].pixx;
        hydro_info[7][1][i-2][j-2][k-2] = (float)q[s].pixy;
      #endif

      #ifdef WTZMU
        hydro_info[8][1][i-2][j-2][k-2] = (float)q[s].WxTz;
        hydro_info[9][1][i-2][j-2][k-2] = (float)q[s].WyTz;
      #endif

      #else                                                           // independent viscous hydro variables

      #ifdef PIMUNU
        hydro_info[4][1][i-2][j-2][k-2] = (float)q[s].pixx;
        hydro_info[5][1][i-2][j-2][k-2] = (float)q[s].pixy;
      #ifndef BOOST_INVARIANT
        hydro_info[6][1][i-2][j-2][k-2] = (float)q[s].pixn;
      #endif
        hydro_info[7][1][i-2][j-2][k-2] = (float)q[s].piyy;
      #ifndef BOOST_INVARIANT
        hydro_info[8][1][i-2][j-2][k-2] = (float)q[s].piyn;
      #endif
      #endif

      #ifdef PI
        hydro_info[9][1][i-2][j-2][k-2] = (float)q[s].Pi;
      #endif

      #endif
      }
    }
  }
}


void freezeout_finder::load_current_grid(hydro_variables * const __restrict__ q, precision * const __restrict__ e, fluid_velocity * const __restrict__ u)
{
  #pragma omp parallel for collapse(3)
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
          hydro_info[n][0][i-2][j-2][k-2] = hydro_info[n][1][i-2][j-2][k-2];
        }

        // set current hydro variables written to first index 1:

        hydro_info[0][1][i-2][j-2][k-2] = (float)u[s].ux;
        hydro_info[1][1][i-2][j-2][k-2] = (float)u[s].uy;

      #ifndef BOOST_INVARIANT
        hydro_info[2][1][i-2][j-2][k-2] = (float)u[s].un;
      #endif

        hydro_info[3][1][i-2][j-2][k-2] = (float)e[s];

      #ifdef ANISO_HYDRO														                  // independent anisotropic hydro variables
        hydro_info[4][1][i-2][j-2][k-2] = (float)q[s].pl;
        hydro_info[5][1][i-2][j-2][k-2] = (float)q[s].pt;

      #ifdef PIMUNU
        hydro_info[6][1][i-2][j-2][k-2] = (float)q[s].pixx;
        hydro_info[7][1][i-2][j-2][k-2] = (float)q[s].pixy;
      #endif

      #ifdef WTZMU
        hydro_info[8][1][i-2][j-2][k-2] = (float)q[s].WxTz;
        hydro_info[9][1][i-2][j-2][k-2] = (float)q[s].WyTz;
      #endif

      #else																	                          // independent viscous hydro variables

      #ifdef PIMUNU
        hydro_info[4][1][i-2][j-2][k-2] = (float)q[s].pixx;
        hydro_info[5][1][i-2][j-2][k-2] = (float)q[s].pixy;
      #ifndef BOOST_INVARIANT
        hydro_info[6][1][i-2][j-2][k-2] = (float)q[s].pixn;
      #endif
        hydro_info[7][1][i-2][j-2][k-2] = (float)q[s].piyy;
      #ifndef BOOST_INVARIANT
        hydro_info[8][1][i-2][j-2][k-2] = (float)q[s].piyn;
      #endif
      #endif

      #ifdef PI
        hydro_info[9][1][i-2][j-2][k-2] = (float)q[s].Pi;
      #endif

      #endif
      }
    }
  }
}


void freezeout_finder::construct_energy_density_cube(float ****energy_density, int ix, int iy)
{
  // cube[it][ix][iy]
  // energy_density[it][ix][iy][iz]

  cube[0][0][0] = (double)energy_density[0][ix  ][iy  ][0];
  cube[1][0][0] = (double)energy_density[1][ix  ][iy  ][0];
  cube[0][1][0] = (double)energy_density[0][ix+1][iy  ][0];
  cube[0][0][1] = (double)energy_density[0][ix  ][iy+1][0];
  cube[1][1][0] = (double)energy_density[1][ix+1][iy  ][0];
  cube[1][0][1] = (double)energy_density[1][ix  ][iy+1][0];
  cube[0][1][1] = (double)energy_density[0][ix+1][iy+1][0];
  cube[1][1][1] = (double)energy_density[1][ix+1][iy+1][0];
}


void freezeout_finder::construct_energy_density_cube_test(double ***cube_test, float ****energy_density, int ix, int iy)
{
  // cube_test[it][ix][iy]
  // energy_density[it][ix][iy][iz]

  cube_test[0][0][0] = (double)energy_density[0][ix  ][iy  ][0];
  cube_test[1][0][0] = (double)energy_density[1][ix  ][iy  ][0];
  cube_test[0][1][0] = (double)energy_density[0][ix+1][iy  ][0];
  cube_test[0][0][1] = (double)energy_density[0][ix  ][iy+1][0];
  cube_test[1][1][0] = (double)energy_density[1][ix+1][iy  ][0];
  cube_test[1][0][1] = (double)energy_density[1][ix  ][iy+1][0];
  cube_test[0][1][1] = (double)energy_density[0][ix+1][iy+1][0];
  cube_test[1][1][1] = (double)energy_density[1][ix+1][iy+1][0];
}


double linear_interpolate_3d(float ****f, int ix, int iy, double x0, double x1, double x2)
{
  // 3d linear interpolation

  return    (1. - x0)  *  (1. - x1)  *  (1. - x2)  *  (double)f[0][ix  ][iy  ][0]
          + (x0)       *  (1. - x1)  *  (1. - x2)  *  (double)f[1][ix  ][iy  ][0]
          + (1. - x0)  *  (x1)       *  (1. - x2)  *  (double)f[0][ix+1][iy  ][0]
          + (1. - x0)  *  (1. - x1)  *  (x2)       *  (double)f[0][ix  ][iy+1][0]
          + (x0)       *  (x1)       *  (1. - x2)  *  (double)f[1][ix+1][iy  ][0]
          + (x0)       *  (1. - x1)  *  (x2)       *  (double)f[1][ix  ][iy+1][0]
          + (1. - x0)  *  (x1)       *  (x2)       *  (double)f[0][ix+1][iy+1][0]
          + (x0)       *  (x1)       *  (x2)       *  (double)f[1][ix+1][iy+1][0];
}



void freezeout_finder::find_2d_freezeout_cells(double t_current, hydro_parameters hydro)
{
  // write freezeout cell's centroid / normal vector and hydro variables to vectors

  lattice_spacing[0] = t_current - t_prev;                                      // update temporal lattice spacing

  Cornelius cornelius;                                                          // initialize cornelius
  cornelius.init(dimension, e_switch, lattice_spacing);                         // todo: moving inside parallelizable loop is so slow, why?

  double conformal_prefactor = hydro.conformal_eos_prefactor;

#ifdef FREEZEOUT_SIZE
  double r_max_call = 0;
#endif

  for(int i = 0; i < nx - 1; i++)
  {
    for(int j = 0; j < ny - 1; j++)
    {
      double cell_t = t_prev;                                                   // lower left corner of freezeout grid
      double cell_x = dx * ((double)i  -  (double)(nx - 1) / 2.);
      double cell_y = dy * ((double)j  -  (double)(ny - 1) / 2.);
      double cell_z = 0;

      construct_energy_density_cube(hydro_info[3], i, j);                       // energy density cube
      cornelius.find_surface_3d(cube);                                          // find centroid and normal vector of each cube
      int freezeout_cells = cornelius.get_Nelements();

      for(int n = 0; n < freezeout_cells; n++)
      {
        double t_frac = cornelius.get_centroid_elem(n, 0) / lattice_spacing[0];
        double x_frac = cornelius.get_centroid_elem(n, 1) / lattice_spacing[1];
        double y_frac = cornelius.get_centroid_elem(n, 2) / lattice_spacing[2];

        double t = cornelius.get_centroid_elem(n, 0) + cell_t;                // centroid position of freezeout cell
        double x = cornelius.get_centroid_elem(n, 1) + cell_x;
        double y = cornelius.get_centroid_elem(n, 2) + cell_y;
        double eta = 0;

      #ifdef FREEZEOUT_SLICE
        if(fabs(y) <= dy)
        {
          tau_slice_x.push_back(t);           // for freezeout slice plot
          x_slice_x.push_back(x);
        }
      #endif

      #ifdef FREEZEOUT_SIZE
        double r = sqrt(x * x  +  y * y);     // for max fireball radius

        if(r > r_max_call)
        {
          r_max_call = r;
        }
        if(r > max_radius)
        {
          tau_coord = t;
          x_coord = x;
          y_coord = y;
          eta_coord = 0;
          max_radius = r;
        }
      #endif

        // uncomment for debugging
        // if(!(t_frac >= 0 && t_frac <= 1) || !(x_frac >= 0 && x_frac <= 1) || !(y_frac >= 0 && y_frac <= 1))
        // {
        //   printf("freezeout_finder::find_2d_freezeout_cells error: freezeout cell outside cube\n");
        //   exit(-1);
        // }

        // get covariant surface normal vector (need to multiply by jacobian factor tau)
        double ds0 = t_current * cornelius.get_normal_elem(n, 0);
        double ds1 = t_current * cornelius.get_normal_elem(n, 1);
        double ds2 = t_current * cornelius.get_normal_elem(n, 2);
        double ds3 = 0;

        // interpolate contravariant fluid velcoity
        double ux = linear_interpolate_3d(hydro_info[0], i, j, t_frac, x_frac, y_frac);
        double uy = linear_interpolate_3d(hydro_info[1], i, j, t_frac, x_frac, y_frac);
     // double un = linear_interpolate_3d(hydro_info[2], i, j, t_frac, x_frac, y_frac); // leave commented for index tracking
        double un = 0;

        // interpolate thermodynamic variables
        double E = linear_interpolate_3d(hydro_info[3], i, j, t_frac, x_frac, y_frac);
        equation_of_state_new eos(E, conformal_prefactor);
        double T = eos.T;
        double P = eos.equilibrium_pressure();


        // interpolate aniso hydro variables
      #ifdef ANISO_HYDRO
        double pl = linear_interpolate_3d(hydro_info[4], i, j, t_frac, x_frac, y_frac);
        double pt = linear_interpolate_3d(hydro_info[5], i, j, t_frac, x_frac, y_frac);

      #ifdef PIMUNU
        double pixx = linear_interpolate_3d(hydro_info[6], i, j, t_frac, x_frac, y_frac);
        double pixy = linear_interpolate_3d(hydro_info[7], i, j, t_frac, x_frac, y_frac);
      #else
        double pixx = 0;
        double pixy = 0;
      #endif

        double WxTz = 0;
        double WyTz = 0;

        // compute pi^\munu and Pi components (using interpolated values)
        double Dxx = -1.  -  ux * ux;                     // \Delta^\munu
        double Dxy = - ux * uy;
        double Dyy = -1.  -  uy * uy;

        double piyy = (- pixx * (1.  +  uy * uy)  +  2. * pixy * ux * uy) / (1.  +  ux * ux);   // piperp^yy (reconstruction)

        // pi^\munu = (pl - pt)/3 . (\Delta^\munu + 3.z^\mu.z^\nu)  +  piperp^\munu     (Wperp = 0)
        pixx = (pl - pt) * Dxx / 3.  +  pixx;             // pixx on rhs is piperp^xx, etc
        pixy = (pl - pt) * Dxy / 3.  +  pixy;
        double pixn = 0;
        piyy = (pl - pt) * Dyy / 3.  +  piyy;
        double piyn = 0;

        double Pi = (pl + 2.*pt)/3. - P;

      // interpolate viscous hydro variables
      #else

      #ifdef PIMUNU
        double pixx = linear_interpolate_3d(hydro_info[4], i, j, t_frac, x_frac, y_frac);
        double pixy = linear_interpolate_3d(hydro_info[5], i, j, t_frac, x_frac, y_frac);
     // double pixn = linear_interpolate_3d(hydro_info[6], i, j, t_frac, x_frac, y_frac); // leave commented for index tracking
        double pixn = 0;
        double piyy = linear_interpolate_3d(hydro_info[7], i, j, t_frac, x_frac, y_frac);
     // double piyn = linear_interpolate_3d(hydro_info[8], i, j, t_frac, x_frac, y_frac); // leave commented for index tracking
        double piyn = 0;
      #else
        double pixx = 0;
        double pixy = 0;
        double pixn = 0;
        double piyy = 0;
        double piyn = 0;
      #endif

      #ifdef PI
        double Pi = linear_interpolate_3d(hydro_info[9], i, j, t_frac, x_frac, y_frac);
      #else
        double Pi = 0;
      #endif

      #endif

        surface.tau.push_back((float)t);                  // append freezeout surface
        surface.x.push_back(  (float)x);
        surface.y.push_back(  (float)y);
        surface.eta.push_back((float)eta);

        surface.dsigma_tau.push_back((float)ds0);
        surface.dsigma_x.push_back(  (float)ds1);
        surface.dsigma_y.push_back(  (float)ds2);
        surface.dsigma_eta.push_back((float)ds3);

        surface.ux.push_back((float)ux);
        surface.uy.push_back((float)uy);
        surface.un.push_back((float)un);

        surface.E.push_back((float)(E * hbarc));          // undo hbarc = 1 units (for JETSCAPE vectors)
        surface.T.push_back((float)(T * hbarc));
        surface.P.push_back((float)(P * hbarc));

        surface.pixx.push_back((float)(pixx * hbarc));
        surface.pixy.push_back((float)(pixy * hbarc));
        surface.pixn.push_back((float)(pixn * hbarc));
        surface.piyy.push_back((float)(piyy * hbarc));
        surface.piyn.push_back((float)(piyn * hbarc));

        surface.Pi.push_back((float)(Pi * hbarc));

      } // n

    } // j
  } // i

  t_prev = t_current;         // set previous time for the next freezeout finder call
}



void freezeout_finder::construct_energy_density_hypercube(float ****energy_density, int ix, int iy, int iz)
{
  // hypercube[it][ix][iy][iz]
  // energy_density[it][ix][iy][iz]

  hypercube[0][0][0][0] = (double)energy_density[0][ix  ][iy  ][iz  ];
  hypercube[1][0][0][0] = (double)energy_density[1][ix  ][iy  ][iz  ];
  hypercube[0][1][0][0] = (double)energy_density[0][ix+1][iy  ][iz  ];
  hypercube[0][0][1][0] = (double)energy_density[0][ix  ][iy+1][iz  ];
  hypercube[0][0][0][1] = (double)energy_density[0][ix  ][iy  ][iz+1];
  hypercube[1][1][0][0] = (double)energy_density[1][ix+1][iy  ][iz  ];
  hypercube[1][0][1][0] = (double)energy_density[1][ix  ][iy+1][iz  ];
  hypercube[1][0][0][1] = (double)energy_density[1][ix  ][iy  ][iz+1];
  hypercube[0][1][1][0] = (double)energy_density[0][ix+1][iy+1][iz  ];
  hypercube[0][1][0][1] = (double)energy_density[0][ix+1][iy  ][iz+1];
  hypercube[0][0][1][1] = (double)energy_density[0][ix  ][iy+1][iz+1];
  hypercube[1][1][1][0] = (double)energy_density[1][ix+1][iy+1][iz  ];
  hypercube[1][1][0][1] = (double)energy_density[1][ix+1][iy  ][iz+1];
  hypercube[1][0][1][1] = (double)energy_density[1][ix  ][iy+1][iz+1];
  hypercube[0][1][1][1] = (double)energy_density[0][ix+1][iy+1][iz+1];
  hypercube[1][1][1][1] = (double)energy_density[1][ix+1][iy+1][iz+1];
}

void freezeout_finder::construct_energy_density_hypercube_test(double ****hypercube_test, float ****energy_density, int ix, int iy, int iz)
{
  // hypercube_test[it][ix][iy][iz]
  // energy_density[it][ix][iy][iz]

  hypercube_test[0][0][0][0] = (double)energy_density[0][ix  ][iy  ][iz  ];
  hypercube_test[1][0][0][0] = (double)energy_density[1][ix  ][iy  ][iz  ];
  hypercube_test[0][1][0][0] = (double)energy_density[0][ix+1][iy  ][iz  ];
  hypercube_test[0][0][1][0] = (double)energy_density[0][ix  ][iy+1][iz  ];
  hypercube_test[0][0][0][1] = (double)energy_density[0][ix  ][iy  ][iz+1];
  hypercube_test[1][1][0][0] = (double)energy_density[1][ix+1][iy  ][iz  ];
  hypercube_test[1][0][1][0] = (double)energy_density[1][ix  ][iy+1][iz  ];
  hypercube_test[1][0][0][1] = (double)energy_density[1][ix  ][iy  ][iz+1];
  hypercube_test[0][1][1][0] = (double)energy_density[0][ix+1][iy+1][iz  ];
  hypercube_test[0][1][0][1] = (double)energy_density[0][ix+1][iy  ][iz+1];
  hypercube_test[0][0][1][1] = (double)energy_density[0][ix  ][iy+1][iz+1];
  hypercube_test[1][1][1][0] = (double)energy_density[1][ix+1][iy+1][iz  ];
  hypercube_test[1][1][0][1] = (double)energy_density[1][ix+1][iy  ][iz+1];
  hypercube_test[1][0][1][1] = (double)energy_density[1][ix  ][iy+1][iz+1];
  hypercube_test[0][1][1][1] = (double)energy_density[0][ix+1][iy+1][iz+1];
  hypercube_test[1][1][1][1] = (double)energy_density[1][ix+1][iy+1][iz+1];
}


double linear_interpolate_4d(float ****f, int ix, int iy, int iz, double x0, double x1, double x2, double x3)
{
  // 4d linear interpolation

  return  (1. - x0)  *  (1. - x1)  *  (1. - x2)  *  (1. - x3)  *  (double)f[0][ix  ][iy  ][iz  ]
        + (x0)       *  (1. - x1)  *  (1. - x2)  *  (1. - x3)  *  (double)f[1][ix  ][iy  ][iz  ]
        + (1. - x0)  *  (x1)       *  (1. - x2)  *  (1. - x3)  *  (double)f[0][ix+1][iy  ][iz  ]
        + (1. - x0)  *  (1. - x1)  *  (x2)       *  (1. - x3)  *  (double)f[0][ix  ][iy+1][iz  ]
        + (1. - x0)  *  (1. - x1)  *  (1. - x2)  *  (x3)       *  (double)f[0][ix  ][iy  ][iz+1]
        + (x0)       *  (x1)       *  (1. - x2)  *  (1. - x3)  *  (double)f[1][ix+1][iy  ][iz  ]
        + (x0)       *  (1. - x1)  *  (x2)       *  (1. - x3)  *  (double)f[1][ix  ][iy+1][iz  ]
        + (x0)       *  (1. - x1)  *  (1. - x2)  *  (x3)       *  (double)f[1][ix  ][iy  ][iz+1]
        + (1. - x0)  *  (x1)       *  (x2)       *  (1. - x3)  *  (double)f[0][ix+1][iy+1][iz  ]
        + (1. - x0)  *  (x1)       *  (1. - x2)  *  (x3)       *  (double)f[0][ix+1][iy  ][iz+1]
        + (1. - x0)  *  (1. - x1)  *  (x2)       *  (x3)       *  (double)f[0][ix  ][iy+1][iz+1]
        + (x0)       *  (x1)       *  (x2)       *  (1. - x3)  *  (double)f[1][ix+1][iy+1][iz  ]
        + (x0)       *  (x1)       *  (1. - x2)  *  (x3)       *  (double)f[1][ix+1][iy  ][iz+1]
        + (x0)       *  (1. - x1)  *  (x2)       *  (x3)       *  (double)f[1][ix  ][iy+1][iz+1]
        + (1. - x0)  *  (x1)       *  (x2)       *  (x3)       *  (double)f[0][ix+1][iy+1][iz+1]
        + (x0)       *  (x1)       *  (x2)       *  (x3)       *  (double)f[1][ix+1][iy+1][iz+1];
}


void freezeout_finder::find_3d_freezeout_cells(double t_current, hydro_parameters hydro)
{
  // write freezeout cell's centroid / normal vector and hydro variable to file

  lattice_spacing[0] = t_current - t_prev;                                      // update temporal lattice spacing

  Cornelius cornelius;                                                          // initialize cornelius (can move to class later)
  cornelius.init(dimension, e_switch, lattice_spacing);

  double conformal_prefactor = hydro.conformal_eos_prefactor;

  for(int i = 0; i < nx - 1; i++)
  {
    for(int j = 0; j < ny - 1; j++)
    {
      for(int k = 0; k < nz - 1; k++)
      {
        double cell_t = t_prev;                                                 // lower left corner of freezeout grid
        double cell_x = dx * ((double)i  -  (double)(nx - 1) / 2.);
        double cell_y = dy * ((double)j  -  (double)(ny - 1) / 2.);
        double cell_z = dz * ((double)k  -  (double)(nz - 1) / 2.);

        construct_energy_density_hypercube(hydro_info[3], i, j, k);             // energy density hyper cube
        cornelius.find_surface_4d(hypercube);                                   // find centroid and normal vector of each hypercube
        int freezeout_cells = cornelius.get_Nelements();

        for(int n = 0; n < freezeout_cells; n++)
        {
          double t_frac = cornelius.get_centroid_elem(n, 0) / lattice_spacing[0];
          double x_frac = cornelius.get_centroid_elem(n, 1) / lattice_spacing[1];
          double y_frac = cornelius.get_centroid_elem(n, 2) / lattice_spacing[2];
          double z_frac = cornelius.get_centroid_elem(n, 3) / lattice_spacing[3];

          double t   = cornelius.get_centroid_elem(n, 0) + cell_t;              // centroid position of freezeout cell
          double x   = cornelius.get_centroid_elem(n, 1) + cell_x;
          double y   = cornelius.get_centroid_elem(n, 2) + cell_y;
          double eta = cornelius.get_centroid_elem(n, 3) + cell_z;

        #ifdef FREEZEOUT_SLICE
        if(fabs(y) <= dy && fabs(eta) < dz)                                     // for freezeout slice
        {
          tau_slice_x.push_back(t);
          x_slice_x.push_back(x);

        #ifndef BOOST_INVARIANT
          tau_slice_z.push_back(t);
          eta_slice_z.push_back(eta);
        #endif
        }
        #endif

        #ifdef FREEZEOUT_SIZE                                                   // for fireball radius
          double r = sqrt(x * x  +  y * y  +  eta * eta);                       // radius in grid (not cartesian radius)

          if(r > max_radius)
          {
            tau_coord = t;
            x_coord = x;
            y_coord = y;
            eta_coord = eta;
            max_radius = r;
          }
      #endif

          // uncomment for debugging
          // if(!(t_frac >= 0 && t_frac <= 1) || !(x_frac >= 0 && x_frac <= 1) || !(y_frac >= 0 && y_frac <= 1) || !(z_frac >= 0 && z_frac <= 1))
          // {
          //   printf("freezeout_finder::find_2d_freezeout_cells error: freezeout cell outside cube\n");
          //   exit(-1);
          // }

          // get covariant surface normal vector (need to multiply by jacobian factor tau)
          double ds0 = t_current * cornelius.get_normal_elem(n, 0);
          double ds1 = t_current * cornelius.get_normal_elem(n, 1);
          double ds2 = t_current * cornelius.get_normal_elem(n, 2);
          double ds3 = t_current * cornelius.get_normal_elem(n, 3);

          // interpolate contravariant fluid velocity
          double ux = linear_interpolate_4d(hydro_info[0], i, j, k, t_frac, x_frac, y_frac, z_frac);
          double uy = linear_interpolate_4d(hydro_info[1], i, j, k, t_frac, x_frac, y_frac, z_frac);
          double un = linear_interpolate_4d(hydro_info[2], i, j, k, t_frac, x_frac, y_frac, z_frac);

          // interpolate thermodynamic variables
          double E = linear_interpolate_4d(hydro_info[3], i, j, k, t_frac, x_frac, y_frac, z_frac);
          equation_of_state_new eos(E, conformal_prefactor);
          double T = eos.T;
          double P = eos.equilibrium_pressure();


          // interpolate aniso hydro variables
        #ifdef ANISO_HYDRO
          double pl = linear_interpolate_4d(hydro_info[4], i, j, k, t_frac, x_frac, y_frac, z_frac);
          double pt = linear_interpolate_4d(hydro_info[5], i, j, k, t_frac, x_frac, y_frac, z_frac);

        #ifdef PIMUNU
          double pixx = linear_interpolate_4d(hydro_info[6], i, j, k, t_frac, x_frac, y_frac, z_frac);
          double pixy = linear_interpolate_4d(hydro_info[7], i, j, k, t_frac, x_frac, y_frac, z_frac);
        #else
          double pixx = 0;
          double pixy = 0;
        #endif

        #ifdef WTZMU
          double WxTz = linear_interpolate_4d(hydro_info[8], i, j, k, t_frac, x_frac, y_frac, z_frac);
          double WyTz = linear_interpolate_4d(hydro_info[9], i, j, k, t_frac, x_frac, y_frac, z_frac);
        #else
          double WxTz = 0;
          double WyTz = 0;
        #endif

          // compute pi^\munu and Pi components (using interpolated values)

          double ut = sqrt(1.  +  ux * ux  +  uy * uy  +  t * t * un * un);
          double utperp = sqrt(1.  +  ux * ux  +  uy * uy);

          double zt = t * un / utperp;
          double zn = ut / (t * utperp);

          double Dxx = -1.  -  ux * ux;                                       // \Delta^\munu
          double Dxy = - ux * uy;
          double Dxn = - ux * un;
          double Dyy = -1.  -  uy * uy;
          double Dyn = - uy * un;

          double piyy = (- pixx * (1.  +  uy * uy)  +  2. * pixy * ux * uy) / (1.  +  ux * ux);   // piperp (reconstruction)
          double pixn = (pixx * ux  +  pixy * uy) * un / (utperp * utperp);
          double piyn = (pixy * ux  +  piyy * uy) * un / (utperp * utperp);


          // pi^\munu = (pl - pt)/3 .(\Delta^\munu + 3.z^\mu.z^\nu)  +  Wperp^\mu.z^\nu  +  Wperp^\nu.z^\mu  +  piperp^\munu
          pixx = (pl - pt) * Dxx / 3.  +  pixx;                               // pixx on rhs is piperp^xx, etc
          pixy = (pl - pt) * Dxy / 3.  +  pixy;
          pixn = (pl - pt) * Dxn / 3.  +  WxTz * zn  +  pixn;
          piyy = (pl - pt) * Dyy / 3.  +  piyy;
          piyn = (pl - pt) * Dyn / 3.  +  WyTz * zn  +  piyn;

          double Pi = (pl + 2.*pt)/3. - P;


        // interpolate viscous hydro variables
        #else

        #ifdef PIMUNU
          double pixx = linear_interpolate_4d(hydro_info[4], i, j, k, t_frac, x_frac, y_frac, z_frac);
          double pixy = linear_interpolate_4d(hydro_info[5], i, j, k, t_frac, x_frac, y_frac, z_frac);
          double pixn = linear_interpolate_4d(hydro_info[6], i, j, k, t_frac, x_frac, y_frac, z_frac);
          double piyy = linear_interpolate_4d(hydro_info[7], i, j, k, t_frac, x_frac, y_frac, z_frac);
          double piyn = linear_interpolate_4d(hydro_info[8], i, j, k, t_frac, x_frac, y_frac, z_frac);
        #else
          double pixx = 0;
          double pixy = 0;
          double pixn = 0;
          double piyy = 0;
          double piyn = 0;
        #endif

        #ifdef PI
          double Pi = linear_interpolate_4d(hydro_info[9], i, j, k, t_frac, x_frac, y_frac, z_frac);
        #else
          double Pi = 0;
        #endif

        #endif

          surface.tau.push_back((float)t);                  // append freezeout surface
          surface.x.push_back(  (float)x);
          surface.y.push_back(  (float)y);
          surface.eta.push_back((float)eta);

          surface.dsigma_tau.push_back((float)ds0);
          surface.dsigma_x.push_back(  (float)ds1);
          surface.dsigma_y.push_back(  (float)ds2);
          surface.dsigma_eta.push_back((float)ds3);

          surface.ux.push_back((float)ux);
          surface.uy.push_back((float)uy);
          surface.un.push_back((float)un);

          surface.E.push_back((float)(E * hbarc));          // undo hbarc = 1 units
          surface.T.push_back((float)(T * hbarc));
          surface.P.push_back((float)(P * hbarc));

          surface.pixx.push_back((float)(pixx * hbarc));
          surface.pixy.push_back((float)(pixy * hbarc));
          surface.pixn.push_back((float)(pixn * hbarc));
          surface.piyy.push_back((float)(piyy * hbarc));
          surface.piyn.push_back((float)(piyn * hbarc));

          surface.Pi.push_back((float)(Pi * hbarc));

        } // n
      } // k
    } // j
  } // i

  t_prev = t_current;         // set previous time for the next freezeout finder call
}


void freezeout_finder::free_finder_memory(int sample)
{
  // deallocate memory in freezeout finder (except freezeout_surface)
	delete [] lattice_spacing;
	free_5d_array(hydro_info, independent_hydro_variables, 2, nx, ny);
	free_4d_array(hypercube, 2, 2, 2);
	free_3d_array(cube, 2, 2);


// output freezeout slices
#ifdef FREEZEOUT_SLICE
  int cells_x = tau_slice_x.size();

  printf("\nNumber of freezeout cells in tau-x slice = %d\n", cells_x);

  FILE * surface_slice_x;
  char fname[255];
  sprintf(fname, "output/surface_slice_x.dat");
  surface_slice_x = fopen(fname, "w");

  for(int i = 0; i < cells_x; i++)
  {
    // fprintf(surface_slice_x, "%.4e\t%.4e\n", tau_slice_x[i], x_slice_x[i]);
    fprintf(surface_slice_x, "%.3f\t%.4e\n", x_slice_x[i], tau_slice_x[i]);
  }

  fclose(surface_slice_x);
  tau_slice_x.clear();
  x_slice_x.clear();

#ifndef BOOST_INVARIANT
  int cells_z = tau_slice_z.size();

  printf("\nNumber of freezeout cells in tau-eta slice = %d\n", cells_z);

  FILE * surface_slice_z;
  sprintf(fname, "output/surface_slice_z.dat");
  surface_slice_z = fopen(fname, "w");

  for(int i = 0; i < cells_z; i++)
  {
    // fprintf(surface_slice_z, "%.4e\t%.4e\n", tau_slice_z[i], eta_slice_z[i]);
    fprintf(surface_slice_z, "%.3f\t%.4e\n", eta_slice_x[i], tau_slice_z[i]);
  }

  fclose(surface_slice_z);
  tau_slice_z.clear();
  eta_slice_z.clear();
#endif

#endif


// output fireball radius
#ifdef FREEZEOUT_SIZE
  printf("Max radius of freezeout surface = %.3f fm\n", max_radius);

  FILE *fireball_radius;

  if(max_radius > 0)
  {
    if(sample == 0)
    {
      fireball_radius = fopen("output/fireball_radius/fireball_radius.dat", "a");
    }
    else
    {
      char fname[255];
      sprintf(fname, "output/fireball_radius/fireball_radius_%d.dat", sample);
      fireball_radius = fopen(fname, "a");
    }

    fprintf(fireball_radius, "%.6g\t%.6g\t%.6g\t%.6g\t%.6g\n", tau_coord, x_coord, y_coord, eta_coord, max_radius);
    fclose(fireball_radius);
  }
#endif
}















