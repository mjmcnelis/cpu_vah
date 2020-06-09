/*
-------------------------------------------------------
| Code          | CPU VAH
-------------------------------------------------------
| Authors       | Dennis Bazow, Mike McNelis
-------------------------------------------------------
| Date created  | 10/12/2015 by Dennis Bazow
-------------------------------------------------------
| Last edited   | 6/5/2020 by Mike McNelis
-------------------------------------------------------
| Description   | A viscous anisotropic hydrodynamic
|               | simulation of a heavy-ion collision
-------------------------------------------------------
 */

#include "HydroWrapper.h"

int main(int argc, char **argv)
{
  // default: cpu_vah can run as a stand-alone program
  // wrapper: HYDRO class can be instantiated in a larger program (JETSCAPE)

  HYDRO vah;

  vah.start_hydro(argc, argv);						// run hydro and store freezeout surface (Derek: need version w/ no command line arguments?)

  vah.free_freezeout_surface();						// free freezeout surface

  return 0;
}


// example code in JETSCAPE:
/*

HYDRO vah;											// make HYDRO class
vah.read_initial_energy_density(init_e);     		// pass energy density vector from trento
vah.start_hydro();									// run hydro (no cli version??)

IS3D particlization;								// make IS3D class and pass freezeout surface from HYDRO (pinn is extraneous)
particlization.read_fo_surf_from_memory(vah.tau, vah.x, vah.y, vah.eta,
										vah.dsigma_tau, vah.dsigma_x, vah.dsigma_y, vah.dsigma_eta,
										vah.ux, vah.uy, vah.un,
										vah.E, vah.T, vah.P,
										vah.pixx, vah.pixy, vah.pixn, vah.piyy, vah.piyn, vah.pinn,
										vah.Pi);

vah.free_freezeout_surface();						// free surface memory in HYDRO

particlization.run_particlization(0);				// run particle sampler (argument = 0 uses surface from memory)

*/




