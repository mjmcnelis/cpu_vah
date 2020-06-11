/*
-------------------------------------------------------
| Code          | CPU VAH
-------------------------------------------------------
| Authors       | Dennis Bazow, Mike McNelis
-------------------------------------------------------
| Date created  | 10/12/2015 by Dennis Bazow
-------------------------------------------------------
| Last edited   | 6/9/2020 by Mike McNelis
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

  if(argc > 1)                          // run hydro and store freezeout surface
  {
    vah.start_hydro(argc, argv);
  }
  else
  {
    vah.start_hydro_no_arguments();
  }

  vah.free_freezeout_surface();         // free freezeout surface

  return 0;
}


// example code in JETSCAPE:
/*

HYDRO vah;                                              // make HYDRO class
vah.read_trento_energy_density_profile(trento_energy);  // pass energy density vector from trento
vah.start_hydro_no_arguments();                         // run hydro (Derek: do I need version with no command line arguments?)

IS3D particlization;                                    // make IS3D class and pass freezeout surface from vah (pinn is extraneous)
particlization.read_fo_surf_from_memory(vah.tau, vah.x, vah.y, vah.eta,
                    vah.dsigma_tau, vah.dsigma_x, vah.dsigma_y, vah.dsigma_eta,
                    vah.ux, vah.uy, vah.un,
                    vah.E, vah.T, vah.P,
                    vah.pixx, vah.pixy, vah.pixn, vah.piyy, vah.piyn, vah.pinn,
                    vah.Pi);

vah.free_freezeout_surface();                           // free surface memory in vah

particlization.run_particlization(0);                   // run particle sampler (argument = 0 uses surface from memory)

*/




