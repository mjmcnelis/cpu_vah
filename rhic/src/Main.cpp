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

  //vah.initialize_ed_from_vector(init_e);      // pass the initial energy density vector (Derek: how would this work?)

                                                // replaced int run_hydro -> void start_hydro
  vah.start_hydro(argc, argv);                  // runs hydro code (Derek: do I need a version with no command line arguments?)

  // Derek: anything else after? (freezeout info)
  // Derek: does HYDRO class need a surface vector to pass to iS3D? (my freezeout_finder class stores a surface vector)

  return 0;
}
