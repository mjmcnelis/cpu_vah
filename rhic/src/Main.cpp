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

#include <iostream>
#include <cpuid.h>
#include <cstring>
#include "HydroWrapper.h"

void print_cpu()
{
    char CPUBrandString[0x40];
    unsigned int CPUInfo[4] = {0,0,0,0};

    __cpuid(0x80000000, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
    unsigned int nExIds = CPUInfo[0];

    memset(CPUBrandString, 0, sizeof(CPUBrandString));

    for (unsigned int i = 0x80000000; i <= nExIds; ++i)
    {
        __cpuid(i, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);

        if (i == 0x80000002)
            memcpy(CPUBrandString, CPUInfo, sizeof(CPUInfo));
        else if (i == 0x80000003)
            memcpy(CPUBrandString + 16, CPUInfo, sizeof(CPUInfo));
        else if (i == 0x80000004)
            memcpy(CPUBrandString + 32, CPUInfo, sizeof(CPUInfo));
    }

    std::cout << "\nCPU Type: " << CPUBrandString << "\n" << std::endl;
}

int main(int argc, char **argv)
{
  // default: cpu_vah can run as a stand-alone program
  // wrapper: HYDRO class can be instantiated in a larger program (JETSCAPE)

  print_cpu();

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

HYDRO vah;                                                          // make HYDRO class
vah.read_trento_energy_density_profile_from_memory(trento_energy);  // pass energy density vector from trento
vah.start_hydro_no_arguments();                                     // run hydro

IS3D particlization;                                                // make IS3D class and pass freezeout surface from vah (pinn is extraneous)

particlization.read_fo_surf_from_memory(vah.tau, vah.x, vah.y, vah.eta,
                    vah.dsigma_tau, vah.dsigma_x, vah.dsigma_y, vah.dsigma_eta,
                    vah.ux, vah.uy, vah.un,
                    vah.E, vah.T, vah.P,
                    vah.pixx, vah.pixy, vah.pixn, vah.piyy, vah.piyn, vah.pinn,
                    vah.Pi);

vah.free_freezeout_surface();                                       // free surface memory in vah

particlization.run_particlization(0);                               // run particle sampler (argument = 0 uses surface from memory)

*/




