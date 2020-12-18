
#ifndef HYDROWRAPPER_H
#define HYDROWRAPPER_H

#include <fstream>
#include "FreezeoutSurface.h"
#include "Macros.h"

class HYDRO
{
    // this is a C++ wrapper for OSU hydro codes (cpu-vh, gpu-vh, cpu_vah, BEShydro, etc)
    // originally written by Derek Everett 2018

    private:

    public:
        HYDRO();
        ~HYDRO();

        // for holding initial energy density input vector from TRENTo:
        std::vector<double> trento_energy_density_profile;      // [GeV/fm^3] (note: no ghost cells included)


        // output freezeout surface vectors to pass to iS3D:
        std::vector<double> tau;                                // contravariant freezeout cell position
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> eta;                                // should change my freezeout surface class to match this (and units also)

        std::vector<double> dsigma_tau;                         // covariant surface normal vector
        std::vector<double> dsigma_x;
        std::vector<double> dsigma_y;
        std::vector<double> dsigma_eta;

        std::vector<double> E;                                  // energy density [GeV/fm^3]
        std::vector<double> T;                                  // temperature [GeV]
        std::vector<double> P;                                  // equilibrium pressure [GeV/fm^3]

        std::vector<double> ux;                                 // contravariant fluid velocity
        std::vector<double> uy;
        std::vector<double> un;
                                                                // contravariant shear stress pi^{\mu\nu}
        std::vector<double> pixx;                               // [GeV/fm^3]
        std::vector<double> pixy;                               // [GeV/fm^3]
        std::vector<double> pixn;                               // [GeV/fm^4]
        std::vector<double> piyy;                               // [GeV/fm^3]
        std::vector<double> piyn;                               // [GeV/fm^4]
        std::vector<double> pinn;                               // [GeV/fm^5]    (extraneous so leave it empty I guess)

        std::vector<double> Pi;                                 // bulk pressure [GeV/fm^3]


        // output freezeout surface to file (output/surface.dat)
    #ifndef JETSCAPE
        std::ofstream freezeout_surface_file;
    #endif

        void read_trento_energy_density_profile_from_memory(std::vector<double> energy_vector);

        void store_freezeout_surface(freezeout_surface surface);
        void free_freezeout_surface();

        void start_hydro(int argc, char **argv);                // run hydro simulation and save freezeout info

        void start_hydro_no_arguments();                        // no command line arguments (can't use auto_grid currently)
};

#endif




