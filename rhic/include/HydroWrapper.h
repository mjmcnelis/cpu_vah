
#ifndef HYDROWRAPPER_H
#define HYDROWRAPPER_H

#include <stdlib.h>
#include "FreezeoutFinder.h"
#include <vector>

// #ifdef _OPENMP       // Derek: my code does not have openmp acceleration, is it required?
// #include <omp.h>
// #endif

using namespace std;

class HYDRO
{
    // this is a C++ wrapper for OSU hydro codes (cpu-vh, gpu-vh, cpu_vah, etc)
    // originally written by Derek Everett 2018

    private:

    public:
        HYDRO();
        ~HYDRO();

        // input initial energy density profile from TRENTo:

        std::vector<double> initial_energy_density_vector;      // initial energy profile [GeV/fm^3] (note: no ghost cells included)


        // output freezeout surface info to pass to iS3D:

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


        void load_initial_energy_density_vector(std::vector<double> energy_vector);

        void store_freezeout_surface(freezeout_surface surface);
        void free_freezeout_surface();

        void start_hydro(int argc, char **argv);                // run hydro simulation and save freezeout info

        // int run_hydro_no_cli();                              // Derek: do I need this version?
};

#endif




