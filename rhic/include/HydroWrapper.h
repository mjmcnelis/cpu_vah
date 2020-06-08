
#ifndef HYDROWRAPPER_H
#define HYDROWRAPPER_H

#include <stdlib.h>
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

        std::vector<double> initial_energy_density_vector;      // only need initial energy density (note: no ghost cells included)

        // Derek: do I need a surface vector here?

        void start_hydro(int argc, char **argv);                // run the hydro code until freezeout (renamed it start_hydro)


        // int run_hydro_no_cli();                              // Derek: do I need this version?

        //save the freezeout surface info
        //void save_fo_history();                               // Derek: what is this?

        // set initial_energy_density from a vector (useful for JETSCAPE)
        // note: argument should be [GeV/fm^3], then we convert to [fm^-4]

        void set_initial_energy_density_vector(std::vector<double> energy_vector);
};

#endif


