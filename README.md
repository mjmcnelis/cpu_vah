**VAH (c) Mike McNelis and Dennis Bazow**

Created on 10/2015 by Dennis Bazow\
Last edited on 12/2020 by Mike McNelis

## Summary
A 3+1d relativistic hydrodynamic simulation for heavy-ion collisions.

The C++ module is based off the hydrodynamic code [GPU VH](https://github.com/bazow/gpu-vh.git).

VAH can run three hydrodynamic models with shear and bulk viscosity:

    VAH = anisotropic hydrodynamics
    VH  = second-order viscous hydrodynamics (quasiparticle)
    VH2 = second-order viscous hydrodynamics (standard)


## References

If you use this code, please cite the following papers:

     D. Bazow, U. Heinz and M. Strickland, Comput. Phys. Commun. 225 (2018) 92-113
     D.P. Bazow, Fluid dynamics for the anisotropically expanding quark gluon plasma, Ph.D. thesis (2017)
     M. McNelis, D. Bazow and U. Heinz, Phys. Rev. C 97, 054912 (2018)


## Setup
The default Makefile uses the gcc compiler but there is another one in `makefiles` that uses the icpc compiler (also assumes OpenMP support). To switch them out, do

    sh makefile.sh icpc     # or gcc

To set up the code on the Ohio Supercomputer Center, for example, login and do

    git clone https://github.com/mjmcnelis/cpu_vah.git
    cd cpu_vah && sh makefile.sh icpc
    module load intel/19.0.3
    module load python/3.6-conda5.2


## Running the code
To compile and run the hydrodynamic simulation once, do

    sh hydro.sh 1

Results from the simulation are stored in `output`.\
Semi-analytic solutions (e.g. Bjorken and Gubser) are stored in `semi_analytic`.


## Parameters

The runtime parameters are located in `parameters`

    hydro.properties
    initial.properties
    lattice.properties

The macro parameters are located in `rhic/include`

    Macros.h
    OpenMP.h


## Tests

You can run the various tests performed in the code documentation paper in `scripts`. Some of these tests require access to multiple computing nodes.

The files and jobs needed to run the tests are located in `tables`, `tests` and `jobs`. You will need to edit the project number and email address in `jobs`.

Results from the test runs are stored in `tests`.
