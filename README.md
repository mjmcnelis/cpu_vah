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


## Running the code
To compile and run the simulation with the default runtime `parameters`, do

    sh hydro.sh 1

Results from the simulation are stored in `output`.\
Semi-analytic solutions (e.g. Bjorken and Gubser) are stored in `semi_analytic`.


## Parameters

The runtime parameters are located in `parameters`

    hydro.properties
    initial.properties
    lattice.properties

You can replace the impact parameter *b* (if `initial_condition_type = 4`) and Bayesian model parameters *P<sub>B</sub>* during runtime. To generate *s* model parameter samples, go to `scripts/auto_grid` and do

    sh sample_model_parameters.sh s

The model parameter samples are stored in `python/model_parameters`.

To run the simulation with `model_parameters_p.dat`  (*p* ∈ [1, *s*]), do

    sh hydro.sh 1 p

or simply run the executable

    ./cpu_vah p

The macro parameters are located in `rhic/include`

    Macros.h
    OpenMP.h


## Tests

You can run the various tests performed in the code documentation paper in `scripts`. Some of these tests require access to multiple computing nodes.

The files and jobs needed to run the tests are located in `tables`, `tests` and `jobs`. You will need to edit the project number and email address in `jobs`.

Results from the test runs are stored in `tests`.
