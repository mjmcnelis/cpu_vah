CPU VAH (c) Mike McNelis and Dennis Bazow

Date created:   10/2015      (Dennis Bazow)\
Last edited:    11/2020     (Mike McNelis)\


## Summary
A 3+1d relativistic hydrodynamic simulation for heavy-ion collisions

Code can run three hydrodynamic models with shear and bulk viscosity

    VAH = anisotropic hydrodynamics
    VH  = quasiparticle viscous hydrodynamics
    VH2 = standard viscous hydrodynamics



## Compile and run
To compile and run `<n>` hydro events with default runtime `parameters`

    sh hydro.sh <n>


Or run `<n>` hydro events with sampled model parameters `<p>` in `python/model_parameters`

To generate model parameter samples run in `scripts/auto_grid`

    sh hydro.sh <n> <p>         


Simulation results are stored in `output`\
Bjorken and Gubser 0+1d solutions are stored in `semi_analytic`



## Runtime parameters

User can edit the runtime `parameters`

    hydro.properties
    initial.properties
    lattice.properties

The impact parameter `b` and Bayesian model parameters `P_B` can be replaced during runtime

To generate `<s>` model parameter samples, run in `scripts/auto_grid`

    sh generateh <n> <p>         (note: model parameters are subset of runtime parameters)
    
The model parameter samples are stored in `python/model_parameters`

Then run `<n>` hydro events with sampled model parameters `<p> \on [1,s]`



## Macro parameters


Source and header files are located in `rhic`



Macro parameters are located in rhic/include/Marcos.h

Option for OpenMP acceleration is located in rhic/include/OpenMP.h

Python files for training automated grid are located in python/
