CPU VAH (c) Mike McNelis and Dennis Bazow

Created on 10/2015 by Dennis Bazow\
Last edited on 11/2020 Mike McNelis


## Summary
A 3+1d relativistic hydrodynamic simulation for heavy-ion collisions

The C++ module is based off the hydrodynamic code GPU VH but runs three hydrodynamic models with shear and bulk viscosity

    VAH = anisotropic hydrodynamics
    VH  = quasiparticle second-order viscous hydrodynamics
    VH2 = standard second-order viscous hydrodynamics


## References

If you use this code, please cite the following papers

     D. Bazow, U. Heinz and M. Strickland, Comput. Phys. Commun. 225 (2018) 92-113    
     M. McNelis, D. Bazow and U. Heinz, Phys. Rev. C 97, 054912 (2018)
     D. P. Bazow, Fluid dynamics for the anisotropically expanding quark gluon plasma, Ph.D. thesis (2017)


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
