CPU VAH (c) Mike McNelis and Dennis Bazow

A 3+1d relativistic hydrodynamic simulation for heavy ion collisions

Code can run three hydrodynamic models with shear and bulk viscosity:

    VAH = anisotropic hydrodynamics
    VH  = quasiparticle viscous hydrodynamics
    VH2 = standard viscous hydrodynamics


## Compile and 
To compile run a `<n>` hydro events with default runtime `parameters`

    sh hydro.sh <n>

Or run `<n>` hydro events with sampled model parameters

Simulation results are stored in `output` 
Bjorken and Gubser mode and `semi_analytic` for 0+1d Bjorken and Gubser solutions)

## Runtime parameters

User can edit runtime `parameters`

        hydro.properties
        initial.properties
        lattice.properties





Source and header files are located in `rhic`

Runtime parameters are located in parameters/ (contains default model parameters)

Macro parameters are located in rhic/include/Marcos.h

Option for OpenMP acceleration is located in rhic/include/OpenMP.h

Python files for training automated grid are located in python/
