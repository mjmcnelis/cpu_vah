CPU VAH (c) Mike McNelis and Dennis Bazow

A 3+1d relativistic hydrodynamic simulation for heavy ion collisions


Code can run three different hydrodynamic models with shear and bulk viscosity:

    VAH = anisotropic hydrodynamics
    VH  = quasiparticle viscous hydrodynamics
    VH2 = standard viscous hydrodynamics


## Compile
To run a single hydro event with default `parameters`

    sh hydro.sh 1

## Usage 

To run iS3D

    ./iS3D

or

    sh runCPU.sh <num_threads>


Results are stored in `output`



Source and header files are located in `rhic`

Runtime parameters are located in parameters/ (contains default model parameters)

Macro parameters are located in rhic/include/Marcos.h

Option for OpenMP acceleration is located in rhic/include/OpenMP.h

Python files for training automated grid are located in python/
