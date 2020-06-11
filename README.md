# cpu_vah
A 3+1d relativistic viscous hydrodynamic simulation for heavy ion collisions

Code can run three different hydrodynamic models with shear and bulk viscosity:

1) Viscous anisotropic hydrodynamics               (VAH)
2) 2nd order viscous hydrodynamics (quasiparticle) (VH)
3) 2nd order viscous hydrodynamics (standard)      (VH2)

To compile and run a single event with default model parameters:  sh hydro.sh 1

C++ source and header files are located in rhic/

Runtime parameters are located in parameters/ (contains default model parameters)

Macro parameters are located in rhic/include/Marcos.h

Option for OpenMP acceleration is located in rhic/include/OpenMP.h

Python files for training auto grid are located in python/
