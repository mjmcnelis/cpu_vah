CPU VAH (c) Mike McNelis and Dennis Bazow

Created on 10/2015 by Dennis Bazow\
Last edited on 11/2020 by Mike McNelis


## Summary
A 3+1d relativistic hydrodynamic simulation for heavy-ion collisions

The C++ module is based off the hydrodynamic code GPU VH\
The code can run three hydrodynamic models with shear and bulk viscosity

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

Simulation results are stored in `output`\
Bjorken and Gubser 0+1d solutions are stored in `semi_analytic`

The above script automatically clears the previous results\
Alternatively, the user can the clear results once by doing

    sh clear_results.sh
    make clean
    make
    for((i = 1; i <= N; i++))   # N = number of events
    do
        ./cpu_vah   # or submit your job
    done
    
This routine is often used by the job submissions in `scripts/`


## Runtime parameters

User can edit the runtime `parameters`

    hydro.properties
    initial.properties
    lattice.properties

The impact parameter `b` and Bayesian model parameters `P_B` can be replaced during runtime\
To generate `<s>` model parameter samples

    sh generate_model_parameter <s>         (scripts/auto_grid)
    
The model parameter samples are stored in `python/model_parameters`

Then run `<n>` hydro events with `model_parameters_<p>.dat  (<p> ∈ [1, <s>])`

    sh hydro.sh <n> <p>    


## Macro parameters


Source and header files are located in `rhic`

Macro parameters are located in rhic/include/Marcos.h

Option for OpenMP acceleration is located in rhic/include/OpenMP.h

Python files for training automated grid are located in python/
