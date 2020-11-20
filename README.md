**CPU VAH** (c) Mike McNelis and Dennis Bazow

Created on 10/2015 by Dennis Bazow\
Last edited on 11/2020 by Mike McNelis

## Summary
A 3+1d relativistic hydrodynamic simulation for heavy-ion collisions

The C++ module is based off the hydrodynamic code GPU VH\
The code can run three hydrodynamic models with shear and bulk viscosity

    VAH = anisotropic hydrodynamics
    VH  = second-order viscous hydrodynamics (quasiparticle)
    VH2 = second-order viscous hydrodynamics (standard)


## References

If you use this code, please cite the following papers

     D. Bazow, U. Heinz and M. Strickland, Comput. Phys. Commun. 225 (2018) 92-113    
     M. McNelis, D. Bazow and U. Heinz, Phys. Rev. C 97, 054912 (2018)
     D. P. Bazow, Fluid dynamics for the anisotropically expanding quark gluon plasma, Ph.D. thesis (2017)


## Compile and run
To compile and run *n* hydro events with default runtime `parameters`

    sh hydro.sh n  

Results from the simulation are stored in `output`\
Semi-analytic solutions (e.g. Bjorken and Gubser) are stored in `semi_analytic`

The above script clears the previous results prior to compiling.\
Alternatively, the user can clear the results once by doing

    sh clear_results.sh
    make clean
    make
    for((i = 1; i <= N; i++))   # N = number of events
    do
        ./cpu_vah   # or submit your job
    done
    
This routine is often used by the job submissions in `scripts/`


## Runtime parameters

The user can edit the runtime `parameters`

    hydro.properties
    initial.properties
    lattice.properties
    
The most important runtime parameters to adjust are

    run_hydro = 1, 2               1 to output hydrodynamic quantites; 2 to construct the freezeout surface
    kinetic_theory_model = 0, 1    0 to run VH2; 1 to run VH (also need comment macro ANISO_HYDRO)
    

The Bayesian model parameters *P<sub>B</sub>* can be replaced during runtime.\
To generate *s* model parameter samples, go to `scripts/auto_grid` and do

    sh sample_model_parameters.sh s        
    
The model parameter samples are stored in `python/model_parameters`

Then run *n* hydro events with `model_parameters_p.dat`  (*p* ∈ [1, *s*])

    sh hydro.sh n p    

or simply run the executable

    ./cpu_vah p


## Macro parameters

The user can edit the macros in `rhic/include`

    Macros.h
    OpenMP.h
    
The most important macros to adjust are
    
    ANISO_HYDRO     (run anisotropic hydro, comment to run second-order viscous hydro)
    BOOST_INVARIANT (run 2+1d hydro, comment to run 3+1d hydro)
    CONFORMAL_EOS   (comment to use QCD equation of state)
    JETSCAPE        (store the freezeout surface)
    
    OPENMP          (accelerate simulation with OpenMP)

## Files

The C++ source and header files are located in `rhic`


## Tests



Python files for training automated grid are located in python/
