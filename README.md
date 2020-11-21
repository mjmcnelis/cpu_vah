**CPU VAH (c) Mike McNelis and Dennis Bazow**

Created on 10/2015 by Dennis Bazow\
Last edited on 11/2020 by Mike McNelis

## Summary
A 3+1d relativistic hydrodynamic simulation for heavy-ion collisions

The C++ module is based off the hydrodynamic code [GPU VH](https://github.com/bazow/gpu-vh.git)
    
CPU VAH can run three hydrodynamic models with shear and bulk viscosity

    VAH = anisotropic hydrodynamics
    VH  = second-order viscous hydrodynamics (quasiparticle)
    VH2 = second-order viscous hydrodynamics (standard)


## References

If you use this code, please cite the following papers

     D. Bazow, U. Heinz and M. Strickland, Comput. Phys. Commun. 225 (2018) 92-113    
     M. McNelis, D. Bazow and U. Heinz, Phys. Rev. C 97, 054912 (2018)
     D.P. Bazow, Fluid dynamics for the anisotropically expanding quark gluon plasma, Ph.D. thesis (2017)


## Compile and run
To compile and run *n* hydro events with the default runtime `parameters`

    sh hydro.sh n  

Results from the simulation are stored in `output`\
Semi-analytic solutions (e.g. Bjorken and Gubser) are stored in `semi_analytic`

The above script clears the previous results prior to compiling. Alternatively, you can clear the results once by doing

    sh clear_results.sh
    make clean
    make
    for((i = 1; i <= n; i++))   # n = number of events (or jobs)
    do
        ./cpu_vah               # or submit your job
    done
    
This routine is often used by the job submissions in `scripts/`


## Compiler

By default, the code uses the `gcc` compiler (assumes no OpenMP support).\
Alternatively, you could use the `icpc` compiler (assumes OpenMP support)

You can switch out the Makefile by doing

    sh makefiles.sh icpc    (or gcc)
    

## Installation

You need to install the GSL libaries `-lgsl` and `-lgslcblas`\
The `icpc` compiler assumes the OpenMP library `-qopenmp` is installed

To use the `icpc` compiler, install `intel/19.0.3`\
To run the python3 scripts, install `python/3.6` 

To run CPU VAH with the Intel compiler on the Ohio State Supercomputer (OSC), do
    
    sh makefile.sh icpc
    module load python/3.6
    module load intel/19.0.3


## Source files

The C++ source and header files are located in `rhic`


## Runtime parameters

You can edit the runtime `parameters`

    hydro.properties
    initial.properties
    lattice.properties
    
The most important runtime parameters to adjust are

    run_hydro                           1 to output hydrodynamic quantites
                                        2 to construct the freezeout surface
                                        
    kinetic_theory_model                0 to run VH2 (comment macro ANISO_HYDRO)
                                        1 to run VH  (comment macro ANISO_HYDRO)
    
    ---------------------------------------------------------------------------------
      
    initial_condition_type              4 to run a smooth or fluctuating Trento initial condition (only for Pb+Pb sqrt(s) = 2.76 TeV)
                                        5 to read in energy density profile from file (tables/e_block.dat)
                                        6 to read in energy density profile from memory (JETSCAPE)
                                        
    trento_average_over_events          0 to run fluctuating Trento event (set initial_condition_type = 4)
                                        1 to run smooth Trento event and store initial energy density profile in tables/ (set initial_condition_type = 4)
                                        
    trento_number_of_average_events     number of events to event-average (e.g. 2000)
    
    trento_fixed_seed                   0 to set seed = clocktime
                                        x > 0 to fix seed = x (e.g. 1000)
    
    ---------------------------------------------------------------------------------
    
    number_of_points_x                  number of physical grid points along x-direction (and similarly for y, eta)
    
    lattice_spacing_x                   lattice spacing (and similarly for y, eta)
    
    resolve_nucleon_width               0 to use customized transverse lattice_spacing_x(y)
                                        1 to set lattice_spacing_x(y) = w/5 (w = nucleon width parameter)
                                                                            (transverse grid size is preserved)
    
    training_grid                       0 to use the customized spatial grid
                                        1 to use the 30 fm x 30 fm transverse grid (for training auto-grid) 
                                        
    auto_grid                           0 to use the customized spatial grid
                                        1 to automate the transverse grid lengths (need to run sh.predict_fireball_radius.sh)
    
    adaptive_time_step                  0 to set time step to fixed_time_step
                                        1 to use adaptive time step
    
    tau_coarse_factor                   number of time steps between freezeout finder calls (e.g. 2)
    

The Bayesian model parameters *P<sub>B</sub>* can be replaced during runtime. To generate *s* model parameter samples, go to `scripts/auto_grid` and do

    sh sample_model_parameters.sh s        
    
The model parameter samples are stored in `python/model_parameters`

Then run *n* hydro events with `model_parameters_p.dat`  (*p* ∈ [1, *s*])

    sh hydro.sh n p    

or simply run the executable

    ./cpu_vah p


## Macro parameters

You can edit the macro parameters in `rhic/include`

    Macros.h
    OpenMP.h
    
The most important macro parameters to adjust are
    
    ANISO_HYDRO         (run anisotropic hydro; comment to run second-order viscous hydro)
    BOOST_INVARIANT     (run 2+1d hydro; comment to run 3+1d hydro)
    CONFORMAL_EOS       (comment to use QCD equation of state)
    FREEZEOUT_VH        (write viscous hydrodymamic variables on freezeout surface for iS3D)
    JETSCAPE            (store the freezeout surface in memory; comment to output surface.dat)
    FLAGS               (print warnings during runtime)
    PRINT_PARAMETERS    (print the runtime parameters)
    OPENMP              (accelerate simulation with OpenMP)


## Tests

You can run the various tests performed in the code documentation paper in `scripts`

The files and jobs needed to run the tests are located in `tests` and `jobs`, respectively. You will need to copy `tables/example/e_block.dat` to the `initial_profile` directories in `tests`\ (only smooth 2d and 3d Trento energy density profiles are provided)

Results from the test runs are stored in `tests`


## Auto grid

The regression models used for the auto-grid are located in `tests/auto_grid/regression_models`. To launch them, go to `scripts/auto_grid` and do

        sh predict_fireball_radius.sh h s

where *h* ∈ [*vah*, *vh*, *vh2*] is hydrodynamic model you wish to run and *s* is the number of model parameter samples (see above). The script generates the `model_parameters` and `fireball_size_predictions` directories in `python`.

Then to use the auto grid, set `auto_grid = 1` in `parameters/lattice.properties` and do (see above)

    sh hydro.sh n p 

If you want to re-train the regression models, follow the steps in `scripts/auto_grid.README.txt`
