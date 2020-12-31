#ifndef MACROS_H_
#define MACROS_H_

// simulation mode
#define ANISO_HYDRO                        // run anisotropic hydro (comment to run 2nd order viscous hydro)
#define BOOST_INVARIANT                    // run 2+1d hydro (comment to run 3+1d hydro)
//#define JETSCAPE                         // store freezeout surface in HYDRO wrapper (comment to write freezeout surface to file)
                                           // todo: in case of fatal error (e.g. nan), return empty surface in run_hydro() w/o exiting program
                                           // places program can crash: initial conditions (read in incorrectly)
                                           //                           reconstruct inferred (nan) and aniso variables (nan, gsl error)
                                           //                           Cornelius
                                           //                           exceed max_time_steps


// equation of state
//#define CONFORMAL_EOS                    // use conformal equation of state (comment to use QCD)

#ifndef CONFORMAL_EOS
	#define LATTICE_QCD                // QCD equation of state (leave defined)
#endif


// hydrodynamic variables
#ifdef ANISO_HYDRO
	#define PIMUNU                     // switch to include transverse shear stress (can comment)

	#ifndef BOOST_INVARIANT
		#define WTZMU              // switch to include longitudinal momentum diffusion current (can comment)
	#endif
	#ifdef LATTICE_QCD
		#define B_FIELD            // mean field (leave defined)
	#endif
#else
	#define PIMUNU                     // switch to include standard shear stress (can comment)

	#ifndef CONFORMAL_EOS
		#define PI                 // switch to include bulk pressure (can comment)
	#endif
#endif
#ifndef BOOST_INVARIANT
	#define VORTICITY                  // switch to include vorticity terms in relaxation equations (can comment)
#endif


// parameters
#define RANDOM_MODEL_PARAMETERS            // switch to use python/model_parameters/model_parameters_p.dat if do ./cpu_vah p
                                           // note: parameter sample {b, P_B} is a subset of all runtime parameters
//#define PRINT_PARAMETERS                 // switch to print parameters during runtime


//#define FLAGS                            // switch to print warnings during runtime for debugging (comment for production runs)


// output (comment for production runs)
//#define FREEZEOUT_SIZE                   // switch to output maximum transverse radius of freezeout surface
//#define FREEZEOUT_SLICE                  // switch to output tau-x (and tau-eta) slice up to tau = 17 fm/c

//#define BENCHMARKS                       // switch to output runtime benchmarks
//#define ADAPTIVE_FILE                    // switch to output adaptive time step

//#define MONITOR_TTAUMU                   // switch to output T^{\tau\mu} violations
//#define MONITOR_REGULATIONS              // switch to output residual shear (or shear + bulk) regulation
#ifdef ANISO_HYDRO
	//#define MONITOR_PLPT             // switch to output pl, pt regulation (i.e. enforce pl, pt > 0 at edge of grid)
#ifdef LATTICE_QCD
	//#define MONITOR_B                // switch to output b regulation
#endif
#endif
#endif




