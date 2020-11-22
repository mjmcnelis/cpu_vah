
#ifndef MACROS_H_
#define MACROS_H_

#define ANISO_HYDRO				// run anisotropic hydro (comment to run 2nd order viscous hydro)

//#define BOOST_INVARIANT 		// run 2+1d hydro (comment to run 3+1d hydro)

#define RANDOM_MODEL_PARAMETERS // option to use python/random_model_parameters/model_parameters_x.dat (x = 1st command line argument)
								// comment to use the fixed values in parameters/
								// note: model parameters is a subset of all parameters

#define JETSCAPE 				// store freezeout surface in vector array for JETSCAPE (no surface.dat)
								// also need to define FREEZEOUT_VH to run JETSCAPE mode

#define FREEZEOUT_VH			// option to convert vah to vh variables on freezeout surface

#define BEST					// use the best eos (I don't think this is used anymore)

//#define CONFORMAL_EOS			// conformal equation of state

#ifndef CONFORMAL_EOS
	#define LATTICE_QCD			// lattice qcd equation of state
#endif


//#define E_CHECK				// separately evolve energy density w/ KT algorithm as a cross check
								// can be commented (still have issue with piperp evolution when turned on)


#ifdef ANISO_HYDRO				// residual shear stress

	#define PIMUNU 				// transverse shear stress

	#ifndef BOOST_INVARIANT
		#define WTZMU 			// longitudinal momentum diffusion current
	#endif

	#ifdef LATTICE_QCD
		#define B_FIELD 		// mean field
	#endif

#else							// viscous stress tensor

	#define PIMUNU				// shear stress

	#ifndef CONFORMAL_EOS
		#define PI 				// bulk pressure
	#endif

#endif

#ifndef BOOST_INVARIANT
	//#define VORTICITY			// include vorticity terms in relaxation equations
#endif



#define PRINT_HYDRO				// option to print current hydro info at center of grid (comment for no printing)
#define PRINT_PERIOD 10			// number of steps per print for run_hydro = 3 (leave this defined!)

#define PRINT_PARAMETERS		// option to print parameters

//#define FLAGS                    // option to print warnings that could repeat during runtime (doesn't include nans or exit(-1) errors)

//#define FREEZEOUT_SIZE			// output maximum radius of freezeout surface

#define FREEZEOUT_SLICE			// output tau-x slice (y = eta = 0) and tau-eta slice (x = y = 0, only for 3+1d) of freezeout surface

#define BENCHMARKS				// output benchmark data (e.g. hydro run time)

#define ADAPTIVE_FILE			// output adaptive time steps

//#define MONITOR_TTAUMU			// output violations of T^{\tau\mu} reproduction in InferredVariables.cpp

#define MONITOR_REGULATIONS	// output viscous or aniso regulations

#ifdef ANISO_HYDRO
	//#define MONITOR_PLPT		// output pl,pt regulation

#ifdef LATTICE_QCD
	#define MONITOR_B			// output b regulation
#endif
#endif

#endif






