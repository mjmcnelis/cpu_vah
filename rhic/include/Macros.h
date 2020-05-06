
#ifndef MACROS_H_
#define MACROS_H_

// macro parameters to fix hydrodynamic variables at compile time

//#define ANISO_HYDRO				// run anisotropic hydro (comment to run 2nd order viscous hydro)
#define BOOST_INVARIANT 		// run 2+1d hydro (comment to run 3+1d hydro)


#define BEST					// use the best eos (might be temporary after settling on one, or use parameter)


// define only one equation of state
//#define CONFORMAL_EOS			// conformal equation of state
#define LATTICE_QCD				// lattice qcd equation of state


//#define E_CHECK					// separately evolve energy density w/ KT algorithm as a cross check
								// can be commented (doesn't affect rest of algorithm)


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

#define MONITOR_TTAUMU			// monitor violations of T^{\tau\mu} reproduction in InferredVariables.cpp

#endif