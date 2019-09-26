
#ifndef MACROS_H_
#define MACROS_H_

// macro parameters to fix equation of state and hydrodynamic variables at compile time

// should I remove CONFORMAL_EOS so I can switch? (worry about it later)
// it's only used for pt matching and pi

#define BOOST_INVARIANT 		// run 2+1d hydro (comment to run 3+1d hydro) (tie with Gubser)

//#define LATTICE_QCD				// use lattice qcd equation of state
#define CONFORMAL_EOS			// use conformal equation of state (if both are defined, will transition eos)

#define ANISO_HYDRO				// run anisotropic hydro (comment to run 2nd order viscous hydro)

#define PIMUNU 					// include transverse shear stress (aniso) or shear stress (viscous hydro)

#ifdef ANISO_HYDRO
	#ifdef CONFORMAL_EOS
		#define PT_MATCHING 0
	#else
		#define PT_MATCHING 1	// may update to facilitate switching?
	#endif
	#ifndef BOOST_INVARIANT
		#define WTZMU 			// include longitudinal momentum diffusion current
	#endif
#else
	#ifndef CONFORMAL_EOS
		#define PI
	#endif
#endif

//#define TEST_TTAUMU 			// test reproduction of t^{\tau\mu} (used in InferredVariables.cpp)

#endif
