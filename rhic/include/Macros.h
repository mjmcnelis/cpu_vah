
#ifndef MACROS_H_
#define MACROS_H_

// macro parameters to fix equation of state and hydrodynamic variables at compile time

// should I remove CONFORMAL_EOS so I can switch? (worry about it later)
// it's only used for pt matching and pi

#define CONFORMAL_EOS			// use conformal equation of state (comment to use lattice QCD)

#define ANISO_HYDRO				// run anisotropic hydro (comment to run 2nd order viscous hydro)

//#define PIMUNU 					// include transverse shear stress (aniso) or shear stress (viscous hydro)

#ifdef ANISO_HYDRO
	#ifdef CONFORMAL_EOS
		#define PT_MATCHING 0
	#else
		#define PT_MATCHING 1	// may update to facilitate switching?
	#endif
	//#define WTZMU 			// include longitudinal momentum diffusion (comment for for 2+1d simulations)
#else
	#ifndef CONFORMAL_EOS
		#define PI
	#endif
#endif

//#define TEST_TTAUMU 			// test reproduction of t^{\tau\mu} (used in InferredVariables.cpp)

#endif
