
#ifndef MACROS_H_
#define MACROS_H_

// macro parameters to fix hydrodynamic variables at compile time

#define ANISO_HYDRO				// run anisotropic hydro (comment to run 2nd order viscous hydro)
#define BOOST_INVARIANT 		// run 2+1d hydro (comment to run 3+1d hydro)


#define BEST					// use the best eos (might be temporary after settling on one, or use parameter)


// I want to get rid of these two macros and use parameters instead
#define LATTICE_QCD				// use lattice qcd equation of state
//#define CONFORMAL_EOS			// use conformal equation of state (if both are defined, will transition eos)



#ifdef ANISO_HYDRO				// residual shear stress

	#define PIMUNU 				// transverse shear stress

	#ifndef BOOST_INVARIANT
		#define WTZMU 			// longitudinal momentum diffusion current
	#endif

#else							// viscous stress tensor

	#define PIMUNU				// shear stress

	#ifndef CONFORMAL_EOS
		#define PI 				// bulk pressure
	#endif

#endif


//#define TEST_TTAUMU 			// test reproduction of T^{\tau\mu} in InferredVariables.cpp


#endif