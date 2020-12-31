#ifndef MACROS_H_
#define MACROS_H_

#define ANISO_HYDRO
//#define BOOST_INVARIANT
//#define JETSCAPE

#define CONFORMAL_EOS

#ifndef CONFORMAL_EOS
	#define LATTICE_QCD
#endif

#ifdef ANISO_HYDRO
	#define PIMUNU

	#ifndef BOOST_INVARIANT
		#define WTZMU
	#endif
	#ifdef LATTICE_QCD
		#define B_FIELD
	#endif
#else
	#define PIMUNU

	#ifndef CONFORMAL_EOS
		#define PI
	#endif
#endif
#ifndef BOOST_INVARIANT
	#define VORTICITY
#endif

#define RANDOM_MODEL_PARAMETERS
#define PRINT_PARAMETERS

//#define FLAGS

//#define FREEZEOUT_SIZE
//#define FREEZEOUT_SLICE

//#define BENCHMARKS
#define ADAPTIVE_FILE

//#define MONITOR_TTAUMU
//#define MONITOR_REGULATIONS
#ifdef ANISO_HYDRO
	//#define MONITOR_PLPT
#ifdef LATTICE_QCD
	//#define MONITOR_B
#endif
#endif
#endif




