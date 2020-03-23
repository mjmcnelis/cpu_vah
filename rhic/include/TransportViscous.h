
#ifndef TRANSPORTVISCOUS_H_
#define TRANSPORTVISCOUS_H_

#include "Precision.h"
#include "Macros.h"
#include "Parameters.h"

class viscous_transport_coefficients
{
	private:
		double T;							// should try to change this to e
		double e;
		double p;
		int kinetic_theory_model;

		double T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11;
		double T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22;
	public:

	#ifdef PIMUNU							// shear transport coefficients
		precision taupi_inverse;
		precision betapi;
		precision delta_pipi;
		precision tau_pipi;
	#ifdef PI
		precision lambda_pibulkPi;
	#endif
	#endif

	#ifdef PI 								// bulk transport coefficients
		precision taubulk_inverse;
		precision betabulk;
		precision delta_bulkPibulkPi;
	#ifdef PIMUNU
		precision lambda_bulkPipi;
	#endif
	#endif

		viscous_transport_coefficients(double T_in, double e_in, double p_in, int kinetic_theory_model_in);
		~viscous_transport_coefficients();

		void compute_shear_transport_coefficients(precision etas);
		void compute_bulk_transport_coefficients(precision zetas, precision third_cs2);
};

#endif