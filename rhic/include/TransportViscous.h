
#ifndef TRANSPORTVISCOUS_H_
#define TRANSPORTVISCOUS_H_

#include "Precision.h"
#include "Macros.h"
#include "Parameters.h"

class viscous_transport_coefficients
{
	private:
		double T;		// should try to change this to e
		double e;
		double p;
	public:

	#ifdef PIMUNU
		precision taupi_inverse;
		precision betapi;
		precision delta_pipi;
		precision tau_pipi;
		precision lambda_pibulkPi;
	#endif

	#ifdef PI
		precision taubulk_inverse;
		precision betabulk;
		precision lambda_bulkPipi;
		precision delta_bulkPibulkPi;
	#endif

		viscous_transport_coefficients(double T_in, double e_in, double p_in);
		~viscous_transport_coefficients();

		void compute_shear_transport_coefficients(precision etas);
		void compute_bulk_transport_coefficients(precision zetas, precision third_cs2);
};

#endif