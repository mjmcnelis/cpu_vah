#include <stdlib.h>
#include "../include/Macros.h"
#include "../include/EquationOfState.h"
#include "../include/TransportViscous.h"
#include "../include/Precision.h"
#include "../include/Parameters.h"

viscous_transport_coefficients::viscous_transport_coefficients(double T_in, double e_in, double p_in)
{
	T = T_in;
	e = e_in;
	p = p_in;
}


viscous_transport_coefficients::~viscous_transport_coefficients()
{

}


void viscous_transport_coefficients::compute_shear_transport_coefficients(precision etas)
{
#ifdef PIMUNU
	taupi_inverse = T / (5. * etas);
	betapi = (e + p) / 5.;
	delta_pipi = 4./3.;
	tau_pipi = 10./7;
	lambda_pibulkPi = 1.2;
#endif
}


void viscous_transport_coefficients::compute_bulk_transport_coefficients(precision zetas, precision third_cs2)
{
#ifdef PIMUNU
	taubulk_inverse	= 15. * third_cs2 * third_cs2 * T / zetas;		// double check this
	betabulk = 15. * third_cs2 * third_cs2 * (e + p);
	lambda_bulkPipi = 1.6 * third_cs2;
	delta_bulkPibulkPi = 2./3.;
#endif
}






