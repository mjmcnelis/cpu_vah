
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>

#include "../include/EquationOfState.h"
#include "../include/TransportCoefficients.h"
#include "../include/AnisoGubser.h"

using namespace std;

double rho_function(double t, double r, double q)
{
	double rho = - asinh((1.0  -  q * q * (t * t  -  r * r)) / (2.0 * q * t));

	if(std::isnan(rho) || std::isinf(rho)) printf("Error: rho = %lf\n", rho);

	return rho;
}


double de_drho(double e, double pl, double rho, double etas)
{
	return (pl  -  3.*e) * tanh(rho);
}


double dpl_drho(double e, double pl, double rho, double etas)
{
	double T = effectiveTemperature(e);
	double taupiInv = 0.2 * T / etas;

	transport_coefficients aniso;
	aniso.compute_transport_coefficients(e, pl, (e - pl)/2.);

	return - taupiInv * (pl - e/3.)  -  (pl + aniso.I_240) * tanh(rho);
}


void gubser_rho_evolution(double * e_hat, double * pl_hat, double * rho_array, int rho_pts, double drho, double etas)
{
	// RK4 time evolution
	for(int i = 1; i < rho_pts; i++)
	{
		double rho = rho_array[i - 1];

		double ep  =  e_hat[i - 1];
		double plp = pl_hat[i - 1];

		double de1  = drho *  de_drho(ep, plp, rho, etas);
		double dpl1 = drho * dpl_drho(ep, plp, rho, etas);

		double de2  = drho *  de_drho(ep + de1/2., plp + dpl1/2., rho + drho/2., etas);
		double dpl2 = drho * dpl_drho(ep + de1/2., plp + dpl1/2., rho + drho/2., etas);

		double de3  = drho *  de_drho(ep + de2/2., plp + dpl2/2., rho + drho/2., etas);
		double dpl3 = drho * dpl_drho(ep + de2/2., plp + dpl2/2., rho + drho/2., etas);

		double de4  = drho *  de_drho(ep + de3, plp + dpl3, rho + drho, etas);
		double dpl4 = drho * dpl_drho(ep + de3, plp + dpl3, rho + drho, etas);

		e_hat[i]  = ep   +  (de1   +  2.0 * de2   +  2.0 * de3   +  de4)  / 6.0;
		pl_hat[i] = plp  +  (dpl1  +  2.0 * dpl2  +  2.0 * dpl3  +  dpl4) / 6.0;
	}

}

