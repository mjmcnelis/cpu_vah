
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include "../include/EquationOfState.h"
#include "../include/TransportCoefficients.h"
#include "../include/Parameters.h"
#include "../include/AnisoGubser.h"

using namespace std;

double rho_function(double t, double r, double q)
{
	return - asinh((1.  -  q * q * (t * t  -  r * r)) / (2. * q * t));
}


double de_drho(double e, double pl, double rho, hydro_parameters hydro)
{
	return (pl  -  3.*e) * tanh(rho);
}


double dpl_drho(double e, double pl, double rho, hydro_parameters hydro)
{
	double conformal_eos_prefactor = hydro.conformal_eos_prefactor;
	double T = effectiveTemperature(e, conformal_eos_prefactor);
	double etas = eta_over_s(T, hydro);
	double taupiInv = T / (5. * etas);

	transport_coefficients aniso;
	aniso.compute_transport_coefficients(e, pl, (e - pl)/2., conformal_eos_prefactor);

	precision zeta_LL = aniso.zeta_LL;

	return - taupiInv * (pl - e/3.)  -  (4. * pl + zeta_LL) * tanh(rho);
	
}


void gubser_rho_evolution(double * e_hat, double * pl_hat, double * rho_array, int rho_pts, double drho, hydro_parameters hydro)
{
	// RK4 time evolution
	for(int i = 1; i < rho_pts; i++)
	{
		double rho = rho_array[i - 1];

		double ep  =  e_hat[i - 1];
		double plp = pl_hat[i - 1];

		double de1  = drho *  de_drho(ep, plp, rho, hydro);
		double dpl1 = drho * dpl_drho(ep, plp, rho, hydro);

		double de2  = drho *  de_drho(ep + de1/2., plp + dpl1/2., rho + drho/2., hydro);
		double dpl2 = drho * dpl_drho(ep + de1/2., plp + dpl1/2., rho + drho/2., hydro);

		double de3  = drho *  de_drho(ep + de2/2., plp + dpl2/2., rho + drho/2., hydro);
		double dpl3 = drho * dpl_drho(ep + de2/2., plp + dpl2/2., rho + drho/2., hydro);

		double de4  = drho *  de_drho(ep + de3, plp + dpl3, rho + drho, hydro);
		double dpl4 = drho * dpl_drho(ep + de3, plp + dpl3, rho + drho, hydro);

		// RK4 update
		e_hat[i]  = ep   +  (de1   +  2. * de2   +  2. * de3   +  de4)  / 6.;
		pl_hat[i] = plp  +  (dpl1  +  2. * dpl2  +  2. * dpl3  +  dpl4) / 6.;
	}

}


