
#ifndef EQUATIONOFSTATE_H_
#define EQUATIONOFSTATE_H_

#include "Precision.h"

precision energy_density_cutoff(precision e_min, precision e);

precision equilibriumEnergyDensity(precision T, precision conformal_prefactor);


class equation_of_state
{
	private:
		double e1, e2, e3, e4, e5, e6, e7, e8, e9, e10;
		double e11, e12, e13, e14, e15, e16, e17, e18, e19, e20;
		double e21, e22, e23, e24;
	public:
		equation_of_state(precision e_in);
		~equation_of_state();

		precision equilibrium_pressure();
		precision speed_of_sound_squared();
		precision effective_temperature(precision conformal_prefactor);
};

#endif
