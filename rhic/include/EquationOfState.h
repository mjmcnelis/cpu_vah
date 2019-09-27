
#ifndef EQUATIONOFSTATE_H_
#define EQUATIONOFSTATE_H_

#include "Precision.h"

precision energy_density_cutoff(precision e_min, precision e);

precision equilibrium_energy_density(precision T, precision conformal_prefactor);


class equation_of_state
{
	private:
		double e1, e2, e3, e4, e5, e6, e7, e8, e9, e10;
		double e11, e12, e13, e14, e15, e16, e17, e18, e19, e20;
		double e21, e22, e23, e24, e25, e26;
	public:
		equation_of_state(precision e_in);
		~equation_of_state();

		precision equilibrium_pressure();
		precision speed_of_sound_squared();
		precision effective_temperature(precision conformal_prefactor);

		// quasiparticle functions (let's not use them yet. just slap them here for future reference)

		precision z_quasi(precision T);		// (temporary, should use e instead of T)
		precision mdmde_quasi();
		precision mdmdT_quasi(precision T);	// is this still used anywhere?
		precision equilibrium_kinetic_pressure(precision T, precision conformal_prefactor);
		precision equilibrium_mean_field(precision T);

		// viscosity coefficients
		precision beta_shear(precision T, precision conformal_prefactor);
		precision beta_bulk(precision T);
};

#endif
