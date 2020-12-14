
#ifndef EQUATIONOFSTATE_H_
#define EQUATIONOFSTATE_H_

#include "Precision.h"

precision energy_density_cutoff(precision e_min, precision e);

precision equilibrium_energy_density_new(precision T, precision conformal_prefactor);

class equation_of_state_new
{
	private:
		double e1;
		double conformal_prefactor;
		double T1, T2, T3, T4, T5, T6, T7, T8, T9, T10;
		double T11, T12, T13, T14, T15, T16, T17, T18, T19, T20;
		double T21, T22;
	public:
		precision T;

		equation_of_state_new(precision e1_in, precision conformal_prefactor_in);
		~equation_of_state_new();

		precision equilibrium_pressure();			// lattice qcd
		precision speed_of_sound_squared();

		precision z_quasi();						// quasiparticle
		precision mdmde_quasi();
		precision equilibrium_mean_field();

		precision beta_shear();						// viscosity / relaxation time
		precision beta_bulk();
};

#endif
