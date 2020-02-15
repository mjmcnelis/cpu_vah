
#ifndef EQUATIONOFSTATE_H_
#define EQUATIONOFSTATE_H_

#include "Precision.h"

precision energy_density_cutoff(precision e_min, precision e);

precision equilibrium_energy_density(precision T, precision conformal_prefactor);


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


		// these need to be updated!!
		precision beta_shear();						// viscosity / relaxation time
		precision beta_bulk();
};


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

		// quasiparticle functions (let's not use them yet)


		// right now I'm replacing these with e dependence instead of T dependence
		// and testing them to make sure nothing has changed

		precision z_quasi(precision T);		// (temporary, should use e instead of T)
		precision mdmde_quasi();
		//precision mdmdT_quasi(precision T);															// is this still used anywhere?
		//precision equilibrium_kinetic_pressure(precision T, precision conformal_prefactor);			// is this still used anywhere?
		precision equilibrium_mean_field(precision T);

		// viscosity coefficients
		precision beta_shear(precision T, precision conformal_prefactor);
		precision beta_bulk(precision T);
};

#endif
