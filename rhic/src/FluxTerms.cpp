#include <stdlib.h>
#include <math.h>
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"


inline precision sign(precision x)
{
	if(x > 0.) return 1.;
	else if(x < 0.) return -1.;
	else return 0.;
}


inline precision minmod(precision x, precision y)
{
	return (sign(x) + sign(y)) * fmin(fabs(x), fabs(y)) / 2.;
}


inline precision minmod3(precision x, precision y, precision z)
{
   return minmod(x, minmod(y, z));
}


precision minmod_derivative(precision qm, precision q, precision qp, precision Theta)
{
	return minmod3(Theta * (q - qm), (qp - qm) / 2., Theta * (qp - q));
}


precision compute_max_local_propagation_speed(const precision * const __restrict__ v_data, precision v, precision Theta)
{
	precision vmm = v_data[0];		// this is for computing the CFL condition
	precision vm  = v_data[1];
	precision vp  = v_data[2];
	precision vpp = v_data[3];

	precision dvp = minmod_derivative(v,   vp, vpp, Theta);
	precision dv  = minmod_derivative(vm,  v,  vp,  Theta);
	precision dvm = minmod_derivative(vmm, vm, v,   Theta);

	// extrapolated spatial speeds
	precision vRp = fabs(vp  -  dvp / 2.);	// |v^{+}_{i+1/2}|
	precision vLp = fabs(v   +  dv  / 2.);	// |v^{-}_{i+1/2}|
	precision vRm = fabs(v   -  dv  / 2.);	// |v^{+}_{i-1/2}|
	precision vLm = fabs(vm  +  dvm / 2.);	// |v^{-}_{i-1/2}|

	// local propagation speeds
	precision sp = fmax(vLp, vRp);	// s_{i+1/2}
	precision sm = fmax(vLm, vRm);	// s_{i-1/2}

	return fmax(sp, sm);			// max local speed for cell index i
}


void flux_terms(precision * const __restrict__ Hp, precision * const __restrict__ Hm, const precision * const __restrict__ q_data, const precision * const __restrict__ q1_data, const precision * const __restrict__ q2_data, const precision * const __restrict__ v_data, precision v, precision Theta)
{
	// neighbor spatial velocities
	precision vmm = v_data[0];
	precision vm  = v_data[1];
	precision vp  = v_data[2];
	precision vpp = v_data[3];

	precision dvp = minmod_derivative(v,   vp, vpp, Theta);
	precision dv  = minmod_derivative(vm,  v,  vp,  Theta);
	precision dvm = minmod_derivative(vmm, vm, v,   Theta);

	// extrapolated spatial velocities
	precision vRp = vp  -  dvp / 2.;	// v^{+}_{i+1/2}
	precision vLp = v   +  dv  / 2.;	// v^{-}_{i+1/2}
	precision vRm = v   -  dv  / 2.;	// v^{+}_{i-1/2}
	precision vLm = vm  +  dvm / 2.;	// v^{-}_{i-1/2}

	// local propagation speeds
	precision sp = fmax(fabs(vLp), fabs(vRp));	// s_{i+1/2}
	precision sm = fmax(fabs(vLm), fabs(vRm));	// s_{}

	int p = 0;

	// compute the flux terms Hp, Hm
	for(int n = 0; n < NUMBER_CONSERVED_VARIABLES; n++)
	{
		// neighbor dynamical variables
		precision q   =  q_data[n];
		precision qm  = q1_data[p];
		precision qp  = q1_data[p + 1];
		precision qmm = q2_data[p];
		precision qpp = q2_data[p + 1];

		// dynamical variable derivatives
		precision dqp = minmod_derivative(q,   qp, qpp, Theta);
		precision dq  = minmod_derivative(qm,  q,  qp,  Theta);
		precision dqm = minmod_derivative(qmm, qm, q,   Theta);

		// extrapolated dynamical variables
		precision qRp = qp  -  dqp / 2.;	// q^{+}_{i+1/2}
		precision qLp = q   +  dq  / 2.;	// q^{-}_{i+1/2}
		precision qRm = q   -  dq  / 2.;	// q^{+}_{i-1/2}
		precision qLm = qm  +  dqm / 2.;	// q^{-}_{i-1/2}

		// neighbor currents
		precision Fm  = qm * vm;
		precision F   = q  * v;
		precision Fp  = qp * vp;

		// extrapolated currents (chain rule)
		precision FRp = Fp  -  (qp * dvp  +  vp * dqp) / 2.;	// F^{+}_{i+1/2}
		precision FLp = F   +  (q  * dv   +  v  * dq)  / 2.;	// F^{-}_{i+1/2}
		precision FRm = F   -  (q  * dv   +  v  * dq)  / 2.;	// F^{+}_{i-1/2}
		precision FLm = Fm  +  (qm * dvm  +  vm * dqm) / 2.;	// F^{-}_{i-1/2}

		// numerical fluxes in KT algorithm
		Hp[n] = (FRp  +  FLp  -  sp * (qRp - qLp)) / 2.;		// H_{i+1/2}
		Hm[n] = (FRm  +  FLm  -  sm * (qRm - qLm)) / 2.;		// H_{i-1/2}

		p += 2;
	}
}


