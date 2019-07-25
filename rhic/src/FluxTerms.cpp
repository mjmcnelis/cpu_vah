
#include <stdlib.h>
#include <math.h>
#include "../include/FluxTerms.h"
#include "../include/Precision.h"
#include "../include/DynamicalVariables.h"


inline precision sign(precision x)
{
	if(x > 0.0) return 1.0;
	else if(x < 0.0) return -1.0;
	else return 0.0;
}

inline precision minmod(precision x, precision y)
{
	return (sign(x) + sign(y)) * fmin(fabs(x), fabs(y)) / 2.0;
}


inline precision minmod3(precision x, precision y, precision z)
{
   return minmod(x, minmod(y, z));
}

precision approximate_derivative(precision qm, precision q, precision qp)
{
	return minmod3(THETA * (q - qm), (qp - qm) / 2.0, THETA * (qp - q));
}


void flux_terms(precision * const __restrict__ Hp, precision * const __restrict__ Hm, const precision * const __restrict__ q_data, const precision * const __restrict__ q1_data, const precision * const __restrict__ q2_data, const precision * const __restrict__ v_data, precision v)
{
	// neighbor spatial velocities
	precision vmm = v_data[0];
	precision vm  = v_data[1];
	precision vp  = v_data[2];
	precision vpp = v_data[3];

	precision dvp = approximate_derivative(v,   vp, vpp);
	precision dv  = approximate_derivative(vm,  v,  vp);
	precision dvm = approximate_derivative(vmm, vm, v);

	// extrapolated spatial velocities
	precision vRp = vp  -  0.5 * dvp;	// q^{+}_{i+1/2}
	precision vLp = v   +  0.5 * dv;	// q^{-}_{i+1/2}
	precision vRm = v   -  0.5 * dv;	// q^{+}_{i-1/2}
	precision vLm = vm  +  0.5 * dvm;	// q^{-}_{i-1/2}

	// local propagation speeds
	precision ap = fmax(fabs(vLp), fabs(vRp));
	precision am = fmax(fabs(vLm), fabs(vRm));

	int p = 0;

	// compute the flux terms Hp, Hm
	for(int n = 0; n < NUMBER_CONSERVED_VARIABLES; n++)
	{
		// neighbor conserved variables
		precision q   =  q_data[n];
		precision qm  = q1_data[p];
		precision qp  = q1_data[p + 1];
		precision qmm = q2_data[p];
		precision qpp = q2_data[p + 1];

		// conserved variable derivatives
		precision dqp = approximate_derivative(q,   qp, qpp);
		precision dq  = approximate_derivative(qm,  q,  qp);
		precision dqm = approximate_derivative(qmm, qm, q);

		// extrapolated conserved variables
		precision qRp = qp  -  0.5 * dqp;	// q^{+}_{i+1/2}	(R = +, p = i + 1/2)
		precision qLp = q   +  0.5 * dq;	// q^{-}_{i+1/2}	(L = -, p = i + 1/2)
		precision qRm = q   -  0.5 * dq;	// q^{+}_{i-1/2}	(R = -, m = i - 1/2)
		precision qLm = qm  +  0.5 * dqm;	// q^{-}_{i-1/2}	(L = -, m = i - 1/2)

		// neighbor fluxes
		precision Fm  = qm * vm;
		precision F   = q  * v;
		precision Fp  = qp * vp;

		// extrapolated fluxes (chain rule)
		precision FRp = Fp  -  0.5 * (qp * dvp  +  vp * dqp);	// F^{+}_{i+1/2}	(R = +, p = i + 1/2)	Eq. (63)
		precision FLp = F   +  0.5 * (q  * dv   +  v  * dq);	// F^{-}_{i+1/2}	(L = -, p = i + 1/2)	Eq. (64)
		precision FRm = F   -  0.5 * (q  * dv   +  v  * dq);	// F^{+}_{i-1/2}	(R = -, m = i - 1/2)	Eq. (65)
		precision FLm = Fm  +  0.5 * (qm * dvm  +  vm * dqm);	// F^{-}_{i-1/2}	(L = -, m = i - 1/2)	Eq. (66)

		// alternative formula (seems noisier overall but works better along x-axis)
		// precision FRp = qRp * vRp;
		// precision FLp = qLp * vLp;
		// precision FRm = qRm * vRm;
		// precision FLm = qLm * vLm;

		// Hp, Hm from Eq.(61)
		Hp[n] = (FRp  +  FLp  -  ap * (qRp - qLp)) / 2.0;
		Hm[n] = (FRm  +  FLm  -  am * (qRm - qLm)) / 2.0;

		p += 2;
	}
}


