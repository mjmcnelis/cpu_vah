
#ifndef HYPERGEOMETRIC_H
#define HYPERGEOMETRIC_H

#include "Precision.h"

const precision delta = 0.01;	// window for taylor expansion of t_nlq functions around z = 0

// hypergeometric functions
precision t_200(precision z, precision t);
precision t_220(precision z, precision t);
precision t_201(precision z, precision t);

precision t_240(precision z, precision t);
precision t_221(precision z, precision t);
precision t_202(precision z, precision t);

#endif