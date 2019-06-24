
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../include/Hypergeometric.h"


precision t_200(precision z, precision t)
{
	if(fabs(z) > delta && z > -1.0)
	{
		return 1.0  +  (1.0 + z) * t;			
	}
	else if(fabs(z) <= delta)
	{
		return 2.0 + z*(0.6666666666666667 + z*(-0.1333333333333333 + z*(0.05714285714285716 + z*(-0.031746031746031744 + z*(0.020202020202020193 + z*(-0.013986013986013984 + (0.010256410256410262 - 0.00784313725490196*z)*z))))));
	}
	else
	{
		printf("Conformal transport coefficients error: z = %lf is out of bounds\n", z);
		exit(-1);
	}
}


precision t_240(precision z, precision t)
{
	if(fabs(z) > delta && z > -1.0)
	{
		return (3.0  +  2.0 * z  -  3.0 * (1.0 + z) * t) / (z * z);
	}
	else if(fabs(z) <= delta)
	{
		return 0.4 + z*(-0.17142857142857149 + z*(0.09523809523809523 + z*(-0.06060606060606058 + z*(0.04195804195804195 + z*(-0.030769230769230785 + z*(0.023529411764705882 + (-0.01857585139318886 + 0.015037593984962405*z)*z))))));
	}
	else
	{
		printf("Hypergeometric function error: z = %lf is out of bounds\n", z);
		exit(-1);
	}
}


precision t_221(precision z, precision t)
{
	if(fabs(z) > delta && z > -1.0)
	{
		return (-3.0  +  (3.0 + z) * t) / (z * z);
	}
	else if(fabs(z) <= delta)
	{
		return 0.2666666666666668 + z*(-0.22857142857142854 + z*(0.19047619047619047 + z*(-0.1616161616161616 +  z*(0.13986013986013987 + z*(-0.12307692307692308 + z*(0.10980392156862744 + (-0.09907120743034055 + 0.09022556390977443*z)*z))))));
	}
	else
	{
		printf("Hypergeometric function error: z = %lf is out of bounds\n", z);
		exit(-1);
	}
}

