
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
		printf("Error: z = %lf is out of bounds\n", z);
		exit(-1);
	}
}


precision t_220(precision z, precision t)
{
	if(fabs(z) > delta && z > -1.0)
	{
		return (-1.0  +  (1.0 + z) * t) / z;
	}
	else if(fabs(z) <= delta)
	{
		precision z2  = z * z;
		precision z3 = z2 * z;
		precision z4 = z3 * z;
		precision z5 = z4 * z;
		precision z6 = z5 * z;
		precision z7 = z6 * z;
		precision z8 = z7 * z;

		return 0.6666666666666667 - 0.1333333333333333*z + 0.05714285714285716*z2 - 0.031746031746031744*z3 + 0.020202020202020193*z4 -
   		0.013986013986013984*z5 + 0.010256410256410262*z6 - 0.00784313725490196*z7 + 0.006191950464396287*z8;

	}
	else
	{
		printf("Error: z = %lf is out of bounds\n", z);
		exit(-1);
	}
}


precision t_201(precision z, precision t)
{
	if(fabs(z) > delta && z > -1.0)
	{
		return (1.0  +  (z - 1.0) * t) / z;
	}
	else if(fabs(z) <= delta)
	{
		precision z2  = z * z;
		precision z3 = z2 * z;
		precision z4 = z3 * z;
		precision z5 = z4 * z;
		precision z6 = z5 * z;
		precision z7 = z6 * z;
		precision z8 = z7 * z;

		return 1.3333333333333333 - 0.5333333333333333*z + 0.34285714285714286*z2 - 0.25396825396825395*z3 + 0.20202020202020202*z4 -
   		0.16783216783216784*z5 + 0.14358974358974358*z6 - 0.12549019607843137*z7 + 0.11145510835913312*z8;
	}
	else
	{
		printf("Error: z = %lf is out of bounds\n", z);
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
		printf("Error: z = %lf is out of bounds\n", z);
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
		printf("Error: z = %lf is out of bounds\n", z);
		exit(-1);
	}
}


precision t_202(precision z, precision t)
{
	if(fabs(z) > delta && z > -1.0)
	{
		return (3.0  +  z  +  (z - 3.0) * (1.0 + z) * t) / (z * z * (1.0 + z));
	}
	else if(fabs(z) <= delta)
	{
		return 1.0666666666666664 + z*(-1.3714285714285712 + z*(1.5238095238095237 + z*(-1.616161616161616 +
            z*(1.6783216783216781 + z*(-1.7230769230769227 + z*(1.756862745098039 + z*(-1.7832817337461297 + 1.8045112781954886*z)))))));
	}
	else
	{
		printf("Error: z = %lf is out of bounds\n", z);
		exit(-1);
	}
}


