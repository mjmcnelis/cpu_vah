/*
 * FileIO.h
 *
 *  Created on: Oct 24, 2015
 *      Author: bazow
 */

#ifndef FILEIO_H_
#define FILEIO_H_

#include "DynamicalVariables.h"

void output(const PRECISION * const var, double t, const char *var_name, int nx, int ny, int nz, double dx, double dy, double dz);
void output_dynamical_variables(double t, int nx, int ny, int nz, double dx, double dy, double dz);

#endif
