#ifndef OPTICALGLAUBER_H_
#define OPTICALGLAUBER_H_

// smooth Glauber initial conditions
double woodsSaxonDistribution(double r, double A);

void Optical_Glauber_energy_density_transverse_profile(double * const __restrict__ energyDensityTransverse, int nx, int ny, double dx, double dy, void * initCondParams);

#endif