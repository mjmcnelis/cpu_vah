#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string>
#include <string.h>
#include <iostream>
#include <iomanip>
using namespace std;
#include <sstream>
#include <fstream>
#include "gauss_integration.hpp"
#include "thermal_integrands.hpp"
#include "qcd.hpp"

// temporary
const int a = 21;
const int gla_pts = 64;
double root_gla[a][gla_pts];
double weight_gla[a][gla_pts];


int load_gauss_laguerre_data()
{
  FILE *fp;

  stringstream laguerre_roots_weights;
  laguerre_roots_weights << "tables/gla_roots_weights_" << gla_pts << "_points.txt";

  if((fp = fopen(laguerre_roots_weights.str().c_str(), "r")) == NULL)
  {
     return 1;
  }
  for(int i = 0; i < a; i++)
  {
   for(int j = 0; j < gla_pts; j++)
   {
      if(fscanf(fp, "%i %lf %lf", &i, &root_gla[i][j], &weight_gla[i][j])!= 3)
      	{
        	printf("error loading roots/weights of Gauss-Laguerre Quadradture at %d %d\n", i, j);
    		return 1;
    	}
   }
  }
  fclose(fp);
  return 0;
}


int main()
{

	const int aJ = 3; // associated Laguerre polynomials a = 3 for I32
	const int aN = 1; // associated Laguerre polynomials a = 1 for I11
	double * pbar_rootJ = (double *)malloc(gla_pts * sizeof(double));
	double * pbar_weightJ = (double *)malloc(gla_pts * sizeof(double));
	double * pbar_rootN = (double *)malloc(gla_pts * sizeof(double));
	double * pbar_weightN = (double *)malloc(gla_pts * sizeof(double));
	// Load gauss laguerre roots-weights
	printf("Start loading gauss data...");
	int num_error;
	if((num_error = load_gauss_laguerre_data()) != 0)
	{
		fprintf(stderr, "Error loading gauss data (%d)!\n", num_error);
		return 1;
	}
	printf("done\n\n");



	// Set momentum bar roots-weights
	for(int i = 0; i < gla_pts; i++)
	{
		  pbar_rootJ[i]   = root_gla[aJ][i];
		pbar_weightJ[i] = weight_gla[aJ][i];

		  pbar_rootN[i]   = root_gla[aN][i];
		pbar_weightN[i] = weight_gla[aN][i];
	}



	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//                                                       ::
	//   	      QUASIPARTICLE RELAXATION MOMENTS           ::
	//                                                       ::
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	int n;
	double T;

	FILE *fp;
	fp = fopen("temperature.dat", "r");
	fscanf(fp, "%d", &n);

	const double nc = 3.;
	const double nf = 3.;
	const double g = pow(M_PI,4) * (4.*(nc*nc - 1.) + 7.*nc*nf) / 180.;  // degeneracy factor (g = 51.4103536012791)

	ofstream betapiplot, betabulkplot;
	betapiplot.open("results/betapi.dat", ios::out);
	betabulkplot.open("results/betabulk.dat", ios::out);

	for(int i = 0; i < n; i++)
	{
		fscanf(fp, "%lf", &T);

		double mbar = z_Quasiparticle(T);

		double factor_betapi = g * pow(T,4) / (30.*M_PI*M_PI);
		double betapi = factor_betapi * Gauss1D(I32_integrand, pbar_rootJ, pbar_weightJ, gla_pts, mbar, 0);

		double factor_I11 = g * pow(T,3) / (6.*M_PI*M_PI);
		double I11 = factor_I11 * Gauss1D(I11_integrand, pbar_rootN, pbar_weightN, gla_pts, mbar, 0);

		// only holds I11 (reconstruct rest of the formula in notebook)
		double betabulk = I11;

		betapiplot << setprecision(15) << T << "\t\t" << setprecision(15) << betapi << endl;
		betabulkplot << setprecision(15) << T << "\t\t" << setprecision(15) << betabulk << endl;
	}

	betapiplot.close();
	betabulkplot.close();
	fclose(fp);
	free(pbar_rootJ);
	free(pbar_weightJ);
	free(pbar_rootN);
	free(pbar_weightN);

	printf("\nDone\n");

	return 0;
}






