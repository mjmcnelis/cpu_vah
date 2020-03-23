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
	double * pbar_root_0 = (double *)malloc(gla_pts * sizeof(double));
	double * pbar_root_1 = (double *)malloc(gla_pts * sizeof(double));
	double * pbar_root_2 = (double *)malloc(gla_pts * sizeof(double));
	double * pbar_root_3 = (double *)malloc(gla_pts * sizeof(double));
	double * pbar_root_4 = (double *)malloc(gla_pts * sizeof(double));


	double * pbar_weight_0 = (double *)malloc(gla_pts * sizeof(double));
	double * pbar_weight_1 = (double *)malloc(gla_pts * sizeof(double));
	double * pbar_weight_2 = (double *)malloc(gla_pts * sizeof(double));
	double * pbar_weight_3 = (double *)malloc(gla_pts * sizeof(double));
	double * pbar_weight_4 = (double *)malloc(gla_pts * sizeof(double));



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
		pbar_root_0[i] = root_gla[0][i];
		pbar_root_1[i] = root_gla[1][i];
		pbar_root_2[i] = root_gla[2][i];
		pbar_root_3[i] = root_gla[3][i];
		pbar_root_4[i] = root_gla[4][i];

		pbar_weight_0[i] = weight_gla[0][i];
		pbar_weight_1[i] = weight_gla[1][i];
		pbar_weight_2[i] = weight_gla[2][i];
		pbar_weight_3[i] = weight_gla[3][i];
		pbar_weight_4[i] = weight_gla[4][i];
	}



	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//                                                       ::
	//   	      QUASIPARTICLE RELAXATION MOMENTS           ::
	//                                                       ::
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	int n;

	FILE *fp;
	fp = fopen("temperature.dat", "r");
	fscanf(fp, "%d", &n);

	const double nc = 3.;
	const double nf = 3.;
	const double g = pow(M_PI,4) * (4.*(nc*nc - 1.) + 7.*nc*nf) / 180.;  // degeneracy factor (g = 51.4103536012791)

	ofstream betapi_data;

	ofstream cE_bar_data;
	ofstream cP_bar_data;
	ofstream cpi_bar_data;

	ofstream I00_data;
	ofstream I01_data;

	ofstream I11_data;

	ofstream I21_data;
	ofstream I22_data;

	ofstream I40_data;
	ofstream I41_data;
	ofstream I42_data;

	betapi_data.open("results/betapi.dat", ios::out);

	cE_bar_data.open("results/cE_bar.dat", ios::out);
	cP_bar_data.open("results/cP_bar.dat", ios::out);
	cpi_bar_data.open("results/cpi_bar.dat", ios::out);

	I00_data.open("results/I00.dat", ios::out);
	I01_data.open("results/I01.dat", ios::out);

	I11_data.open("results/I11.dat", ios::out);

	I21_data.open("results/I21.dat", ios::out);
	I22_data.open("results/I22.dat", ios::out);

	I40_data.open("results/I40.dat", ios::out);
	I41_data.open("results/I41.dat", ios::out);
	I42_data.open("results/I42.dat", ios::out);


	for(int i = 0; i < n; i++)
	{
		double T;
		fscanf(fp, "%lf", &T);

		double T2 = T * T;
		double T3 = T2 * T;
		double T4 = T3 * T;
		double T5 = T4 * T;
		double T6 = T5 * T;

		double mbar = z_Quasiparticle(T);


		double I00 = g * T2 / (2. * M_PI * M_PI)  * Gauss1D(I00_integrand, pbar_root_0, pbar_weight_0, gla_pts, mbar);
		double I01 = g * T2 / (6. * M_PI * M_PI)  * Gauss1D(I01_integrand, pbar_root_0, pbar_weight_0, gla_pts, mbar);


		double I11 = g * T3 / (6. * M_PI * M_PI)  * Gauss1D(I11_integrand, pbar_root_1, pbar_weight_1, gla_pts, mbar);


		double I21 = g * T4 / (6. * M_PI * M_PI)  * Gauss1D(I21_integrand, pbar_root_2, pbar_weight_2, gla_pts, mbar);
		double I22 = g * T4 / (30. * M_PI * M_PI) * Gauss1D(I22_integrand, pbar_root_2, pbar_weight_2, gla_pts, mbar);


		double I32 = g * T5 / (30. * M_PI * M_PI) * Gauss1D(I32_integrand, pbar_root_3, pbar_weight_3, gla_pts, mbar);


		double I40 = g * T6 / (2. * M_PI * M_PI)  * Gauss1D(I40_integrand, pbar_root_4, pbar_weight_4, gla_pts, mbar);
		double I41 = g * T6 / (6. * M_PI * M_PI)  * Gauss1D(I41_integrand, pbar_root_4, pbar_weight_4, gla_pts, mbar);
		double I42 = g * T6 / (30. * M_PI * M_PI) * Gauss1D(I42_integrand, pbar_root_4, pbar_weight_4, gla_pts, mbar);


		double cE_bar = - I41 / (5./3. * I40 * I42  -  I41 * I41);
		double cP_bar =   I40 / (5./3. * I40 * I42  -  I41 * I41);
		double cpi_bar = 1. / I42;

		betapi_data << setprecision(15) << T << "\t\t" << setprecision(15) << I32 / T << endl;

		cE_bar_data << setprecision(15) << T << "\t\t" << setprecision(15) << cE_bar << endl;
		cP_bar_data << setprecision(15) << T << "\t\t" << setprecision(15) << cP_bar << endl;
		cpi_bar_data << setprecision(15) << T << "\t\t" << setprecision(15) << cpi_bar << endl;

		I00_data 	<< setprecision(15) << T << "\t\t" << setprecision(15) << I00 << endl;
		I01_data 	<< setprecision(15) << T << "\t\t" << setprecision(15) << I01 << endl;

		I11_data 	<< setprecision(15) << T << "\t\t" << setprecision(15) << I11 << endl;

		I21_data 	<< setprecision(15) << T << "\t\t" << setprecision(15) << I21 << endl;
		I22_data 	<< setprecision(15) << T << "\t\t" << setprecision(15) << I22 << endl;

		I40_data 	<< setprecision(15) << T << "\t\t" << setprecision(15) << I40 << endl;
		I41_data 	<< setprecision(15) << T << "\t\t" << setprecision(15) << I41 << endl;
		I42_data 	<< setprecision(15) << T << "\t\t" << setprecision(15) << I42 << endl;
	}

	betapi_data.close();

	cE_bar_data.close();
	cP_bar_data.close();
	cpi_bar_data.close();

	I00_data.close();
	I01_data.close();

	I11_data.close();

	I21_data.close();
	I22_data.close();

	I40_data.close();
	I41_data.close();
	I42_data.close();

	fclose(fp);

	printf("\nDone\n");

	return 0;
}






