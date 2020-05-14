
#ifndef MEMORY_H_
#define MEMORY_H_

double *** calloc_3d_array(double ***array, int dim1, int dim2, int dim3);
double **** calloc_4d_array(double ****array, int dim1, int dim2, int dim3, int dim4);
float ***** calloc_5d_array(float *****array, int dim1, int dim2, int dim3, int dim4, int dim5);

void free_3d_array(double ***array, int dim1, int dim2);
void free_4d_array(double ****array, int dim1, int dim2, int dim3);
void free_5d_array(float *****array, int dim1, int dim2, int dim3, int dim4);


#endif