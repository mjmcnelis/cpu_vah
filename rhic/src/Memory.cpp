
#include <stdlib.h> 


double *** calloc_3d_array(double ***array, int dim1, int dim2, int dim3)
{
  array = (double ***)calloc(dim1, sizeof(double **));

  for(int i = 0; i < dim1; i++)
  {
    array[i] = (double **)calloc(dim2, sizeof(double *));

    for(int j = 0; j < dim2; j++)
    {
      array[i][j] = (double *)calloc(dim3, sizeof(double));
    }
  }
  return array;
}


double **** calloc_4d_array(double ****array, int dim1, int dim2, int dim3, int dim4)
{
  array = (double ****)calloc(dim1, sizeof(double ***));

  for(int i = 0; i < dim1; i++)
  {
    array[i] = (double ***)calloc(dim2, sizeof(double **));

    for(int j = 0; j < dim2; j++)
    {
      array[i][j] = (double **)calloc(dim3, sizeof(double *));

      for(int k = 0; k < dim3; k++)
      {
        array[i][j][k] = (double *)calloc(dim4, sizeof(double));
      }
    }
  }

  return array;
}


double ***** calloc_5d_array(double *****array, int dim1, int dim2, int dim3, int dim4, int dim5)
{
  array = (double *****)calloc(dim1, sizeof(double ****));

  for(int i = 0; i < dim1; i++)
  {
    array[i] = (double ****)calloc(dim2, sizeof(double ***));

    for (int j = 0; j < dim2; j++)
    {
      array[i][j] = (double ***)calloc(dim3, sizeof(double **));

      for(int k = 0; k < dim3; k++)
      {
        array[i][j][k] = (double **)calloc(dim4, sizeof(double *));

        for (int l = 0; l < dim4; l++)
        {
          array[i][j][k][l] = (double *)calloc(dim5, sizeof(double));
        }
      }
    }
  }

  return array;
}



void free_3d_array(double ***array, int dim1, int dim2)
{
  for(int i = 0; i < dim1; i++)
  {
    for(int j = 0; j < dim2; j++)
    {
      free(array[i][j]);
    }

    free(array[i]);
  }

  free(array);
}


void free_4d_array(double ****array, int dim1, int dim2, int dim3)
{
  for(int i = 0; i < dim1; i++)
  {
    for(int j = 0; j < dim2; j++)
    {
      for(int k = 0; k < dim3; k++)
      {
        free(array[i][j][k]);
      }

      free(array[i][j]);
    }

    free(array[i]);
  }

  free(array);
}



void free_5d_array(double *****array, int dim1, int dim2, int dim3, int dim4)
{
  for(int i = 0; i < dim1; i++)
  {
    for(int j = 0; j < dim2; j++)
    {
      for(int k = 0; k < dim3; k++)
      {
        for(int l = 0; l < dim4; l++)
        {
          free(array[i][j][k][l]);
        }

        free(array[i][j][k]);
      }

      free(array[i][j]);
    }

    free(array[i]);
  }

  free(array);
}





