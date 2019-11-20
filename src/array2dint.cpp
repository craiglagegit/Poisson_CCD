/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, UC Davis
  Date: Oct 23, 2019

  Standalone cpp Poisson solver

*/

//****************** array2dint.cpp **************************

#include "array2dint.h"

Array2DInt::Array2DInt(double Xmin, double Xmax, int Nx, double Ymin, double Ymax, int Ny)
{
  int i, j;
  nx = Nx; ny = Ny;
  xmin = Xmin; xmax = Xmax; ymin = Ymin; ymax = Ymax;
  dx = (xmax - xmin) / (double) nx;
  dy = (ymax - ymin) / (double) ny;
  x = new double[nx]; y = new double[ny];
  data = new int[nx * ny];
  for (j=0; j<ny; j++)
    {
      y[j] = ymin + dy/2.0 + (double) j * dy;
    }
  for (i=0; i<nx; i++)
    {
      x[i] = xmin + dx/2.0 + (double) i * dx;	  
      for (j=0; j<ny; j++)
	{
	  data[i + j * nx] = 0.0;
	}
    }
}

Array2DInt::~Array2DInt()
{
  delete[] x;
  delete[] y;
  delete[] data;
}

int Array2DInt::XIndex(double x)
{
  return (int)((x - xmin) / dx);
}

int Array2DInt::YIndex(double y)
{
  return (int)((y - ymin) / dy);
}

