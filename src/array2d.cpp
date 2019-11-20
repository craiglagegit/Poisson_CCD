/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, UC Davis
  Date: Oct 23, 2019

  Standalone cpp Poisson solver

*/

//****************** array2d.cpp **************************

#include "array2d.h"

Array2D::Array2D(double Xmin, double Xmax, int Nx, double Ymin, double Ymax, int Ny)
{
  int i, j;
  nx = Nx; ny = Ny;
  xmin = Xmin; xmax = Xmax; ymin = Ymin; ymax = Ymax;
  dx = (xmax - xmin) / (double) nx;
  dy = (ymax - ymin) / (double) ny;
  x = new double[nx]; y = new double[ny];
  data = new double[nx * ny];
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

Array2D::~Array2D()
{
  delete[] x;
  delete[] y;
  delete[] data;
}

int Array2D::XIndex(double x)
{
  return (int)((x - xmin) / dx);
}

int Array2D::YIndex(double y)
{
  return (int)((y - ymin) / dy);
}

