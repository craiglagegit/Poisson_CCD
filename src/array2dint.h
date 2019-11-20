/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, UC Davis
  Date: Oct 23, 2019

  Standalone cpp Poisson solver

*/

//****************** array2dint.h **************************

#include <stdio.h>       
#include <stdlib.h>      
#include <math.h>        

#include <globals.h>


class Array2DInt //This packages the 2D data sets
{
 public:
  int nx, ny, *data;
  double xmin, xmax, ymin, ymax, dx, dy, dz, *x, *y;
  Array2DInt() {};
  Array2DInt(double, double, int, double, double, int);
  ~Array2DInt();
  int XIndex(double);
  int YIndex(double);

};
