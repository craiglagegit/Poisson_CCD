/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, UC Davis
  Date: Oct 23, 2019

  Standalone cpp Poisson solver

*/

//****************** array3d.h **************************

#include <algorithm>
#include <stdio.h>       
#include <stdlib.h>      
#include <math.h>        
#include <globals.h>

using namespace std;

class Array3D //This packages the 3D data sets
{
 public:
  int nx, ny, nz;
  double xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dzp, nzexp, *x, *y, *z, *zp, *zpz, *zmz, *zw, *zplus, *zminus, *dzpdz, *data;
  double sensorthickness;
  int ChannelCkmin, ChannelStopCkmin, ChannelVkmin, ChannelStopVkmin;
  Array3D() {};
  Array3D(double, double, int, double, double, int, double, double, int, double, double);
  ~Array3D();
  double PyramidalKernel3D(double, double, double);
  double DataInterpolate3D(double, double, double);
  int XIndex(double);
  int YIndex(double);
  int ZIndex(double);  
  double ZP(double);
  double DZPDz(double);
  double D2ZPDz2(double);
  double Z(double);
  double ZPlus(double);
  double ZMinus(double);
};
