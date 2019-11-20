/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, UC Davis
  Date: Oct 23, 2019

  Standalone cpp Poisson solver

*/

//****************** polygon.h **************************

#include <stdio.h>       
#include <stdlib.h>      
#include <algorithm>            // for min_element
#include <math.h>        

using namespace std;

#include <globals.h>

class Point
{
 public:
  double x,y,theta;
  void* owner;
  Point() {};
  Point(double, double, double);
};

class Polygon
{
 public:
  int npoints;
  bool sorted;
  double area;
  Point** pointlist;
  Polygon() {};
  Polygon(int);// Constructor
  ~Polygon();  //Destructor
  void AddPoint(Point*);
  void Sort();
  double Area();
  bool PointInside(Point*);
};

