/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, UC Davis
  Date: Oct 23, 2019

  Standalone cpp Poisson solver

*/
//****************** fileio.h **************************

#include <stdio.h> 
#include <string.h>
#include <string>
#include <fstream>

using namespace std;
double GetDoubleParam(string, string, double, int);
double* GetDoubleList(string, string, int, double*, int);
int GetIntParam(string, string, int, int);
int* GetIntList(string, string, int, int*, int);
string GetStringParam(string, string, string, int);
string ReadPar(string, string); 
