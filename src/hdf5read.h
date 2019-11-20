/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, UC Davis
  Date: Oct 23, 2019

  Standalone cpp Poisson solver

*/

//****************** hdf5read.h **************************

#include <string>

#include "H5Cpp.h"
using namespace std;

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

int ReadHDF5File3(string, string, int, int, int, double*); // Rank 3 files 



