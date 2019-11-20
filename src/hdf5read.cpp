/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and Copyright.html.  COPYING can be found at the root   *
 * of the source code distribution tree; Copyright.html can be found at the  *
 * root level of an installed copy of the electronic HDF5 document set and   *
 * is linked from the top-level documents page.  It can also be found at     *
 * http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have          *
 * access to either file, you may request a copy from help@hdfgroup.org.     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, UC Davis
  Date: Oct 23, 2019

  Standalone cpp Poisson solver

*/

//****************** hdf5read.cpp **************************

#include "hdf5read.h"

int ReadHDF5File3(string filename , string datasetname, int NX, int NY, int NZ, double* data_out) 
{
  int RANK_OUT = 3;
   /*
    * Try block to detect exceptions raised by any of the calls inside it
    */
   try
   {
      /*
       * Turn off the auto-printing when failure occurs so that we can
       * handle the errors appropriately
       */
      Exception::dontPrint();
      /*
       * Open the specified file and the specified dataset in the file.
       */
      H5File file( filename, H5F_ACC_RDONLY );
      DataSet dataset = file.openDataSet( datasetname );
      //H5T_class_t type_class = dataset.getTypeClass();

      /*
       * Get dataspace of the dataset.
       */
      DataSpace dataspace = dataset.getSpace();
      hsize_t     dimsm[3];              /* memory space dimensions */
      dimsm[0] = NX;
      dimsm[1] = NY;
      dimsm[2] = NZ;
      DataSpace memspace( RANK_OUT, dimsm );
      dataset.read( data_out, PredType::NATIVE_DOUBLE, memspace, dataspace );
   }  // end of try block
   // catch failure caused by the H5File operations
   catch( FileIException error )
   {
      error.printErrorStack();
      return -1;
   }
   // catch failure caused by the DataSet operations
   catch( DataSetIException error )
   {
      error.printErrorStack();
      return -1;
   }
   // catch failure caused by the DataSpace operations
   catch( DataSpaceIException error )
   {
      error.printErrorStack();
      return -1;
   }
   // catch failure caused by the DataSpace operations
   catch( DataTypeIException error )
   {
      error.printErrorStack();
      return -1;
   }
   printf("HDF5 file %s successfully read.\n", filename.c_str());
   return 0;  // successfully terminated
}

