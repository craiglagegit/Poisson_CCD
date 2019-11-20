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

//****************** hdf5write.cpp **************************

#include "hdf5write.h"

int WriteHDF5File3(string filename , string datasetname, int NX, int NY, int NZ, double* data) 
{
  int RANK = 3;
   // Try block to detect exceptions raised by any of the calls inside it
   try
   {
      /*
       * Turn off the auto-printing when failure occurs so that we can
       * handle the errors appropriately
       */
      Exception::dontPrint();

      /*
       * Create a new file using H5F_ACC_TRUNC access,
       * default file creation properties, and default file
       * access properties.
       */
      H5File file(filename, H5F_ACC_TRUNC );

      /*
       * Define the size of the array and create the data space for fixed
       * size dataset.
       */
      hsize_t     dimsf[3];              // dataset dimensions
      dimsf[0] = NX;
      dimsf[1] = NY;
      dimsf[2] = NZ;
      DataSpace dataspace( RANK, dimsf );

      /*
       * Define datatype for the data in the file.
       * We will store float numbers
       */
      FloatType datatype( PredType::NATIVE_DOUBLE );
      datatype.setOrder( H5T_ORDER_LE );

      /*
       * Create a new dataset within the file using defined dataspace and
       * datatype and default dataset creation properties.
       */
      DataSet dataset = file.createDataSet(datasetname, datatype, dataspace );

      /*
       * Write the data to the dataset using default memory space, file
       * space, and transfer properties.
       */
      dataset.write( data, PredType::NATIVE_DOUBLE );
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

   //printf("HDF5 data successfully written.\n");
   printf("HDF5 file %s successfully written.\n", filename.c_str());
   return 0;  // successfully terminated
}

int WriteHDF5File2(string filename , string datasetname, int NX, int NY, double* data) 
{
  int RANK = 2;
   // Try block to detect exceptions raised by any of the calls inside it
   try
   {
      /*
       * Turn off the auto-printing when failure occurs so that we can
       * handle the errors appropriately
       */
      Exception::dontPrint();

      /*
       * Create a new file using H5F_ACC_TRUNC access,
       * default file creation properties, and default file
       * access properties.
       */
      H5File file(filename, H5F_ACC_TRUNC );

      /*
       * Define the size of the array and create the data space for fixed
       * size dataset.
       */
      hsize_t     dimsf[2];              // dataset dimensions
      dimsf[0] = NX;
      dimsf[1] = NY;
      DataSpace dataspace( RANK, dimsf );

      /*
       * Define datatype for the data in the file.
       * We will store float numbers
       */
      FloatType datatype( PredType::NATIVE_DOUBLE );
      datatype.setOrder( H5T_ORDER_LE );

      /*
       * Create a new dataset within the file using defined dataspace and
       * datatype and default dataset creation properties.
       */
      DataSet dataset = file.createDataSet(datasetname, datatype, dataspace );

      /*
       * Write the data to the dataset using default memory space, file
       * space, and transfer properties.
       */
      dataset.write( data, PredType::NATIVE_DOUBLE );

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

   //printf("HDF5 data successfully written.\n");
   printf("HDF5 file %s successfully written.\n", filename.c_str());
   return 0;  // successfully terminated
}

int WriteHDF5File1(string filename , string datasetname, int NX, double* data) 
{
  int RANK = 1;
   // Try block to detect exceptions raised by any of the calls inside it
   try
   {
      /*
       * Turn off the auto-printing when failure occurs so that we can
       * handle the errors appropriately
       */
      Exception::dontPrint();

      /*
       * Create a new file using H5F_ACC_TRUNC access,
       * default file creation properties, and default file
       * access properties.
       */
      H5File file(filename, H5F_ACC_TRUNC );

      /*
       * Define the size of the array and create the data space for fixed
       * size dataset.
       */
      hsize_t     dimsf[1];              // dataset dimensions
      dimsf[0] = NX;
      DataSpace dataspace( RANK, dimsf );

      /*
       * Define datatype for the data in the file.
       * We will store float numbers
       */
      FloatType datatype( PredType::NATIVE_DOUBLE );
      datatype.setOrder( H5T_ORDER_LE );

      /*
       * Create a new dataset within the file using defined dataspace and
       * datatype and default dataset creation properties.
       */
      DataSet dataset = file.createDataSet(datasetname, datatype, dataspace );

      /*
       * Write the data to the dataset using default memory space, file
       * space, and transfer properties.
       */
      dataset.write( data, PredType::NATIVE_DOUBLE );
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

   //printf("HDF5 data successfully written.\n");
   printf("HDF5 file %s successfully written.\n", filename.c_str());
   return 0;  // successfully terminated
}

int WriteHDF5IntAttribute(string filename , string datasetname, string attr_name, int attr_dims, int* attr_data) 
{
  hsize_t dims[1] = { (hsize_t)attr_dims };

   // Try block to detect exceptions raised by any of the calls inside it
   try
   {
      
      // Turn off the auto-printing when failure occurs so that we can
      // handle the errors appropriately
       
      Exception::dontPrint();

      // Open an existing file and dataset.

      H5File file(filename, H5F_ACC_RDWR );
      DataSet dataset = file.openDataSet(datasetname);

      // Create the data space for the attribute.

      DataSpace attr_dataspace = DataSpace (1, dims );

      // Create a dataset attribute. 

      Attribute attribute = dataset.createAttribute(attr_name, PredType::STD_I32BE, 
                                                attr_dataspace, PropList::DEFAULT);
     
      // Write the attribute data. 

      attribute.write( PredType::NATIVE_INT, attr_data);

   }  // end of try block


   // catch failure caused by the H5File operations
   catch( DataSpaceIException error )
   {
      error.printErrorStack();
      return -1;
   }

   
   // catch failure caused by the H5File operations
   catch( AttributeIException error )
   {
      error.printErrorStack();
      return -1;
   }


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

   return 0;  // successfully terminated
}

int WriteHDF5DoubleAttribute(string filename , string datasetname, string attr_name, int attr_dims, double* attr_data) 
{
  hsize_t dims[1] = { (hsize_t)attr_dims };

   // Try block to detect exceptions raised by any of the calls inside it
   try
   {
      
      // Turn off the auto-printing when failure occurs so that we can
      // handle the errors appropriately
       
      Exception::dontPrint();

      // Open an existing file and dataset.

      H5File file(filename, H5F_ACC_RDWR );
      DataSet dataset = file.openDataSet(datasetname);

      // Create the data space for the attribute.

      DataSpace attr_dataspace = DataSpace (1, dims );

      // Create a dataset attribute. 

      Attribute attribute = dataset.createAttribute(attr_name, PredType::NATIVE_DOUBLE, 
                                                attr_dataspace, PropList::DEFAULT);
     
      // Write the attribute data. 

      attribute.write( PredType::NATIVE_DOUBLE, attr_data);

   }  // end of try block


   // catch failure caused by the H5File operations
   catch( DataSpaceIException error )
   {
      error.printErrorStack();
      return -1;
   }

   
   // catch failure caused by the H5File operations
   catch( AttributeIException error )
   {
      error.printErrorStack();
      return -1;
   }


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

   return 0;  // successfully terminated
}
