
#include <iostream>
#include <string>
#include "H5Cpp.h"

const H5std_string FILE_NAME("test_cpp.h5");                                       // file name setting
const H5std_string DATASET_NAME("IntArray");                                       // example dataset
const int RANK = 2;                                                                // rank for data set. It means, the data set is 2-dimensional array
const int NX = 5;                                                                  // Nrows
const int NY = 6;                                                                  // Ncols

int main(void)
{
  int i, j;
  int data[NX][NY];
  for (j=0; j<NX; j++)                                                             // initialize the data array
    {
      for(i=0; i<NX; i++)
        data[j][i] = i + j;
    }
  try
    {
      H5::Exception::dontPrint();                                                  // turn of auto-printing in the case of error
      H5::H5File file(FILE_NAME, H5F_ACC_TRUNC);                                   // open file.
                                                                                   // check the properties for H5F_ACC_TRUNC 

      hsize_t dimsf[2];                                                            // dataset dimensions
      dimsf[0] = NX;
      dimsf[1] = NY;

      H5::DataSpace dataspace(RANK, dimsf);                                        // generate dataspace
      H5::IntType datatype(H5::PredType::NATIVE_INT);                              // set datatype as Integer type
      datatype.setOrder(H5T_ORDER_LE); 
      H5::DataSet dataset = file.createDataSet(DATASET_NAME, datatype, dataspace); // create dataset
      dataset.write(data, H5::PredType::NATIVE_INT);                               // write
    }

  catch(H5::FileIException error)                                                  // catch failure caused by the H5File operation
    {
      error.printError();
      return -1;
    }
  catch(H5::DataSetIException error)                                                // catch failure caused by the DataSet operations
    {
      error.printError();
      return -1;
    }
  catch(H5::DataSpaceIException error)                                             // catch failure caused by the DataSpace operations
    {
      error.printError();
      return -1;
    }
  catch(H5::DataTypeIException error)                                              // catch failure caused by the DataType operations
    {
      error.printError();
      return -1;
    }
  return 0;
}

/*
 * Local variables:
 * compile-command: "icpc -o exam_test exam_test.cpp -I/usr/include -I/usr/local/opt/szip/include -L/usr/local/Cellar/hdf5/1.8.16_1/lib /usr/local/Cellar/hdf5/1.8.16_1/lib/libhdf5_hl_cpp.a /usr/local/Cellar/hdf5/1.8.16_1/lib/libhdf5_cpp.a /usr/local/Cellar/hdf5/1.8.16_1/lib/libhdf5_hl.a /usr/local/Cellar/hdf5/1.8.16_1/lib/libhdf5.a -L/usr/lib -L/usr/local/opt/szip/lib -lsz -lz -ldl -lm"
 * End:
 */
