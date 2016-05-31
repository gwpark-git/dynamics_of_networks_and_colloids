
/*
  #################################################################
  ##   Private MATRIX Library ver 0.3                            ##
  ##   Made by "Gun Woo Park", Email : gunwoo.park@unina.it      ##
  ##   Fellow of SUPOLEN, Marie Curie Actions                    ##
  ##   DICMaPI, University of Naples Federico II                 ##
  #################################################################
*/

// all the integer is converted to MKL_LONG

// @ denoted for tested function
// Lots of functions are based on Call by Reference.
// Notice that for matrix operation, we have to put &.
// Vector can use by rows = 1.

// Note that I seperated library for linear algebra and matrix.
// If some error occurs during compiling, please check to include linear algebra code.

#ifndef MATRIX_H
#define MATRIX_H
#define TRUE 1
#define FALSE 0
#define BIT 64
#include <iostream>
#include <fstream>
#include <mkl.h>
#include <math.h>
//#include <omp.h>

extern "C" {
#include <gsl/gsl_sort_double.h>
}

class MATRIX
{
 public:
  MKL_LONG rows;
  MKL_LONG cols;
  MKL_LONG size;
  double *data;  // 1-dimensional array
  bool VECTOR; 
  bool INITIALIZATION;
  bool DIAGONALIZATION;
  double *eigen_value;
  MKL_LONG DIAGONALIZATION_INITIAL()                     //@ 
  {
    /* eigen_value = new double [rows]; */
    /* eigen_value = (double*) mkl_malloc(rows*sizeof(double), BIT); */
    eigen_value = new double [rows];
    DIAGONALIZATION = TRUE;
    return 0;
  }
  MKL_LONG print_eigen();
  MKL_LONG print_eigen(MKL_LONG n_ele);                                //@
  // Public Member Function
  MKL_LONG p_self();                                     //@
  MKL_LONG initial();
  MKL_LONG initial(MKL_LONG N_r, MKL_LONG N_c);                    //@
  MKL_LONG initial(MKL_LONG N_r, MKL_LONG N_c, double x);
  MKL_LONG print();                                      //@                     
  MKL_LONG print(MKL_LONG n_row, MKL_LONG n_col);                  //@

  MKL_LONG fprint_skip(std::ofstream& file, MKL_LONG N_skip);
  MKL_LONG fprint_skip_transpose(std::ofstream& file, MKL_LONG N_skip);
  MKL_LONG fprint_LONG_skip(std::ofstream& file, MKL_LONG N_skip);
  MKL_LONG fprint_LONG_skip_transpose(std::ofstream& file, MKL_LONG N_skip);
  MKL_LONG fprint_LONG_skip_transpose_LIMROWS(std::ofstream& file, MKL_LONG N_skip, MKL_LONG N_lim_rows);
  MKL_LONG fprint(std::ofstream& file);                             // ???
  MKL_LONG fprint_transpose(std::ofstream& file);                             // ???
  MKL_LONG fprint_LONG(std::ofstream& file);                             // ???
  MKL_LONG fprint_LONG_transpose(std::ofstream& file);                             // ??? 
  MKL_LONG fprint_row(std::ofstream& file, MKL_LONG given_row);
  MKL_LONG fprint_out_skip(std::ofstream& file, MKL_LONG N_skip);
  MKL_LONG fprint_out(std::ofstream& file);                             // ??? 

  MKL_LONG test_arr();                                   //@
  MKL_LONG test_arr_symm();                              //@
  MKL_LONG nonzero();
  MKL_LONG ABS_cond(double x);
  double ABS_sum();
  MKL_LONG force_clean(double TOLERANCE);
  double sum();
  double average();
  MKL_LONG add(MATRIX& given_MAT);


  

  // Constructor
  MATRIX();                                         //@
  MATRIX(MKL_LONG N_r, MKL_LONG N_c);                         //@
  MATRIX(MKL_LONG N_r, MKL_LONG N_c, double x);               //@
  // Copy Constructor
  MKL_LONG copy_obj(const MATRIX& Mat);
  MATRIX(const MATRIX& Mat);                        //@
  /* int connect_copy_constructor(const MATRIX& Mat); */
  // Destructor
  virtual ~MATRIX();                                        //@

  MKL_LONG data_delete();                                //@
  MKL_LONG undefined_error();                            //@
  /* MKL_LONG index(MKL_LONG i, MKL_LONG j);                          //@ */ // inlined
  MKL_LONG sort();
  MKL_LONG sort2(MATRIX& given_index);
  /* MKL_LONG sort2(size_t* index); */
  // inlined functions
  // note that function inside class definition is inlined
  // these are frequently called function that have small amount of computation

  double norm()
  {
    double result = 0.0;
    for(MKL_LONG  i=0; i<size; i++)
      {
        result += pow(data[i],2.0);
      }
    return sqrt(result);
  }

  
  MKL_LONG set_value(double x)
  {
    if(INITIALIZATION)
      {
        for(MKL_LONG  i=0; i<size; i++)
          {
            data[i] = x;
          } // i
      }
    else
      {
        undefined_error();
      }
    return 0;
  }

  
  MKL_LONG index(MKL_LONG i, MKL_LONG j)
  {
    return i*cols+j;
  }
  
  // operator
  // inlined operator overloading
  MATRIX& operator=(const MATRIX &Mat)
    {
      copy_obj(Mat);
      return *this;
    }

  MATRIX& operator+=(const MATRIX &Mat)
    {
      // std::cout << "the operator += for MATRIX class is in buggy status\n";
      // std::cout << "if this operator is called, to check the internal problems\n";
      if (size != Mat.size)
        std::cout << "The += operator must happen when two size are the same\n";
      // copy_obj(Mat);
      // for(MKL_LONG i=0; i<size; i++)
      //   {
      //     data[i] += Mat.data[i];
      //   }
      cblas_daxpy(size, 1.0, Mat.data, 1, data, 1);
      return *this;
    }

  double& operator()(MKL_LONG  i, MKL_LONG  j)
  {
    return data[i*cols+j];
  }

  double& operator()(MKL_LONG  i)
  {
    return data[i];
  }

  
};


inline MKL_LONG  make_unit_vector(MATRIX& given_vec)
{
  double norm = 0.;
  for(MKL_LONG  i=0; i<given_vec.size; i++)
    {
      norm += pow(given_vec.data[i], 2.0);
    }
  norm = sqrt(norm);
  for(MKL_LONG  i=0; i<given_vec.size; i++)
    {
      given_vec.data[i] /= norm;
    }
  return norm;
}

inline MKL_LONG matrix_mul(MATRIX& given_vec, double val)
{
  cblas_dscal(given_vec.size, val, given_vec.data, 1);
  return 0;
}


MKL_LONG nonzero(const MATRIX &A);

#endif 
