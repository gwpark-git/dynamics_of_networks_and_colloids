
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
    eigen_value = (double*) mkl_malloc(rows*sizeof(double), BIT);
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
  /* MKL_LONG set_value(double x);                          //@ */ // inlined
  MKL_LONG print();                                      //@                     
  MKL_LONG print(MKL_LONG n_row, MKL_LONG n_col);                  //@
  MKL_LONG fprint_skip(const char *fn, MKL_LONG N_skip);
  MKL_LONG fprint_skip_transpose(const char *fn, MKL_LONG N_skip);
  MKL_LONG fprint_LONG_skip(const char *fn, MKL_LONG N_skip);
  MKL_LONG fprint_LONG_skip_transpose(const char *fn, MKL_LONG N_skip);
  MKL_LONG fprint_LONG_skip_transpose_LIMROWS(const char *fn, MKL_LONG N_skip, MKL_LONG N_lim_rows);
  MKL_LONG fprint(const char *fn);                             // ???
  MKL_LONG fprint_transpose(const char *fn);                             // ???
  MKL_LONG fprint_LONG(const char *fn);                             // ???
  MKL_LONG fprint_LONG_transpose(const char *fn);                             // ??? 
  MKL_LONG fprint_row(const char *fn, MKL_LONG given_row);
  MKL_LONG fprint_out_skip(const char *fn, MKL_LONG N_skip);
  MKL_LONG fprint_out(const char *fn);                             // ??? 
  MKL_LONG test_arr();                                   //@
  MKL_LONG test_arr_symm();                              //@
  MKL_LONG nonzero();
  MKL_LONG ABS_cond(double x);
  double ABS_sum();
  /* double norm(); */ // inlined
  MKL_LONG force_clean(double TOLERANCE);
  double sum();
  double average();
  MKL_LONG add(MATRIX& given_MAT);


  
  /* MKL_LONG ROW(const MATRIX &ROW_A, MKL_LONG i); */
  /* MKL_LONG COL(const MATRIX &COL_A, MKL_LONG j); */

  // Matrix Operator
  /* MATRIX ROW(MKL_LONG i);                                //@  */
  /* MATRIX COL(MKL_LONG j);                                //@ */

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
  MKL_LONG sort2(MATRIX& index);

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
  /* the following declaration will be changed to definition in order to inlining */
  /* double& operator()(MKL_LONG i, MKL_LONG j);                 //@ */
  /* double& operator()(MKL_LONG i); */
  /* MATRIX& operator=(const MATRIX &Mat); */
  /* MATRIX& operator+=(const MATRIX &Mat); */

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

MKL_LONG nonzero(const MATRIX &A);
/* MKL_LONG matrix_mul(MATRIX& given_vec, double val); */
/* MKL_LONG make_unit_vector(MATRIX& given_vec); */

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
  // for(MKL_LONG  k=0; k<given_vec.size; k++)
  //   {
  //     given_vec.data[k] *= val;
  //   }
  return 0;
}



/* /\* from here, the inline operator overloading is used *\/ */
/* /\* note that if the operator is defined inside class bracket, it automatically inlined *\/ */
/* /\* here is just following recommendation for readability *\/ */
/*  MATRIX& MATRIX::operator=(const MATRIX &Mat) */
/* { */
/*   copy_obj(Mat); */
/*   return *this; */
/* } */

/*  MATRIX& MATRIX::operator+=(const MATRIX &Mat) */
/* { */
/*   // std::cout << "the operator += for MATRIX class is in buggy status\n"; */
/*   // std::cout << "if this operator is called, to check the internal problems\n"; */
/*   if (size != Mat.size) */
/*     std::cout << "The += operator must happen when two size are the same\n"; */
/*   // copy_obj(Mat); */
/*   // for(MKL_LONG i=0; i<size; i++) */
/*   //   { */
/*   //     data[i] += Mat.data[i]; */
/*   //   } */
/*   cblas_daxpy(size, 1.0, Mat.data, 1, data, 1); */
/*   return *this; */
/* } */

/* // Operator Overloading */
/*  double& MATRIX::operator()(MKL_LONG  i, MKL_LONG  j) */
/* { */
/*   return data[i*cols+j]; */
/* } */

/*  double& MATRIX::operator()(MKL_LONG  i) */
/* { */
/*   return data[i]; */
/* } */


/* MATRIX partition(const MATRIX &A, MKL_LONG st_row, MKL_LONG end_row, MKL_LONG st_col, MKL_LONG end_col); */
/* MATRIX identity(MKL_LONG RANK);                          //@ */
/* MATRIX diagonal(double *vals, MKL_LONG RANK);            //@  */
// Basic MATRIX Operation

/* MKL_LONG get_unit_vector(MATRIX& given_vec); */
/* MATRIX return_unit_vector(MATRIX& given_vec); */


// MATRIX Addition : C = A+B
/* MATRIX operator+(const MATRIX &A, const MATRIX &B); //@ */

/* MATRIX operator-(const MATRIX &A, const MATRIX &B); //@ */
// Scalar Multiplification : C = a*A
/* MATRIX operator*(const double a, const MATRIX &A);  //@ */
// MATRIX Multiplification : C = A*B
/* MATRIX operator*(const MATRIX &A, const MATRIX &B); //@ */

// Unary operator
/* MATRIX operator-(const MATRIX &A);                        //@ */


#endif 
