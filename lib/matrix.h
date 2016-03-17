
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
MKL_LONG set_value(double x);                          //@
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
double norm();
MKL_LONG force_clean(double TOLERANCE);
double sum();
double average();
MKL_LONG add(MATRIX& given_MAT);

// Overloading
double& operator()(MKL_LONG i, MKL_LONG j);                 //@
double& operator()(MKL_LONG i);
MATRIX& operator=(const MATRIX &Mat);
MATRIX& operator+=(const MATRIX &Mat);
// Matrix Operator
/* MATRIX ROW(MKL_LONG i);                                //@  */
/* MATRIX COL(MKL_LONG j);                                //@ */
MKL_LONG ROW(const MATRIX &ROW_A, MKL_LONG i);
MKL_LONG COL(const MATRIX &COL_A, MKL_LONG j);

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
MKL_LONG index(MKL_LONG i, MKL_LONG j);                          //@
MKL_LONG sort();
MKL_LONG sort2(MATRIX& index);

};

MKL_LONG nonzero(const MATRIX &A);

/* MATRIX partition(const MATRIX &A, MKL_LONG st_row, MKL_LONG end_row, MKL_LONG st_col, MKL_LONG end_col); */
/* MATRIX identity(MKL_LONG RANK);                          //@ */
/* MATRIX diagonal(double *vals, MKL_LONG RANK);            //@  */
// Basic MATRIX Operation

MKL_LONG make_unit_vector(MATRIX& given_vec);
/* MKL_LONG get_unit_vector(MATRIX& given_vec); */
/* MATRIX return_unit_vector(MATRIX& given_vec); */

MKL_LONG matrix_mul(MATRIX& given_vec, double val);

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
