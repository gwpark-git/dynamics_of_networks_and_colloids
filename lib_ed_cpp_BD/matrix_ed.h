
/*
#################################################################
##   Private MATRIX Library ver 0.3                            ##
##   Made by "Gun Woo Park", Email : gunwoo.park@unina.it      ##
##   Fellow of SUPOLEN, Marie Curie Actions                    ##
##   DICMaPI, University of Naples Federico II                 ##
#################################################################
*/

// all the integer is converted to long

// @ denoted for tested function
// Lots of functions are based on Call by Reference.
// Notice that for matrix operation, we have to put &.
// Vector can use by rows = 1.

// Note that I seperated library for linear algebra and matrix.
// If some error occurs during compiling, please check to include linear algebra code.

#ifndef MATRIX_ED_H
#define MATRIX_ED_H
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
long rows;
long cols;
long size;
double *data;  // 1-dimensional array
bool VECTOR; 
bool INITIALIZATION;
bool DIAGONALIZATION;
double *eigen_value;
long DIAGONALIZATION_INITIAL()                     //@ 
{
/* eigen_value = new double [rows]; */
  eigen_value = (double*) mkl_malloc(rows*sizeof(double), BIT);
  DIAGONALIZATION = TRUE;
return 0;
}
long print_eigen();
long print_eigen(long n_ele);                                //@
// Public Member Function
long p_self();                                     //@
long initial();
long initial(long N_r, long N_c);                    //@
long initial(long N_r, long N_c, double x);
long set_value(double x);                          //@
long print();                                      //@                     
long print(long n_row, long n_col);                  //@
long fprint_skip(const char *fn, long N_skip);
long fprint_skip_transpose(const char *fn, long N_skip);
long fprint_LONG_skip(const char *fn, long N_skip);
long fprint_LONG_skip_transpose(const char *fn, long N_skip);
long fprint(const char *fn);                             // ???
long fprint_transpose(const char *fn);                             // ???
long fprint_LONG(const char *fn);                             // ???
long fprint_LONG_transpose(const char *fn);                             // ??? 
long fprint_row(const char *fn, long given_row);
long fprint_out_skip(const char *fn, long N_skip);
long fprint_out(const char *fn);                             // ??? 
long test_arr();                                   //@
long test_arr_symm();                              //@
long nonzero();
long ABS_cond(double x);
double ABS_sum();
double norm();
long force_clean(double TOLERANCE);
double sum();
double average();
long add(MATRIX& given_MAT);

// Overloading
double& operator()(long i, long j);                 //@
double& operator()(long i);
MATRIX& operator=(const MATRIX &Mat);
MATRIX& operator+=(const MATRIX &Mat);
// Matrix Operator
/* MATRIX ROW(long i);                                //@  */
/* MATRIX COL(long j);                                //@ */
long ROW(const MATRIX &ROW_A, long i);
long COL(const MATRIX &COL_A, long j);

// Constructor
MATRIX();                                         //@
MATRIX(long N_r, long N_c);                         //@
MATRIX(long N_r, long N_c, double x);               //@
// Copy Constructor
long copy_obj(const MATRIX& Mat);
MATRIX(const MATRIX& Mat);                        //@
/* int connect_copy_constructor(const MATRIX& Mat); */
// Destructor
virtual ~MATRIX();                                        //@

long data_delete();                                //@
long undefined_error();                            //@
long index(long i, long j);                          //@
long sort();
long sort2(MATRIX& index);

};

long nonzero(const MATRIX &A);

/* MATRIX partition(const MATRIX &A, long st_row, long end_row, long st_col, long end_col); */
/* MATRIX identity(long RANK);                          //@ */
/* MATRIX diagonal(double *vals, long RANK);            //@  */
// Basic MATRIX Operation

long make_unit_vector(MATRIX& given_vec);
/* long get_unit_vector(MATRIX& given_vec); */
/* MATRIX return_unit_vector(MATRIX& given_vec); */

long matrix_mul(MATRIX& given_vec, double val);

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
