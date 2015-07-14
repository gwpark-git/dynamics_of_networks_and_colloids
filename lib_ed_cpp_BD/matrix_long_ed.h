
#ifndef MATRIX_LONG_ED_H
#define MATRIX_LONG_ED_H

#define TRUE 1
#define FALSE 0
#define BIT 64
#include <iostream>
#include <fstream>
#include <mkl.h>
#include <math.h>

class MATRIX_LONG
{
 public:
  MKL_LONG rows;
  MKL_LONG cols;
  MKL_LONG size;
  MKL_LONG *data;
  bool INITIALIZATION;

  MATRIX_LONG& operator=(const MATRIX_LONG &Mat);
  MATRIX_LONG& operator+=(const MATRIX_LONG &Mat);
  MKL_LONG& operator()(MKL_LONG i, MKL_LONG j);
  MKL_LONG& operator()(MKL_LONG i);
  MATRIX_LONG ROW(MKL_LONG i);
  MATRIX_LONG COL(MKL_LONG j);
  MKL_LONG ROW(const MATRIX_LONG &ROW_A, MKL_LONG i);
  MKL_LONG COL(const MATRIX_LONG &COL_A, MKL_LONG j);
  MATRIX_LONG();
  MATRIX_LONG(MKL_LONG  N_r, MKL_LONG  N_c);
  MATRIX_LONG(MKL_LONG  N_r, MKL_LONG  N_c, MKL_LONG x);
  MKL_LONG copy_obj(const MATRIX_LONG& Mat);
  MATRIX_LONG(const MATRIX_LONG& Mat);
  ~MATRIX_LONG();
  MKL_LONG set_value(MKL_LONG x);
  MKL_LONG print();
  MKL_LONG print(MKL_LONG  n_row, MKL_LONG  n_col);
  MKL_LONG fprint_skip(const char *fn, MKL_LONG  N_skip);
  MKL_LONG fprint(const char *fn);
  MKL_LONG fprint_row(const char *fn, MKL_LONG  given_row);
  MKL_LONG fprint_out_skip(const char *fn, MKL_LONG  N_skip);
  MKL_LONG fprint_out(const char *fn);
  /* MKL_LONG nonzero(); */
  /* MKL_LONG ABS_cond(double x); */
  MKL_LONG initial();
  MKL_LONG initial(MKL_LONG  N_r, MKL_LONG  N_c);
  MKL_LONG initial(MKL_LONG  N_r, MKL_LONG  N_c, MKL_LONG x);
  MKL_LONG data_delete();
  MKL_LONG undefined_error();
  MKL_LONG index(MKL_LONG  i, MKL_LONG  j);
  MKL_LONG norm();
  MKL_LONG sum();
  double average();
};

/* MKL_LONG nonzero(const MATRIX_LONG &A); */
/* MATRIX_LONG identity(MKL_LONG RANK); */
/* MATRIX_LONG diagonal(double *vals, MKL_LONG  RANK); */
MATRIX_LONG transpose(const MATRIX_LONG &A);
MATRIX_LONG operator-(const MATRIX_LONG &A);
MATRIX_LONG operator+(const MATRIX_LONG &A, const MATRIX_LONG &B);
MATRIX_LONG operator-(const MATRIX_LONG& A, const MATRIX_LONG &B);
MATRIX_LONG operator*(const MKL_LONG a, const MATRIX_LONG &A);
MATRIX_LONG operator*(const MATRIX_LONG &A, const MATRIX_LONG &B);



#endif
