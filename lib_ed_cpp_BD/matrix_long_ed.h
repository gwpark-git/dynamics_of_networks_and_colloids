
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
  long rows;
  long cols;
  long size;
  long *data;
  bool INITIALIZATION;

  MATRIX_LONG& operator=(const MATRIX_LONG &Mat);
  /* MATRIX_LONG& operator+=(const MATRIX_LONG &Mat); */ // deleted because of bug
  long& operator()(long i, long j);
  long& operator()(long i);
  MATRIX_LONG ROW(long i);
  MATRIX_LONG COL(long j);
  long ROW(const MATRIX_LONG &ROW_A, long i);
  long COL(const MATRIX_LONG &COL_A, long j);
  MATRIX_LONG();
  MATRIX_LONG(long  N_r, long  N_c);
  MATRIX_LONG(long  N_r, long  N_c, long x);
  long copy_obj(const MATRIX_LONG& Mat);
  MATRIX_LONG(const MATRIX_LONG& Mat);
  virtual ~MATRIX_LONG();
  long set_value(long x);
  long print();
  long print(long  n_row, long  n_col);
  long fprint_skip(const char *fn, long  N_skip);
  long fprint(const char *fn);
  long fprint_row(const char *fn, long  given_row);
  long fprint_out_skip(const char *fn, long  N_skip);
  long fprint_out(const char *fn);
  /* long nonzero(); */
  /* long ABS_cond(double x); */
  long initial();
  long initial(long  N_r, long  N_c);
  long initial(long  N_r, long  N_c, long x);
  long data_delete();
  long undefined_error();
  long index(long  i, long  j);
  long norm();
  long sum();
  double average();
};

/* long nonzero(const MATRIX_LONG &A); */
/* MATRIX_LONG identity(long RANK); */
/* MATRIX_LONG diagonal(double *vals, long  RANK); */
MATRIX_LONG transpose(const MATRIX_LONG &A);
MATRIX_LONG operator-(const MATRIX_LONG &A);
MATRIX_LONG operator+(const MATRIX_LONG &A, const MATRIX_LONG &B);
MATRIX_LONG operator-(const MATRIX_LONG& A, const MATRIX_LONG &B);
MATRIX_LONG operator*(const long a, const MATRIX_LONG &A);
MATRIX_LONG operator*(const MATRIX_LONG &A, const MATRIX_LONG &B);



#endif
