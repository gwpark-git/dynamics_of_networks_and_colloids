
#include "matrix_ed.h"
#include <iomanip>
#include <math.h>

long nonzero(const MATRIX &A)
{
  long  cnt=0;
  for(long  i=0; i<A.size; i++)
    {
      if(A.data[i] != 0.0)
	{
	  cnt++;
	}
    }
  return cnt;
}

MATRIX& MATRIX::operator=(const MATRIX &Mat)
{
  copy_obj(Mat);
  return *this;
}

MATRIX& MATRIX::operator+=(const MATRIX &Mat)
{
  // std::cout << "the operator += for MATRIX class is in buggy status\n";
  // std::cout << "if this operator is called, to check the internal problems\n";
  if (size != Mat.size)
    std::cout << "The += operator must happen when two size are the same\n";
  // copy_obj(Mat);
  // for(long i=0; i<size; i++)
  //   {
  //     data[i] += Mat.data[i];
  //   }
  cblas_daxpy(size, 1.0, Mat.data, 1, data, 1);
  return *this;
}


// #### Unary Operator
// MATRIX operator-(const MATRIX &A)
// {
//   std::cout << "Unary operator - has overhead. Avoiding this operator for performance\n";
//   double x = -1.0;
//   MATRIX C = x*A;
//   return C;
// }


// #### Binary Operator
// MATRIX Addition
// MATRIX operator+(const MATRIX &A, const MATRIX &B)
// {
//   MATRIX C;
//   if(!(A.rows == B.rows && A.cols == B.cols))
//     {
//       std::cout << "Dimension does NOT match.\n";
//     }
//   else
//     {
//       C.initial(A.rows, A.cols);
//       for(long  i=0; i<A.size; i++)
// 	{
// 	  C.data[i] = A.data[i] + B.data[i];
// 	} // i 
//     } 
//   return C;
// }

// MATRIX operator-(const MATRIX& A, const MATRIX &B)
// {
//   MATRIX C;
//   if(!(A.rows == B.rows && A.cols == B.cols))
//     {
//       std::cout << "Dimension does NOT match.\n";
//     }
//   else
//     {
//       C.initial(A.rows, A.cols);
//       for(long  i=0; i<A.size; i++)
// 	{
// 	  C.data[i] = A.data[i] - B.data[i];
// 	} // i
//     }
//   return C;
// }

// Scalar Multiplification
// MATRIX operator*(const double a, const MATRIX &A)
// {
//   MATRIX C = A;
//   for(long  i=0; i<C.size; i++)
//     {
//       C.data[i] *= a;
//     }
//   return C;
// }

// MATRIX Multiplification
// MATRIX operator*(const MATRIX &A, const MATRIX &B)
// {
//   long  C_rows = A.rows;
//   long  C_cols = B.cols;
//   long  cal_index = A.cols;
//   MATRIX C(C_rows, C_cols, 0.0);
//   if(A.cols == B.rows)
//     {
//       long  i, j, k;
//       //  // here // #pragma omp parallel for private(i,j,k) schedule(guided)
//       for(i=0; i<A.rows; i++)
// 	{
// 	  for(j=0; j<B.cols; j++)
// 	    {
// 	      for(k=0; k<A.cols; k++)
// 		{
// 		  C.data[i*C.cols+j] += A.data[i*A.cols+k]*B.data[k*B.cols+j];
// 		  //	      C(i,j) += A(i,k)*B(k,j);
// 		} // k
// 	    }// j
// 	} // i
//     }
//   else
//     {
//       std::cout << "DIMENSION NOT MATCHED DURING MATRIX MULTIPLIFICATION\n";
//     }
//   return C;
// }
// From Here. CLASS MATRIX

// Operator Overloading
double& MATRIX::operator()(long  i, long  j)
{
  return data[i*cols+j];
}

double& MATRIX::operator()(long  i)
{
  return data[i];
}


// Matrix Operator
// MATRIX MATRIX::ROW(long  i)
// {
//   MATRIX C;
//   C.initial(1, cols);
//   for(long  j=0; j<cols; j++)
//     {
//       C.data[j] = data[i*cols+j];
//     }
//   return C;
// }

// MATRIX MATRIX::COL(long  j)
// {
//   MATRIX C;
//   C.initial(rows, 1);
//   for(long  i=0; i<rows; i++)
//     {
//       C.data[i] = data[i*cols+j];
//     }
//   return C;
// }


// long  MATRIX::ROW(const MATRIX &ROW_A, long  i)
// {
//   if(INITIALIZATION)
//     {
//       for(long  j=0; j<ROW_A.cols; j++)
// 	{
// 	  data[i*cols + j] = ROW_A.data[j];
// 	}
//     }
//   else
//     {
//       initial(1,ROW_A.cols);
//       for(long   j=0; j<ROW_A.cols; j++)
// 	{
// 	  data[j] = ROW_A.data[j];
// 	}
//     }
//   return 0;
// }

// long  MATRIX::COL(const MATRIX &COL_A, long  j)
// {
//   if(INITIALIZATION)
//     {
//       for(long  i=0; i<COL_A.rows; i++)
// 	{
// 	  data[i*cols + j] = COL_A.data[i];
// 	}
//     }
//   else
//     {
//       initial(COL_A.rows, 1);
//       for(long  i=0; i<COL_A.rows; i++)
// 	{
// 	  data[i] = COL_A.data[i];
// 	}
//     }
//   return 0;
// }


// Constructor
MATRIX::MATRIX()
{
  initial();  
  //std::cout << "Undefined dimension\n";
}

MATRIX::MATRIX(long  N_r, long  N_c)
{
  initial(N_r, N_c);
}


MATRIX::MATRIX(long  N_r, long  N_c, double x) // initilization with value x
{
  initial(N_r, N_c);
  set_value(x);
}


// Copy-Constructor
long MATRIX::copy_obj(const MATRIX& Mat)
{
  initial(Mat.rows, Mat.cols);
  for(long i=0; i<size; i++)
    {
      data[i] = Mat.data[i];
    }
  return 0;
}

MATRIX::MATRIX(const MATRIX& Mat)
{
  copy_obj(Mat);
}
  //std::cout << "COPY\n";
  // initial(Mat.rows, Mat.cols);
  // for(long  i=0; i<Mat.size; i++)
  //   {
  //     data[i] = Mat.data[i];
  //   }    
// }

// int MATRIX::connect_copy_constructor(const MATRIX& Mat)
// {
//   MATRIX(Mat);
//   return 0;
// }

// Destructor
MATRIX::~MATRIX()
{
  //std::cout << "Destructor" << this << std::endl;
  // std::cout << "MATRIX destructor is called" << std::endl;
  // destructor is tested
  if (INITIALIZATION)
    {
      data_delete();
    }
  if (DIAGONALIZATION)
    {
      mkl_free(eigen_value);
      // delete[] eigen_value;
    }
}


// Public Member Function

long  MATRIX::set_value(double x)
{
  if(INITIALIZATION)
    {
      for(long  i=0; i<size; i++)
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
long  MATRIX::print_eigen()
{
  print_eigen(rows);
  return 0;
}


long  MATRIX::print_eigen(long  n_ele)
{
  if(DIAGONALIZATION)
    {
      std::cout << "Latent Roots = ";
      for(long  i=0; i<n_ele; i++)
	{
	  std::cout << eigen_value[i] << '\t';
	}
      std::cout << std::endl;
    }
  else
    {
      std::cout << "MATRIX DOES NOT CALCULATED LATENT ROOTS\n";
      return -1;
    }
  return 0;
}

long  MATRIX::print()
{
  print(rows, cols);
  return 0;
}

long  MATRIX::print(long  n_row, long  n_col)
{
  long  cnt = 0;
  //std::cout << "MATRIX " << this << std::endl;
  for(long  i=0; i<n_row; i++)
    {
      for(long  j=0; j<n_col; j++)
	{
	  printf("%8.5lf\t", data[i*cols+j]);
	}
      printf("\n");
    } // i
  printf("\n");
  return 0;
}

long MATRIX::fprint_skip(const char *fn, long  N_skip)
{
  std::ofstream FILE1;
  FILE1.open(fn, std::ios_base::app);
  if(!FILE1)
    {
      std::cout << "ERROR OCCURS DURING FILE OUT\n";
      return -1;
    }
  long  cnt = 0;
  for(long  i=0; i<rows; i+=N_skip)
    {
      for(long  j=0; j<cols; j++)
        {
          FILE1 << std::scientific << data[index(j, i)] << '\t';
        }
      FILE1 << std::endl;
    } // i
  FILE1.close();
  return 0;
}

long  MATRIX::fprint(const char *fn)
{
  fprint_skip(fn, 1);
  return 0;
}

long MATRIX::fprint_skip_transpose(const char *fn, long  N_skip)
{
  std::ofstream FILE1;
  FILE1.open(fn, std::ios_base::app);
  if(!FILE1)
    {
      std::cout << "ERROR OCCURS DURING FILE OUT\n";
      return -1;
    }
  long  cnt = 0;
  for(long  j=0; j<cols; j++)
    {
      for(long  i=0; i<rows; i+=N_skip)
        {
          FILE1 << std::scientific << data[index(i, j)] << '\t';
        }
      FILE1 << std::endl;
    } // i
  FILE1.close();
  return 0;
}

long  MATRIX::fprint_transpose(const char *fn)
{
  fprint_skip_transpose(fn, 1);
  return 0;
}


long  MATRIX::fprint_LONG_skip(const char *fn, long  N_skip)
{
  std::ofstream FILE1;
  FILE1.open(fn, std::ios_base::app);
  if(!FILE1)
    {
      std::cout << "ERROR OCCURS DURING FILE OUT\n";
      return -1;
    }
  long  cnt = 0;
  for(long  i=0; i<rows; i+=N_skip)
    {
      for(long  j=0; j<cols; j++)
        {
          FILE1 << (long)data[index(i, j)] << '\t';
        }
      FILE1 << std::endl;
    } // i
  FILE1.close();
  return 0;
}

long  MATRIX::fprint_LONG(const char *fn)
{
  fprint_LONG_skip(fn, 1);
  return 0;
}

long  MATRIX::fprint_LONG_skip_transpose(const char *fn, long  N_skip)
{
  std::ofstream FILE1;
  FILE1.open(fn, std::ios_base::app);
  if(!FILE1)
    {
      std::cout << "ERROR OCCURS DURING FILE OUT\n";
      return -1;
    }
  long  cnt = 0;
  for(long  j=0; j<cols; j++)
    {
      for(long  i=0; i<rows; i+=N_skip)
        {
          FILE1 << (long)data[index(i, j)] << '\t';
        }
      FILE1 << std::endl;
    } // i
  FILE1.close();
  return 0;
}

long  MATRIX::fprint_LONG_transpose(const char *fn)
{
  fprint_LONG_skip_transpose(fn, 1);
  return 0;
}


long  MATRIX::fprint_row(const char *fn, long  given_row)
{
  std::ofstream FILE1;
  FILE1.open(fn, std::ios_base::app);
  if(!FILE1)
    {
      std::cout << "ERROR OCCURS DURING FILE OUT\n";
      return -1;
    }
  long  cnt = 0;
  long  i=given_row;
  for(long  j=0; j<cols; j++)
    {
      FILE1 << std::scientific << data[index(i, j)] << '\t';
    }
  FILE1 << std::endl;
  FILE1.close();
  return 0;
}

long  MATRIX::fprint_out_skip(const char *fn, long  N_skip)
{
  std::ofstream FILE1;
  FILE1.open(fn, std::ios_base::out);
  if(!FILE1)
    {
      std::cout << "ERROR OCCURS DURING FILE OUT\n";
      return -1;
    }
  long  cnt = 0;
  for(long  i=0; i<rows; i+=N_skip)
    {
      for(long  j=0; j<cols; j++)
	{
	  FILE1 << std::scientific << data[cnt++] << '\t';
	}
      FILE1 << std::endl;
    } // i
  FILE1 << std::endl;
  FILE1.close();
  return 0;

}

long  MATRIX::fprint_out(const char *fn)
{
  fprint_out_skip(fn, 1);
  return 0;
}



long  MATRIX::p_self()
{
  std::cout << this;
  return 0;
}

long  MATRIX::nonzero()
{
  for(long  i=0; i<size; i++)
    {
      if(data[i]!=0)
	{
	  printf("nonzero (i,j) = (%ld, %ld) = %6.3e\n",(int)i/cols,(int)i%cols, data[i]);
	}
    }
  return 0;
}

long  MATRIX::ABS_cond(double x)
{
  for(long  i=0; i<size; i++)
    {
      if(fabs(data[i]) > fabs(x))
	{
	  printf("nonzero (i,j) = (%ld, %ld) = %6.3e\n",(int)i/cols,(int)i%cols, data[i]);
	}
    }
  return 0;
}


// Private Member Function
long MATRIX::initial()
{
  INITIALIZATION = FALSE;
  DIAGONALIZATION = FALSE;
  rows = 0;
  cols = 0;
  size = 0;
  data = NULL;
  return 0;
}

long  MATRIX::initial(long  N_r, long  N_c)
{
  // std::cout << "CONSTRUCTOR " << this << std::endl;
  if (size == N_r*N_c && INITIALIZATION == TRUE)
    {
      rows = N_r;
      cols = N_c;
      return 0;
    }
  else
    {
      // this is of importance to prevent memory leak
      if (INITIALIZATION == TRUE)
        data_delete();
      INITIALIZATION = TRUE;
      DIAGONALIZATION = FALSE;
      rows = N_r;
      cols = N_c;
      size = rows*cols;
      // note that the dynamic allocation should be defined callocation
      // otherwise, it should be set with the valid initial variable
      data = (double*) mkl_calloc(size, sizeof(double), BIT);
      // data = new double [size];
    }
  return 0;
}

long  MATRIX::initial(long  N_r, long  N_c, double x)
{
  initial(N_r, N_c);
  for(long  i=0; i<size; i++)
    {
      data[i] = x;
    }
  return 0;
}


long MATRIX::data_delete()
{
  // delete[] data;
  mkl_free(data);
  INITIALIZATION = FALSE;
  data = NULL;
  return 0;
}

long  MATRIX::undefined_error()
{
  std::cout << "MATRIX is NOT constructed\n" << std::endl;
  return -1;
}


long  MATRIX::index(long  i, long  j)
{
  return i*cols+j;
}


long  MATRIX::test_arr()
{
  for(long  i=0; i<size; i++)
    {
      data[i] = i+1;
    } // i
  return 0;
}

long  MATRIX::test_arr_symm()
{
  for(long  i=0; i<rows; i++)
    {
      for(long  j=0; j<cols; j++)
        {
          data[i*cols+j] = 1.0/(1.0+i+j);
        }
    }
  return 0;
}

double MATRIX::ABS_sum()
{
  double result = 0.0;
  double tmp_x;
  for(long  i=0; i<size; i++)
    {
      tmp_x = data[i];
      if(tmp_x < 0)
	result -= tmp_x;
      else
	result += tmp_x;
    }
  return result;
}

double MATRIX::norm()
{
  double result = 0.0;
  for(long  i=0; i<size; i++)
    {
      result += pow(data[i],2.0);
    }
  return sqrt(result);
}



double MATRIX::sum()
{
  double result = 0.0;
  for(long  i=0; i<size; i++)
    {
      result += data[i];
    }
  return result;
}

double MATRIX::average()
{
  double result = sum();
  return result/(double)size;
}

long  make_unit_vector(MATRIX& given_vec)
{
  double norm = 0.;
  for(long  i=0; i<given_vec.size; i++)
    {
      norm += pow(given_vec.data[i], 2.0);
    }
  norm = sqrt(norm);
  for(long  i=0; i<given_vec.size; i++)
    {
      given_vec.data[i] /= norm;
    }
  return norm;
}

// MATRIX return_unit_vector(MATRIX& given_vec)
// {
//   MATRIX re = given_vec;
//   double norm = 0.;
//   for(long  i=0; i<given_vec.size; i++)
//     {
//       norm += pow(given_vec.data[i], 2.0);
//     }
//   norm = sqrt(norm);
//   for(long  i=0; i<given_vec.size; i++)
//     {
//       re.data[i] /= norm;
//     }
//   return re;
// }



long matrix_mul(MATRIX& given_vec, double val)
{
  cblas_dscal(given_vec.size, val, given_vec.data, 1);
  // for(long  k=0; k<given_vec.size; k++)
  //   {
  //     given_vec.data[k] *= val;
  //   }
  return 0;
}


long MATRIX::add(MATRIX& given_MAT)
{
  if (size != given_MAT.size)
    {
      std::cout << "ERR: MATRIX::add different size\n";
      return -1;
    }
  cblas_daxpy(size, 1.0, given_MAT.data, 1, data, 1);
  // for(long  k=0; k<size; k++)
  //   {
  //     data[k] += given_MAT.data[k];
  //   }
  return 0;
}

long MATRIX::sort()
{
  gsl_sort(data, 1, size);
  return 0;
}

long MATRIX::sort2(MATRIX& index)
{
  gsl_sort2(data, 1, index.data, 1, size);
  return 0;
}
