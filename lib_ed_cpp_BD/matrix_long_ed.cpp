
#include "matrix_long_ed.h"

// MKL_LONG nonzero(const MATRIX_LONG &A)
// {
//   MKL_LONG  cnt=0;
//   for(MKL_LONG  i=0; i<A.size; i++)
//     {
//       if(A.data[i] != 0)
//         {
//           cnt++;
//         }
//     }
//   return cnt;
// }


// MATRIX_LONG identity(MKL_LONG RANK)
// {
//   MATRIX_LONG C(RANK, RANK, 0);
//   for(MKL_LONG  i=0; i<RANK; i++)
//     {
//       C(i,i) = 1.0;
//     }
//   return C;
// }

// MATRIX_LONG diagonal(double *vals, MKL_LONG  RANK)
// {
//   MATRIX_LONG C(RANK, RANK, 0);
//   for(MKL_LONG  i=0; i<RANK; i++)
//     {
//       C(i,i) = vals[i];
//     }
//   return C;
// }

MATRIX_LONG transpose(const MATRIX_LONG &A)
{
  MATRIX_LONG C(A.cols, A.rows);
  for(MKL_LONG  i=0; i<A.rows; i++)
    {
      for(MKL_LONG  j=0; j < A.cols; j++)
        {
          C.data[j*A.rows+i] = A.data[i*A.cols+j];
        } // j
    } // i
  return C;
}


MATRIX_LONG& MATRIX_LONG::operator=(const MATRIX_LONG &Mat)
{
  copy_obj(Mat);
  return *this;
}

MATRIX_LONG& MATRIX_LONG::operator+=(const MATRIX_LONG &Mat)
{
  copy_obj(Mat);
  return *this;
}


// #### Unary Operator
MATRIX_LONG operator-(const MATRIX_LONG &A)
{
  double x = -1.0;
  MATRIX_LONG C = x*A;
  return C;
}


// #### Binary Operator
// MATRIX_LONG Addition
MATRIX_LONG operator+(const MATRIX_LONG &A, const MATRIX_LONG &B)
{
  MATRIX_LONG C;
  if(!(A.rows == B.rows && A.cols == B.cols))
    {
      std::cout << "Dimension does NOT match.\n";
    }
  else
    {
      C.initial(A.rows, A.cols);
      for(MKL_LONG  i=0; i<A.size; i++)
        {
          C.data[i] = A.data[i] + B.data[i];
        } // i 
    } 
  return C;
}

MATRIX_LONG operator-(const MATRIX_LONG& A, const MATRIX_LONG &B)
{
  MATRIX_LONG C;
  if(!(A.rows == B.rows && A.cols == B.cols))
    {
      std::cout << "Dimension does NOT match.\n";
    }
  else
    {
      C.initial(A.rows, A.cols);
      for(MKL_LONG  i=0; i<A.size; i++)
        {
          C.data[i] = A.data[i] - B.data[i];
        } // i
    }
  return C;
}

// Scalar Multiplification
MATRIX_LONG operator*(const MKL_LONG a, const MATRIX_LONG &A)
{
  MATRIX_LONG C = A;
  for(MKL_LONG  i=0; i<C.size; i++)
    {
      C.data[i] *= a;
    }
  return C;
}

// MATRIX_LONG Multiplification
MATRIX_LONG operator*(const MATRIX_LONG &A, const MATRIX_LONG &B)
{
  MKL_LONG  C_rows = A.rows;
  MKL_LONG  C_cols = B.cols;
  MKL_LONG  cal_index = A.cols;
  MATRIX_LONG C(C_rows, C_cols, 0);
  if(A.cols == B.rows)
    {
      MKL_LONG  i, j, k;
      //  // here // #pragma omp parallel for private(i,j,k) schedule(guided)
      for(i=0; i<A.rows; i++)
        {
          for(j=0; j<B.cols; j++)
            {
              for(k=0; k<A.cols; k++)
                {
                  C.data[i*C.cols+j] += A.data[i*A.cols+k]*B.data[k*B.cols+j];
                  //	      C(i,j) += A(i,k)*B(k,j);
                } // k
            }// j
        } // i
    }
  else
    {
      std::cout << "DIMENSION NOT MATCHED DURING MATRIX_LONG MULTIPLIFICATION\n";
    }
  return C;
}
// From Here. CLASS MATRIX_LONG

// Operator Overloading
MKL_LONG& MATRIX_LONG::operator()(MKL_LONG i, MKL_LONG j)
{
  return data[i*cols+j];
}

MKL_LONG& MATRIX_LONG::operator()(MKL_LONG i)
{
  return data[i];
}


// Matrix Operator
MATRIX_LONG MATRIX_LONG::ROW(MKL_LONG i)
{
  MATRIX_LONG C;
  C.initial(1, cols);
  for(MKL_LONG  j=0; j<cols; j++)
    {
      C.data[j] = data[i*cols+j];
    }
  return C;
}

MATRIX_LONG MATRIX_LONG::COL(MKL_LONG j)
{
  MATRIX_LONG C;
  C.initial(rows, 1);
  for(MKL_LONG  i=0; i<rows; i++)
    {
      C.data[i] = data[i*cols+j];
    }
  return C;
}


MKL_LONG MATRIX_LONG::ROW(const MATRIX_LONG &ROW_A, MKL_LONG i)
{
  if(INITIALIZATION)
    {
      for(MKL_LONG j=0; j<ROW_A.cols; j++)
        {
          data[i*cols + j] = ROW_A.data[j];
        }
    }
  else
    {
      initial(1,ROW_A.cols);
      for(MKL_LONG j=0; j<ROW_A.cols; j++)
        {
          data[j] = ROW_A.data[j];
        }
    }
  return 0;
}

MKL_LONG MATRIX_LONG::COL(const MATRIX_LONG &COL_A, MKL_LONG j)
{
  if(INITIALIZATION)
    {
      for(MKL_LONG  i=0; i<COL_A.rows; i++)
        {
          data[i*cols + j] = COL_A.data[i];
        }
    }
  else
    {
      initial(COL_A.rows, 1);
      for(MKL_LONG  i=0; i<COL_A.rows; i++)
        {
          data[i] = COL_A.data[i];
        }
    }
  return 0;
}


// Constructor
MATRIX_LONG::MATRIX_LONG()
{
  initial();  
  //std::cout << "Undefined dimension\n";
}

MATRIX_LONG::MATRIX_LONG(MKL_LONG  N_r, MKL_LONG  N_c)
{
  initial(N_r, N_c);
}


MATRIX_LONG::MATRIX_LONG(MKL_LONG  N_r, MKL_LONG  N_c, MKL_LONG x) // initilization with value x
{
  initial(N_r, N_c);
  set_value(x);
}

// Copy-Constructor
MKL_LONG MATRIX_LONG::copy_obj(const MATRIX_LONG& Mat)
{
  initial(Mat.rows, Mat.cols);
  for(MKL_LONG i=0; i<size; i++)
    {
      data[i] = Mat.data[i];
    }
  return 0;
}

MATRIX_LONG::MATRIX_LONG(const MATRIX_LONG& Mat)
{
  copy_obj(Mat);
}
//std::cout << "COPY\n";
// initial(Mat.rows, Mat.cols);
// for(MKL_LONG  i=0; i<Mat.size; i++)
//   {
//     data[i] = Mat.data[i];
//   }    
// }

// int MATRIX_LONG::connect_copy_constructor(const MATRIX_LONG& Mat)
// {
//   MATRIX_LONG(Mat);
//   return 0;
// }

// Destructor
MATRIX_LONG::~MATRIX_LONG()
{
  //std::cout << "Destructor" << this << std::endl;
  // std::cout << "MATRIX_LONG destructor is called" << std::endl;
  // destructor is tested
  if (INITIALIZATION)
    {
      data_delete();
      data = NULL;
    }
}


// Public Member Function

MKL_LONG MATRIX_LONG::set_value(MKL_LONG x)
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

MKL_LONG  MATRIX_LONG::print()
{
  print(rows, cols);
  return 0;
}

MKL_LONG  MATRIX_LONG::print(MKL_LONG  n_row, MKL_LONG  n_col)
{
  MKL_LONG  cnt = 0;
  //std::cout << "MATRIX_LONG " << this << std::endl;
  for(MKL_LONG  i=0; i<n_row; i++)
    {
      for(MKL_LONG  j=0; j<n_col; j++)
        {
          printf("%16ld\t", data[i*cols+j]);
        }
      printf("\n");
    } // i
  printf("\n");
  return 0;
}

MKL_LONG  MATRIX_LONG::fprint_skip(const char *fn, MKL_LONG  N_skip)
{
  std::ofstream FILE1;
  FILE1.open(fn, std::ios_base::app);
  if(!FILE1)
    {
      std::cout << "ERROR OCCURS DURING FILE OUT\n";
      return -1;
    }
  MKL_LONG  cnt = 0;
  for(MKL_LONG  i=0; i<rows; i+=N_skip)
    {
      for(MKL_LONG  j=0; j<cols; j++)
        {
          FILE1 << data[index(i, j)] << '\t';
        }
      FILE1 << std::endl;
    } // i
  FILE1.close();
  return 0;
}

MKL_LONG  MATRIX_LONG::fprint(const char *fn)
{
  fprint_skip(fn, 1);
  return 0;
}

MKL_LONG  MATRIX_LONG::fprint_row(const char *fn, MKL_LONG  given_row)
{
  std::ofstream FILE1;
  FILE1.open(fn, std::ios_base::app);
  if(!FILE1)
    {
      std::cout << "ERROR OCCURS DURING FILE OUT\n";
      return -1;
    }
  MKL_LONG  cnt = 0;
  MKL_LONG  i=given_row;
  for(MKL_LONG  j=0; j<cols; j++)
    {
      FILE1 << data[index(i, j)] << '\t';
    }
  FILE1 << std::endl;
  FILE1.close();
  return 0;
}

MKL_LONG  MATRIX_LONG::fprint_out_skip(const char *fn, MKL_LONG  N_skip)
{
  std::ofstream FILE1;
  FILE1.open(fn, std::ios_base::out);
  if(!FILE1)
    {
      std::cout << "ERROR OCCURS DURING FILE OUT\n";
      return -1;
    }
  MKL_LONG  cnt = 0;
  for(MKL_LONG  i=0; i<rows; i+=N_skip)
    {
      for(MKL_LONG  j=0; j<cols; j++)
        {
          FILE1 << data[cnt++] << '\t';
        }
      FILE1 << std::endl;
    } // i
  FILE1 << std::endl;
  FILE1.close();
  return 0;

}

MKL_LONG MATRIX_LONG::fprint_out(const char *fn)
{
  fprint_out_skip(fn, 1);
  return 0;
}


// MKL_LONG MATRIX_LONG::nonzero()
// {
//   for(MKL_LONG  i=0; i<size; i++)
//     {
//       if(data[i]!=0)
//         {
//           printf("nonzero (i,j) = (%ld, %ld) = %ld\n",(int)i/cols,(int)i%cols, data[i]);
//         }
//     }
//   return 0;
// }

// MKL_LONG MATRIX_LONG::ABS_cond(double x)
// {
//   for(MKL_LONG  i=0; i<size; i++)
//     {
//       if(fabs(data[i]) > fabs(x))
//         {
//           printf("nonzero (i,j) = (%ld, %ld) = %6.3e\n",(int)i/cols,(int)i%cols, data[i]);
//         }
//     }
//   return 0;
// }


// Private Member Function
MKL_LONG MATRIX_LONG::initial()
{
  INITIALIZATION = FALSE;
  rows = 0;
  cols = 0;
  size = 0;
  data = NULL;
  return 0;
}

MKL_LONG MATRIX_LONG::initial(MKL_LONG  N_r, MKL_LONG  N_c)
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
      if (INITIALIZATION == TRUE)
        data_delete();
      INITIALIZATION = TRUE;
      rows = N_r;
      cols = N_c;
      size = rows*cols;
      data = new MKL_LONG [size];
    }
  return 0;
}

MKL_LONG MATRIX_LONG::initial(MKL_LONG  N_r, MKL_LONG  N_c, MKL_LONG x)
{
  initial(N_r, N_c);
  for(MKL_LONG  i=0; i<size; i++)
    {
      data[i] = x;
    }
  return 0;
}


MKL_LONG MATRIX_LONG::data_delete()
{
  delete[] data;
  INITIALIZATION = FALSE;
  return 0;
}

MKL_LONG MATRIX_LONG::undefined_error()
{
  std::cout << "MATRIX_LONG is NOT constructed\n" << std::endl;
  return -1;
}


MKL_LONG  MATRIX_LONG::index(MKL_LONG  i, MKL_LONG  j)
{
  return i*cols+j;
}


MKL_LONG MATRIX_LONG::norm()
{
  MKL_LONG result = 0;
  for(MKL_LONG  i=0; i<size; i++)
    {
      result += pow(data[i],2);
    }
  return sqrt(result);
}



MKL_LONG MATRIX_LONG::sum()
{
  MKL_LONG result = 0;
  for(MKL_LONG  i=0; i<size; i++)
    {
      result += data[i];
    }
  return result;
}

double MATRIX_LONG::average()
{
  double result = sum();
  return result/(double)size;
}


// MKL_LONG matrix_mul(MATRIX_LONG& given_vec, double val)
// {
//   for(MKL_LONG  k=0; k<given_vec.size; k++)
//     {
//       given_vec.data[k] *= val;
//     }
//   return 0;
// }


// MKL_LONG MATRIX_LONG::add(MATRIX_LONG& given_MAT)
// {
//   if (size == given_MAT.size)
//     {
//       for(MKL_LONG  k=0; k<size; k++)
//         {
//           data[k] += given_MAT.data[k];
//         }
//     }
//   return 0;
// }

