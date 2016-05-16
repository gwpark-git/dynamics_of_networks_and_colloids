
#include "matrix.h"
#include <iomanip>
#include <math.h>

MKL_LONG nonzero(const MATRIX &A)
{
  MKL_LONG  cnt=0;
  for(MKL_LONG  i=0; i<A.size; i++)
    {
      if(A.data[i] != 0.0)
	{
	  cnt++;
	}
    }
  return cnt;
}






// Constructor
MATRIX::MATRIX()
{
  initial();  
  //std::cout << "Undefined dimension\n";
}

MATRIX::MATRIX(MKL_LONG  N_r, MKL_LONG  N_c)
{
  initial(N_r, N_c);
}


MATRIX::MATRIX(MKL_LONG  N_r, MKL_LONG  N_c, double x) // initilization with value x
{
  initial(N_r, N_c);
  set_value(x);
}


// Copy-Constructor
MKL_LONG MATRIX::copy_obj(const MATRIX& Mat)
{
  initial(Mat.rows, Mat.cols);
  for(MKL_LONG i=0; i<size; i++)
    {
      data[i] = Mat.data[i];
    }
  return 0;
}

MATRIX::MATRIX(const MATRIX& Mat)
{
  copy_obj(Mat);
}

// Destructor
MATRIX::~MATRIX()
{
  if (INITIALIZATION)
    {
      data_delete();
    }
  if (DIAGONALIZATION)
    {
      mkl_free(eigen_value);
    }
}


// Public Member Function

// inlined
// MKL_LONG  MATRIX::set_value(double x)
// {
//   if(INITIALIZATION)
//     {
//       for(MKL_LONG  i=0; i<size; i++)
//         {
//           data[i] = x;
//         } // i
//     }
//   else
//     {
//       undefined_error();
//     }
//   return 0;
// }


MKL_LONG  MATRIX::print_eigen()
{
  print_eigen(rows);
  return 0;
}


MKL_LONG  MATRIX::print_eigen(MKL_LONG  n_ele)
{
  if(DIAGONALIZATION)
    {
      std::cout << "Latent Roots = ";
      for(MKL_LONG  i=0; i<n_ele; i++)
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

MKL_LONG  MATRIX::print()
{
  print(rows, cols);
  return 0;
}

MKL_LONG  MATRIX::print(MKL_LONG  n_row, MKL_LONG  n_col)
{
  MKL_LONG  cnt = 0;
  //std::cout << "MATRIX " << this << std::endl;
  for(MKL_LONG  i=0; i<n_row; i++)
    {
      for(MKL_LONG  j=0; j<n_col; j++)
	{
	  printf("%8.5lf\t", data[i*cols+j]);
	}
      printf("\n");
    } // i
  printf("\n");
  return 0;
}

MKL_LONG MATRIX::fprint_skip(std::ofstream& file, MKL_LONG  N_skip)
{
  // std::ofstream FILE1;
  // FILE1.open(fn, std::ios_base::app);
  // if(!FILE1)
  //   {
  //     std::cout << "ERROR OCCURS DURING FILE OUT\n";
  //     return -1;
  //   }
  MKL_LONG  cnt = 0;
  for(MKL_LONG  i=0; i<rows; i+=N_skip)
    {
      for(MKL_LONG  j=0; j<cols; j++)
        {
          file << std::scientific << data[index(j, i)] << '\t';
        }
      file << std::endl;
    } // i
  // FILE1.close();
  return 0;
}

MKL_LONG  MATRIX::fprint(std::ofstream& file)
{
  fprint_skip(file, 1);
  return 0;
}

MKL_LONG MATRIX::fprint_skip_transpose(std::ofstream& file, MKL_LONG  N_skip)
{
  // std::ofstream FILE1;
  // FILE1.open(fn, std::ios_base::app);
  // if(!FILE1)
  //   {
  //     std::cout << "ERROR OCCURS DURING FILE OUT\n";
  //     return -1;
  //   }
  MKL_LONG  cnt = 0;
  for(MKL_LONG  j=0; j<cols; j++)
    {
      for(MKL_LONG  i=0; i<rows; i+=N_skip)
        {
          file << std::scientific << data[index(i, j)] << '\t';
        }
      file << std::endl;
    } // i
  // FILE1.close();
  return 0;
}

MKL_LONG  MATRIX::fprint_transpose(std::ofstream& file)
{
  fprint_skip_transpose(file, 1);
  return 0;
}


MKL_LONG  MATRIX::fprint_LONG_skip(std::ofstream& file, MKL_LONG  N_skip)
{
  // std::ofstream FILE1;
  // FILE1.open(fn, std::ios_base::app);
  // if(!FILE1)
  //   {
  //     std::cout << "ERROR OCCURS DURING FILE OUT\n";
  //     return -1;
  //   }
  MKL_LONG  cnt = 0;
  for(MKL_LONG  i=0; i<rows; i+=N_skip)
    {
      for(MKL_LONG  j=0; j<cols; j++)
        {
          file << (MKL_LONG)data[index(i, j)] << '\t';
        }
      file << std::endl;
    } // i
  // FILE1.close();
  return 0;
}

MKL_LONG  MATRIX::fprint_LONG(std::ofstream& file)
{
  fprint_LONG_skip(file, 1);
  return 0;
}

MKL_LONG MATRIX::fprint_LONG_skip_transpose_LIMROWS(std::ofstream& file, MKL_LONG N_skip, MKL_LONG N_lim_rows)
{
  // std::ofstream FILE1;
  // FILE1.open(fn, std::ios_base::app);
  // if(!FILE1)
  //   {
  //     std::cout << "ERROR OCCURS DURING FILE OUT\n";
  //     return -1;
  //   }
  MKL_LONG  cnt = 0;
  for(MKL_LONG  j=0; j<cols; j++)
    {
      for(MKL_LONG  i=0; i<N_lim_rows; i+=N_skip)
        {
          file << (MKL_LONG)data[index(i, j)] << '\t';
        }
      file << std::endl;
    } // i
  // FILE1.close();
  return 0;
}

MKL_LONG MATRIX::fprint_LONG_skip_transpose(std::ofstream& file, MKL_LONG  N_skip)
{
  MATRIX::fprint_LONG_skip_transpose_LIMROWS(file, N_skip, rows); // count all rows
  return 0;
}

MKL_LONG  MATRIX::fprint_LONG_transpose(std::ofstream& file)
{
  fprint_LONG_skip_transpose(file, 1);
  return 0;
}


MKL_LONG  MATRIX::fprint_row(std::ofstream& file, MKL_LONG  given_row)
{
  // std::ofstream FILE1;
  // FILE1.open(fn, std::ios_base::app);
  // if(!FILE1)
  //   {
  //     std::cout << "ERROR OCCURS DURING FILE OUT\n";
  //     return -1;
  //   }
  MKL_LONG  cnt = 0;
  MKL_LONG  i=given_row;
  for(MKL_LONG  j=0; j<cols; j++)
    {
      file << std::scientific << data[index(i, j)] << '\t';
    }
  file << std::endl;
  // FILE1.close();
  return 0;
}

MKL_LONG  MATRIX::fprint_out_skip(std::ofstream& file, MKL_LONG  N_skip)
{
  // std::ofstream FILE1;
  // FILE1.open(fn, std::ios_base::out);
  // if(!FILE1)
  //   {
  //     std::cout << "ERROR OCCURS DURING FILE OUT\n";
  //     return -1;
  //   }
  MKL_LONG  cnt = 0;
  for(MKL_LONG  i=0; i<rows; i+=N_skip)
    {
      for(MKL_LONG  j=0; j<cols; j++)
        {
          file << std::scientific << data[cnt++] << '\t';
        }
      file << std::endl;
    } // i
  file << std::endl;
  // FILE1.close();
  return 0;

}

MKL_LONG  MATRIX::fprint_out(std::ofstream& file)
{
  fprint_out_skip(file, 1);
  return 0;
}



MKL_LONG  MATRIX::p_self()
{
  std::cout << this;
  return 0;
}

MKL_LONG  MATRIX::nonzero()
{
  for(MKL_LONG  i=0; i<size; i++)
    {
      if(data[i]!=0)
        {
          printf("nonzero (i,j) = (%ld, %ld) = %6.3e\n",(int)i/cols,(int)i%cols, data[i]);
        }
    }
  return 0;
}

MKL_LONG  MATRIX::ABS_cond(double x)
{
  for(MKL_LONG  i=0; i<size; i++)
    {
      if(fabs(data[i]) > fabs(x))
	{
	  printf("nonzero (i,j) = (%ld, %ld) = %6.3e\n",(int)i/cols,(int)i%cols, data[i]);
	}
    }
  return 0;
}


// Private Member Function
MKL_LONG MATRIX::initial()
{
  INITIALIZATION = FALSE;
  DIAGONALIZATION = FALSE;
  rows = 0;
  cols = 0;
  size = 0;
  data = NULL;
  return 0;
}

MKL_LONG  MATRIX::initial(MKL_LONG  N_r, MKL_LONG  N_c)
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

MKL_LONG  MATRIX::initial(MKL_LONG  N_r, MKL_LONG  N_c, double x)
{
  initial(N_r, N_c);
  for(MKL_LONG  i=0; i<size; i++)
    {
      data[i] = x;
    }
  return 0;
}


MKL_LONG MATRIX::data_delete()
{
  // delete[] data;
  mkl_free(data);
  INITIALIZATION = FALSE;
  data = NULL;
  return 0;
}

MKL_LONG  MATRIX::undefined_error()
{
  std::cout << "MATRIX is NOT constructed\n" << std::endl;
  return -1;
}


// inlined
// MKL_LONG  MATRIX::index(MKL_LONG  i, MKL_LONG  j)
// {
//   return i*cols+j;
// }


MKL_LONG  MATRIX::test_arr()
{
  for(MKL_LONG  i=0; i<size; i++)
    {
      data[i] = i+1;
    } // i
  return 0;
}

MKL_LONG  MATRIX::test_arr_symm()
{
  for(MKL_LONG  i=0; i<rows; i++)
    {
      for(MKL_LONG  j=0; j<cols; j++)
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
  for(MKL_LONG  i=0; i<size; i++)
    {
      tmp_x = data[i];
      if(tmp_x < 0)
	result -= tmp_x;
      else
	result += tmp_x;
    }
  return result;
}

// inlined
// double MATRIX::norm()
// {
//   double result = 0.0;
//   for(MKL_LONG  i=0; i<size; i++)
//     {
//       result += pow(data[i],2.0);
//     }
//   return sqrt(result);
// }



double MATRIX::sum()
{
  double result = 0.0;
  for(MKL_LONG  i=0; i<size; i++)
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

MKL_LONG MATRIX::add(MATRIX& given_MAT)
{
  if (size != given_MAT.size)
    {
      std::cout << "ERR: MATRIX::add different size\n";
      return -1;
    }
  cblas_daxpy(size, 1.0, given_MAT.data, 1, data, 1);
  // for(MKL_LONG  k=0; k<size; k++)
  //   {
  //     data[k] += given_MAT.data[k];
  //   }
  return 0;
}

MKL_LONG MATRIX::sort()
{
  gsl_sort(data, 1, size);
  return 0;
}

MKL_LONG MATRIX::sort2(MATRIX& index)
{
  gsl_sort2(data, 1, index.data, 1, size);
  return 0;
}


// inlined
// MKL_LONG  make_unit_vector(MATRIX& given_vec)
// {
//   double norm = 0.;
//   for(MKL_LONG  i=0; i<given_vec.size; i++)
//     {
//       norm += pow(given_vec.data[i], 2.0);
//     }
//   norm = sqrt(norm);
//   for(MKL_LONG  i=0; i<given_vec.size; i++)
//     {
//       given_vec.data[i] /= norm;
//     }
//   return norm;
// }

// MATRIX return_unit_vector(MATRIX& given_vec)
// {
//   MATRIX re = given_vec;
//   double norm = 0.;
//   for(MKL_LONG  i=0; i<given_vec.size; i++)
//     {
//       norm += pow(given_vec.data[i], 2.0);
//     }
//   norm = sqrt(norm);
//   for(MKL_LONG  i=0; i<given_vec.size; i++)
//     {
//       re.data[i] /= norm;
//     }
//   return re;
// }



// inlined
// MKL_LONG matrix_mul(MATRIX& given_vec, double val)
// {
//   cblas_dscal(given_vec.size, val, given_vec.data, 1);
//   // for(MKL_LONG  k=0; k<given_vec.size; k++)
//   //   {
//   //     given_vec.data[k] *= val;
//   //   }
//   return 0;
// }



// Operator Overloading


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
//       for(MKL_LONG  i=0; i<A.size; i++)
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
//       for(MKL_LONG  i=0; i<A.size; i++)
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
//   for(MKL_LONG  i=0; i<C.size; i++)
//     {
//       C.data[i] *= a;
//     }
//   return C;
// }

// MATRIX Multiplification
// MATRIX operator*(const MATRIX &A, const MATRIX &B)
// {
//   MKL_LONG  C_rows = A.rows;
//   MKL_LONG  C_cols = B.cols;
//   MKL_LONG  cal_index = A.cols;
//   MATRIX C(C_rows, C_cols, 0.0);
//   if(A.cols == B.rows)
//     {
//       MKL_LONG  i, j, k;
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



// Matrix Operator
// MATRIX MATRIX::ROW(MKL_LONG  i)
// {
//   MATRIX C;
//   C.initial(1, cols);
//   for(MKL_LONG  j=0; j<cols; j++)
//     {
//       C.data[j] = data[i*cols+j];
//     }
//   return C;
// }

// MATRIX MATRIX::COL(MKL_LONG  j)
// {
//   MATRIX C;
//   C.initial(rows, 1);
//   for(MKL_LONG  i=0; i<rows; i++)
//     {
//       C.data[i] = data[i*cols+j];
//     }
//   return C;
// }


// MKL_LONG  MATRIX::ROW(const MATRIX &ROW_A, MKL_LONG  i)
// {
//   if(INITIALIZATION)
//     {
//       for(MKL_LONG  j=0; j<ROW_A.cols; j++)
// 	{
// 	  data[i*cols + j] = ROW_A.data[j];
// 	}
//     }
//   else
//     {
//       initial(1,ROW_A.cols);
//       for(MKL_LONG   j=0; j<ROW_A.cols; j++)
// 	{
// 	  data[j] = ROW_A.data[j];
// 	}
//     }
//   return 0;
// }

// MKL_LONG  MATRIX::COL(const MATRIX &COL_A, MKL_LONG  j)
// {
//   if(INITIALIZATION)
//     {
//       for(MKL_LONG  i=0; i<COL_A.rows; i++)
// 	{
// 	  data[i*cols + j] = COL_A.data[i];
// 	}
//     }
//   else
//     {
//       initial(COL_A.rows, 1);
//       for(MKL_LONG  i=0; i<COL_A.rows; i++)
// 	{
// 	  data[i] = COL_A.data[i];
// 	}
//     }
//   return 0;
// }

