
#include "matrix.h"
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <math.h>

MKL_LONG
MATRIX::
read_exist_data(const char* fn_given_array)
{
  std::ifstream GIVEN_FILE;
  GIVEN_FILE.open(fn_given_array);
  MKL_LONG cnt = 0;
  std::string line;
  MKL_LONG cols = 0;
  while(getline(GIVEN_FILE, line))
    {
      if(cnt ==0)
        {
          std::stringstream is(line);
          double tmp;
          while(is >> tmp)
            {cols ++;}
        }
      cnt ++;
    }
  MKL_LONG rows = cnt;
  GIVEN_FILE.clear();
  GIVEN_FILE.seekg(0);
  
  initial(rows, cols, 0.);
  // printf("rows = %ld, cols = %ld\n", rows, cols);
  
  for(MKL_LONG i=0; i<rows; i++)
    {
      for(MKL_LONG j=0; j<cols; j++)
        {
          GIVEN_FILE >> (*this)(i, j);
          // printf("%ld\t", (*this)(i,j));
        }
      // printf("\n");
    }
  GIVEN_FILE.close();
  return 0;
}

MKL_LONG
nonzero
(const MATRIX &A)
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
MATRIX::
MATRIX()
{
  initial();  
  //std::cout << "Undefined dimension\n";
}

MATRIX::
MATRIX
(MKL_LONG  N_r, MKL_LONG  N_c)
{
  initial(N_r, N_c);
}


MATRIX::
MATRIX
(MKL_LONG  N_r, MKL_LONG  N_c, double x) // initilization with value x
{
  initial(N_r, N_c);
  set_value(x);
}

// Copy-Constructor
MKL_LONG
MATRIX::
copy_obj
(const MATRIX& Mat)
{
  initial(Mat.rows, Mat.cols);
  for(MKL_LONG i=0; i<size; i++)
    {
      data[i] = Mat.data[i];
    }
  return 0;
}

MATRIX::
MATRIX
(const MATRIX& Mat)
{
  copy_obj(Mat);
}

// Destructor
MATRIX::
~MATRIX()
{
  if (INITIALIZATION)
    {
      data_delete();
    }
  if (DIAGONALIZATION)
    {
      // mkl_free(eigen_value);
      delete[] eigen_value;
    }
}


// Public Member Function

MKL_LONG
MATRIX::
print_eigen()
{
  print_eigen(rows);
  return 0;
}

MKL_LONG
MATRIX::
print_eigen
(MKL_LONG  n_ele)
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

MKL_LONG
MATRIX::
print()
{
  print(rows, cols);
  return 0;
}

MKL_LONG
MATRIX::
print
(MKL_LONG  n_row, MKL_LONG  n_col)
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

MKL_LONG
MATRIX::
fprint_skip
(std::ofstream& file, MKL_LONG  N_skip)
{
  MKL_LONG  cnt = 0;
  for(MKL_LONG  i=0; i<rows; i+=N_skip)
    {
      for(MKL_LONG  j=0; j<cols; j++)
        {
          file << std::scientific << data[index(j, i)] << '\t';
        }
      file << std::endl;
    } // i
  return 0;
}

MKL_LONG
MATRIX::
fprint
(std::ofstream& file)
{
  fprint_skip(file, 1);
  return 0;
}

MKL_LONG
MATRIX::
fprint_skip_transpose
(std::ofstream& file, MKL_LONG  N_skip)
{
  MKL_LONG  cnt = 0;
  for(MKL_LONG  j=0; j<cols; j++)
    {
      for(MKL_LONG  i=0; i<rows; i+=N_skip)
        {
          file << std::scientific << data[index(i, j)] << '\t';
        }
      file << std::endl;
    } // i
  return 0;
}

MKL_LONG
MATRIX::
fprint_transpose
(std::ofstream& file)
{
  fprint_skip_transpose(file, 1);
  return 0;
}

MKL_LONG
MATRIX::
fprint_LONG_skip
(std::ofstream& file, MKL_LONG  N_skip)
{
  MKL_LONG  cnt = 0;
  for(MKL_LONG  i=0; i<rows; i+=N_skip)
    {
      for(MKL_LONG  j=0; j<cols; j++)
        {
          file << (MKL_LONG)data[index(i, j)] << '\t';
        }
      file << std::endl;
    } // i
  return 0;
}

MKL_LONG
MATRIX::
fprint_LONG
(std::ofstream& file)
{
  fprint_LONG_skip(file, 1);
  return 0;
}

MKL_LONG
MATRIX::
fprint_LONG_skip_transpose_LIMROWS
(std::ofstream& file, MKL_LONG N_skip, MKL_LONG N_lim_rows)
{
  MKL_LONG  cnt = 0;
  for(MKL_LONG  j=0; j<cols; j++)
    {
      for(MKL_LONG  i=0; i<N_lim_rows; i+=N_skip)
        {
          file << (MKL_LONG)data[index(i, j)] << '\t';
        }
      file << std::endl;
    } // i
  return 0;
}

MKL_LONG
MATRIX::
fprint_LONG_skip_transpose
(std::ofstream& file, MKL_LONG  N_skip)
{
  MATRIX::fprint_LONG_skip_transpose_LIMROWS(file, N_skip, rows); // count all rows
  return 0;
}

MKL_LONG
MATRIX::
fprint_LONG_transpose
(std::ofstream& file)
{
  fprint_LONG_skip_transpose(file, 1);
  return 0;
}


MKL_LONG
MATRIX::
fprint_row
(std::ofstream& file, MKL_LONG  given_row)
{
  MKL_LONG  cnt = 0;
  MKL_LONG  i=given_row;
  for(MKL_LONG  j=0; j<cols; j++)
    {
      file << std::scientific << data[index(i, j)] << '\t';
    }
  file << std::endl;
  return 0;
}

MKL_LONG
MATRIX::
fprint_out_skip
(std::ofstream& file, MKL_LONG  N_skip)
{
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
  return 0;

}

MKL_LONG
MATRIX::
fprint_out
(std::ofstream& file)
{
  fprint_out_skip(file, 1);
  return 0;
}

MKL_LONG
MATRIX::
p_self()
{
  std::cout << this;
  return 0;
}

MKL_LONG
MATRIX::
nonzero()
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

MKL_LONG
MATRIX::
ABS_cond
(double x)
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
MKL_LONG
MATRIX::
initial()
{
  INITIALIZATION = FALSE;
  DIAGONALIZATION = FALSE;
  rows = 0;
  cols = 0;
  size = 0;
  data = NULL;
  return 0;
}

MKL_LONG
MATRIX::
initial
(MKL_LONG  N_r, MKL_LONG  N_c)
{
  // std::cout << "CONSTRUCTOR " << this << std::endl;
  // if (INITIALIZATION == TRUE)
  //   {
  //     rows = N_r;
  //     cols = N_c;
  //     size = rows*cols;
  //     delete[] data;
  //     data = new double [size];
  //     for(MKL_LONG i=0; i<size; i++)
  //       data[i] = 0.;
  //     return 0;
  //   }
  // else
  //   {
      // this is of importance to prevent memory leak
  // if (INITIALIZATION == TRUE)
  //   data_delete();
  INITIALIZATION = TRUE;
  DIAGONALIZATION = FALSE;
  rows = N_r;
  cols = N_c;
  size = rows*cols;
  // note that the dynamic allocation should be defined callocation
  // otherwise, it should be set with the valid initial variable
  // data = (double*) mkl_calloc(size, sizeof(double), BIT);
  data = new double [size];
  for(MKL_LONG i=0; i<size; i++)
    data[i] = 0.;
    //   // data = new double [size];
    // }
  return 0;
}

MKL_LONG
MATRIX::
initial
(MKL_LONG  N_r, MKL_LONG  N_c, double x)
{
  initial(N_r, N_c);
  for(MKL_LONG  i=0; i<size; i++)
    {
      data[i] = x;
    }
  return 0;
}


MKL_LONG
MATRIX::
data_delete()
{
  // delete[] data;
  // mkl_free(data);
  delete[] data;
  INITIALIZATION = FALSE;
  data = NULL;
  return 0;
}

MKL_LONG
MATRIX::
undefined_error()
{
  std::cout << "MATRIX is NOT constructed\n" << std::endl;
  return -1;
}


MKL_LONG
MATRIX::
test_arr()
{
  for(MKL_LONG  i=0; i<size; i++)
    {
      data[i] = i+1;
    } // i
  return 0;
}

MKL_LONG
MATRIX::
test_arr_symm()
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

double
MATRIX::
ABS_sum()
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

double
MATRIX::
sum()
{
  double result = 0.0;
  for(MKL_LONG  i=0; i<size; i++)
    {
      result += data[i];
    }
  return result;
}

double
MATRIX::
average()
{
  double result = sum();
  return result/(double)size;
}

MKL_LONG
MATRIX::
add
(MATRIX& given_MAT)
{
  if (size != given_MAT.size)
    {
      std::cout << "ERR: MATRIX::add different size\n";
      return -1;
    }
  cblas_daxpy(size, 1.0, given_MAT.data, 1, data, 1);
  return 0;
}

MKL_LONG
MATRIX::
sort()
{
  gsl_sort(data, 1, size);
  return 0;
}


// MKL_LONG MATRIX::sort2(MATRIX& index)
// {

// MKL_LONG MATRIX::sort2(MATRIX& index)
// {
//   gsl_sort2(data, 1, index.data, 1, size);
//   return 0;
// }

// MKL_LONG MATRIX::sort2(size_t* index)
// {
//   // it is of importance that the gsl_sort2 is not supported by the SCOPE_GRID server.
//   // hence, we are applying new scheme
//   // gsl_sort2(data, 1, index.data, 1, size);
  
//   gsl_sort_index(index, data, 1, size);
//   gsl_sort(data, 1, size);
//   return 0;
// }

MKL_LONG
MATRIX::
sort2
(MATRIX& index)
{

  gsl_sort2(data, 1, index.data, 1, size);
  return 0;
}
