
#include "linalg_ed.h"
/*
MATRIX direct_cosine_3d(MATRIX& A, MATRIX& B)
{
  MATRIX C(3, 3, 0.0);
  MATRIX unit_A = (1.0/A.norm())*A;
  MATRIX unit_B = (1.0/B.norm())*B;
  
  // \lambda_{ij} = \cos(x'_i, x_j)


}
*/

MATRIX solver_PLU(MATRIX A, MATRIX b)
{
  MATRIX x(A.cols,1,0.0);
  gsl_matrix_view MAT = gsl_matrix_view_array(A.data, A.rows, A.cols);
  gsl_vector_view VEC = gsl_vector_view_array(b.data, b.size);
  gsl_permutation *p = gsl_permutation_alloc(A.rows);
  gsl_vector *sol = gsl_vector_alloc(A.rows);
  long s;
  gsl_linalg_LU_decomp(&MAT.matrix, p, &s);
  gsl_linalg_LU_solve(&MAT.matrix, p, &VEC.vector, sol);
  for(long i=0; i<b.size; i++)
    {
      x.data[i] = gsl_vector_get(sol,i);
    }
  gsl_permutation_free(p);
  gsl_vector_free(sol);

  return x;
}

MATRIX inverse(const MATRIX A) // Call By VALUE (GSL library acting on matrix A. That is reasion why I used call by value.
// By using PLU decomposition.
{
  MATRIX C;
  if(A.rows!=A.cols)
    {
      std::cout << "CANNOT CALCULATION INVERSE MATRIX when COLS!=ROWS\n";
    }
  else
    {
      //      MATRIX A = A; // gsl_matrix_inverse ACTING ON A. 
      gsl_matrix_view MAT = gsl_matrix_view_array(A.data, A.rows, A.cols);
      gsl_permutation *p = gsl_permutation_alloc(A.rows);
      gsl_matrix *inv_MAT = gsl_matrix_alloc(A.rows, A.cols);
      long s;
      gsl_linalg_LU_decomp(&MAT.matrix, p, &s);
      gsl_linalg_LU_invert(&MAT.matrix, p, inv_MAT);
      C.initial(A.rows, A.cols);
      for(long i=0; i<C.rows; i++)
	{
	  for(long j=0; j<C.cols; j++)
	    {
	      C(i,j) = gsl_matrix_get(inv_MAT, i, j);
	    } // j
	} // i

      gsl_matrix_free(inv_MAT);
      gsl_permutation_free(p);
    }
  return C;
}


MATRIX diagonalization_symm(const MATRIX A) // Call by Value
// It can apply symmetry matrix.
// PDinv(P)
{
  MATRIX C;
  if(A.rows != A.cols)
    {
      std::cout << "CANNOT CALCULATION INVERSE MATRIX when COLS!=ROWS\n";
    }
  else
    {
      C.initial(A.rows, A.cols);
      gsl_matrix_view MAT = gsl_matrix_view_array(A.data, A.rows, A.cols);
      gsl_vector *eval_MAT = gsl_vector_alloc(A.rows);
      gsl_matrix *evec_MAT = gsl_matrix_alloc(A.rows, A.cols);
      gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(A.rows);
      gsl_eigen_symmv(&MAT.matrix, eval_MAT, evec_MAT, w);
      gsl_eigen_symmv_free(w);
      C.DIAGONALIZATION_INITIAL();
      for(long i=0; i<C.rows; i++)
	{
	  C.eigen_value[i] = gsl_vector_get(eval_MAT, i);
	  for(long j=0; j<C.cols; j++)
	    {
	      C(i,j) = gsl_matrix_get(evec_MAT, i, j);
	    }
	}
    }
  return C;
}


MATRIX dot_mat(MATRIX& A, MATRIX& B)
{
  MATRIX C(A.rows, B.cols, 0.0);
  if(A.cols != B.rows)
    {
      std::cout << "Dimension does NOT match in dot_mat function\n";
    }
  else
    {
      for(long i=0; i<A.rows; i++)
	{
	  for(long j=0; j<B.cols; j++)
	    {
	      for(long k=0; k<A.cols; k++)
		{
		  C(i,j) += A(i,k)*B(k,j);
		}
	    }
	}
    }
  return C;
}

long dot_mat_ref(MATRIX& A, MATRIX& B, MATRIX& C)
{
  C.initial(A.rows, B.cols, 0.0);
  if(A.cols != B.rows)
    {
      std::cout << "Dimension does NOT match in dot_mat function\n";
    }
  else
    {
      for(long i=0; i<A.rows; i++)
	{
	  for(long j=0; j<B.cols; j++)
	    {
	      for(long k=0; k<A.cols; k++)
		{
		  C(i,j) += A(i,k)*B(k,j);
		}
	    }
	}
    }
  return 0;
}


double dot_vec(MATRIX& A, MATRIX& B)
{
  double result = 0.0;
  if(A.size != B.size)
    {
      std::cout << "Dimension does NOT match in dot_vec function\n";
    }
  else
    {
      for(long i=0; i<A.size; i++)
	{
	  result += A(i)*B(i);
	}
    }
  return result;
}

MATRIX cross_3d(MATRIX& A, MATRIX& B)
{
  MATRIX C(3, 1, 0.0);
  C(0) = A(1)*B(2) - A(2)*B(1);
  C(1) = A(2)*B(0) - A(0)*B(2);
  C(2) = A(0)*B(1) - A(1)*B(0);
  return C;
}

long cross_3d_ref(MATRIX& A, MATRIX& B, MATRIX& C)
{
  C.initial(3, 1, 0.0);
  C(0) = A(1)*B(2) - A(2)*B(1);
  C(1) = A(2)*B(0) - A(0)*B(2);
  C(2) = A(0)*B(1) - A(1)*B(0);
  return 0;
}

double det_3d(MATRIX& A)
{
  double result = A(0,0)*(A(1,1)*A(2,2)-A(1,2)*A(2,1)) - A(0,1)*(A(1,0)*A(2,2) - A(1,2)*A(2,0)) + A(0,2)*(A(1,0)*A(2,1) - A(1,1)*A(2,0));
  return result;
}
