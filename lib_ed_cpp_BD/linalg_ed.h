
/*
#################################################################
##   Private Linear Algebra Library ver 0.1                    ##
##   Made by "Gun Woo Park", Email : goraion@msn.com           ##
##   Lab. Polymer Dynamics                                     ##
##   Department of Polymer Science and Engineering             ##  
##   Kyungpook National University                             ##
##   Korea                                                     ##
#################################################################
*/

#ifndef LINALG_ED_H
#define LINALG_ED_H

#include <iostream>
#include "matrix_ed.h"
extern "C"
{
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
}

MATRIX solver_PLU(MATRIX A, MATRIX B);
MATRIX inverse(const MATRIX A);                    //@
MATRIX diagonalization_symm(const MATRIX A);       //@
MATRIX transpose(const MATRIX &A);                  //@
MATRIX dot_mat(MATRIX& A, MATRIX& B);
long dot_mat_ref(MATRIX& A, MATRIX& B, MATRIX& C);
double dot_vec(MATRIX& A, MATRIX& B);
MATRIX cross_3d(MATRIX& A, MATRIX& B);
long cross_3d_ref(MATRIX& A, MATRIX& B, MATRIX& C);
double det_3d(MATRIX& A);
#endif
