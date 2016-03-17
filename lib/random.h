
#ifndef RANDOM_H
#define RANDOM_H

extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_const_mksa.h>
}

#include <math.h>
#include <mkl.h>
#include <iostream>
#include "matrix.h"

// note that the RANDOM package should not has the dependency with lib_traj.h

namespace RANDOM
{
  MKL_LONG random_vector_generator(MATRIX& R_VEC_TRANS);
  MKL_LONG single_random_vector_generator(MATRIX& given_vec);
  MKL_LONG single_random_vector_generator_variance(MATRIX& given_vec, double s_2);
  MKL_LONG single_random_vector_generator_boost(MATRIX& given_vec, gsl_rng* r_boost);
  MKL_LONG single_random_vector_generator_variance_boost(MATRIX& given_vec, double s_2, gsl_rng* r_boost);
  MKL_LONG single_unit_random_vector_generator(MATRIX& given_vec);
  MKL_LONG unit_random_vector_generator(MATRIX& R_VEC_TRANS);
  MKL_LONG unit_random_vector_generator_2D(MATRIX& R_VEC_TRANS);

  MKL_LONG return_LONG_INT_rand(MKL_LONG SUP);
  MKL_LONG return_LONG_INT_rand_boost(gsl_rng* r, MKL_LONG SUP);
  double return_double_rand_SUP1();
  double return_double_rand_SUP1_boost(gsl_rng* r);
  MKL_LONG get_LONG_ARR_rand_boost(gsl_rng* r, MKL_LONG SUP, MKL_LONG* given_long_arr, MKL_LONG N_arr);
  MKL_LONG get_DOUBLE_ARR_rand_boost(gsl_rng* r, double* given_double_arr, MKL_LONG N_arr);
  
}


#endif 
