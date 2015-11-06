
#ifndef LIB_RANDOM_H
#define LIB_RANDOM_H

extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_const_mksa.h>
}

#include <math.h>
#include <mkl.h>
#include <iostream>
#include "matrix_ed.h"

// note that the RANDOM package should not has the dependency with lib_traj.h

namespace RANDOM
{
  /* long dumbbell_connectivity(TRAJECTORY& TRAJ); */
  /* long random_position_dumbbell_generator(TRAJECTORY& TRAJ); */

  /* long fill_time(TRAJECTORY& TRAJ, double dt); */
  long random_vector_generator(MATRIX& R_VEC_TRANS);
  long single_random_vector_generator(MATRIX& given_vec);
  long single_random_vector_generator_variance(MATRIX& given_vec, double s_2);
  /* long random_position_generator_REF(TRAJECTORY& TRAJ, MATRIX& R_VEC_TRANS); */
  long single_unit_random_vector_generator(MATRIX& given_vec);
  long unit_random_vector_generator(MATRIX& R_VEC_TRANS);
  long unit_random_vector_generator_2D(MATRIX& R_VEC_TRANS);

  long return_LONG_INT_rand(long SUP);
  long return_LONG_INT_rand_boost(gsl_rng* r, long SUP);
  double return_double_rand_SUP1();
  double return_double_rand_SUP1_boost(gsl_rng* r);
  long get_LONG_ARR_rand_boost(gsl_rng* r, long SUP, long* given_long_arr, long N_arr);
  long get_DOUBLE_ARR_rand_boost(gsl_rng* r, double* given_double_arr, long N_arr);
  /* long gsl_rng_uniform_arr(gsl_rng* r, double*given_double_arr, long N_arr); */
  
}


#endif 
