
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
  /* MKL_LONG dumbbell_connectivity(TRAJECTORY& TRAJ); */
  /* MKL_LONG random_position_dumbbell_generator(TRAJECTORY& TRAJ); */

  /* MKL_LONG fill_time(TRAJECTORY& TRAJ, double dt); */
  MKL_LONG random_vector_generator(MATRIX& R_VEC_TRANS);
  MKL_LONG single_random_vector_generator(MATRIX& given_vec);
  MKL_LONG single_random_vector_generator_variance(MATRIX& given_vec, double s_2);
  /* MKL_LONG random_position_generator_REF(TRAJECTORY& TRAJ, MATRIX& R_VEC_TRANS); */
  MKL_LONG single_unit_random_vector_generator(MATRIX& given_vec);
  MKL_LONG unit_random_vector_generator(MATRIX& R_VEC_TRANS);
  MKL_LONG unit_random_vector_generator_2D(MATRIX& R_VEC_TRANS);

  MKL_LONG return_LONG_INT_rand(MKL_LONG SUP);
  MKL_LONG return_LONG_INT_rand_boost(gsl_rng* r, MKL_LONG SUP);

  double return_double_rand_SUP1();
  double return_double_rand_SUP1_boost(gsl_rng* r);

  
}


#endif 
