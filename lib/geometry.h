
#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <iostream>
#include "trajectory.h"

namespace GEOMETRY
{
  double get_minimum_distance(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j, MATRIX& given_vec);
  double get_minimum_distance_for_particle(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_particle, MATRIX& R_minimum_boost_particle, MATRIX** R_minimum_vec_boost);  

  double return_minimum_distance(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j);
  MKL_LONG minimum_image_convention(TRAJECTORY& TRAJ, MKL_LONG target_t);
  MKL_LONG get_minimum_distance_pos_vector(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec);
  MKL_LONG get_minimum_distance_pos_vector(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec);
  MKL_LONG get_minimum_distance_rel_vector(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec);

  double get_simple_distance(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG i);
}

namespace UTIL_ARR
{
  double get_minimum_image_k_from_x(double x, double k, double dimension);
  MKL_LONG get_index_minimum_abs(double *k, MKL_LONG N);
}


#endif
