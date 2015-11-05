
#ifndef LIB_GEOMETRY_H
#define LIB_GEOMETRY_H

#include <iostream>
#include "lib_traj.h"

namespace GEOMETRY
{
  double get_minimum_distance(TRAJECTORY& TRAJ, long index_t, long index_i, long index_j, MATRIX& given_vec);
  double return_minimum_distance(TRAJECTORY& TRAJ, long index_t, long index_i, long index_j);
  long minimum_image_convention(TRAJECTORY& TRAJ, long target_t);
  long get_minimum_distance_pos_vector(TRAJECTORY& TRAJ, long index_t, long given_index, long target_index, MATRIX& given_vec);
  /* MATRIX return_minimum_distance_pos_vector(TRAJECTORY& TRAJ, long index_t, long given_index, long target_index); */
  long get_minimum_distance_pos_vector(TRAJECTORY& TRAJ, long index_t, long given_index, long target_index, MATRIX& given_vec);
  long get_minimum_distance_rel_vector(TRAJECTORY& TRAJ, long index_t, long given_index, long target_index, MATRIX& given_vec);

  /* MATRIX return_minimum_distance_rel_vector(TRAJECTORY& TRAJ, long index_t, long given_index, long target_index); */
  double get_simple_distance(TRAJECTORY& TRAJ, long index_t, long index_i, long i);
}

namespace UTIL_ARR
{
  double get_minimum_image_k_from_x(double x, double k, double dimension);
  long get_index_minimum_abs(double *k, long N);
}


#endif
