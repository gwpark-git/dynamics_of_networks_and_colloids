
#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <iostream>
#include "trajectory.h"
/* #include "read_file_condition.h" */
#include "file_IO.h"
#include "cell_list.h"

class RDIST : public CLIST
{
 public:

  /* double** Rmin; */
  MATRIX** Rvec; // Rvec[i][j] will be the relative vector between i- and j-th particles
  MATRIX* Rsca; // Rsca[i](j) will be the relative distance between i- and j-th particles. The form is differ from the original.
  
  RDIST()
    {
      std::cout << "There is no empty constructor for RDIST class\n" << std::endl;
    }
  RDIST(COND& given_condition);
  ~RDIST()
    {
      mkl_free(Rsca);
      for(MKL_LONG i=0; i<Np; i++)
        {
          mkl_free(Rvec[i]);
        }
      mkl_free(Rvec);
    }
};

namespace GEOMETRY
{
  double get_minimum_distance(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j, MATRIX& given_vec);
  double get_minimum_distance_cell_list(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j, MATRIX& given_vec, MKL_LONG* beyond_box_check);
  double get_minimum_distance_for_particle(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_particle, MATRIX& R_minimum_boost_particle, MATRIX** R_minimum_vec_boost);  

  double return_minimum_distance(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG index_j);
  MKL_LONG minimum_image_convention(TRAJECTORY& TRAJ, MKL_LONG target_t);
  MKL_LONG get_minimum_distance_pos_vector(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec);
  MKL_LONG get_minimum_distance_pos_vector(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec);
  MKL_LONG get_minimum_distance_rel_vector(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec);
  MKL_LONG get_minimum_distance_rel_vector_cell_list(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG given_index, MKL_LONG target_index, MATRIX& given_vec, MKL_LONG* beyond_box_check);

  double get_simple_distance(TRAJECTORY& TRAJ, MKL_LONG index_t, MKL_LONG index_i, MKL_LONG i);
}

namespace UTIL_ARR
{
  double get_minimum_image_k_from_x(double x, double k, double dimension);
  MKL_LONG get_index_minimum_abs(double *k, MKL_LONG N);
}


#endif
